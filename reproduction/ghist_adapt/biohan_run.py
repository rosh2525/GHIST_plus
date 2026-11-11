"""v9: BioHAN, a biology-aware hurdle adapter network trained on core cells (protocol config_v9.json).

Per cell, inputs (all from GHIST+ outputs and coordinates; no RNA outside cores):
  standardised [S, logit q] for all genes                       (as HAC)
  niche: mean S over same-predicted-type neighbours within 30 um (cell state follows its microenvironment)
  niche: predicted cell-type composition within 50 um           (tissue context, e.g. immune cells near tumour)
  predicted cell-type probabilities (9)
Architecture: MLP trunk -> 64-unit gene-programme bottleneck -> two heads over genes:
  occurrence  q_hat = sigmoid(logit q + d_occ)                  "is the gene on?"
  magnitude   m_hat = 0.01 + softplus(c_g + w_g * relu(mu) + d_mag)   "how much, if on?"
  output      q_hat * m_hat  (strictly positive: no zero ties, the failure mode of the v2 MLP)
Loss: BCE(q_hat, 1[y>0]) + MSE(m_hat, y) on expressing entries; rare cell types up-weighted by 1/sqrt(frequency).
Heads start at zero, so training starts from GHIST+'s own q and a per-gene magnitude fit.
Epoch count chosen by holding out one core group; final model refit on all cores for that many epochs.
Ablation BioHAN-noniche: same network without the two niche inputs.

Usage:  python biohan_run.py breast5k [--smoke]
"""
from __future__ import annotations

import argparse
import gc
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from scipy import sparse
from scipy.spatial import cKDTree

import adapt
from anchors import evaluation_mask, place_cores, prediction_pcs
from hurdle import load_hurdle
from metrics import pcc_columns
from p1_run import CONFIGS, _hac_predict, evaluate, hac_fit, slide_stats
from run_v3 import CFG, gene_groups
from slides import load

HERE = Path(__file__).resolve().parent
DEV = torch.device("cuda" if torch.cuda.is_available() else "cpu")
HID, K, EPOCHS_MAX, PATIENCE, BATCH, LR, WD, DROP = 512, 64, 60, 8, 512, 1e-3, 1e-4, 0.1


def neighbour_mean(xy, r, ct=None):
    pairs = cKDTree(xy).query_pairs(r, output_type="ndarray")
    if ct is not None:
        pairs = pairs[ct[pairs[:, 0]] == ct[pairs[:, 1]]]
    n = len(xy)
    i = np.concatenate([pairs[:, 0], pairs[:, 1]]); j = np.concatenate([pairs[:, 1], pairs[:, 0]])
    W = sparse.csr_matrix((np.ones(len(i), np.float32), (i, j)), shape=(n, n))
    deg = np.asarray(W.sum(1)).ravel()
    return sparse.diags(np.where(deg > 0, 1.0 / np.maximum(deg, 1), 0).astype(np.float32)) @ W


class BioHAN(nn.Module):
    def __init__(self, d_in, G, c0):
        super().__init__()
        self.trunk = nn.Sequential(nn.Linear(d_in, HID), nn.GELU(), nn.Dropout(DROP),
                                   nn.Linear(HID, HID), nn.GELU(), nn.Dropout(DROP), nn.Linear(HID, K), nn.GELU())
        self.occ = nn.Linear(K, G); self.mag = nn.Linear(K, G)
        for h in (self.occ, self.mag):
            nn.init.zeros_(h.weight); nn.init.zeros_(h.bias)
        self.c = nn.Parameter(torch.as_tensor(c0)); self.w = nn.Parameter(torch.zeros(G))

    def forward(self, x, lq, r):
        z = self.trunk(x)
        qh = torch.sigmoid(lq + self.occ(z))
        mh = 0.01 + nn.functional.softplus(self.c + self.w * r + self.mag(z))
        return qh, mh


def inv_softplus(v):
    v = np.maximum(v, 0.02)
    return np.log(np.expm1(v)).astype(np.float32)


def train(X, LQ, R, Y, wts, epochs=None, Xv=None, LQv=None, Rv=None, Yv=None, seed=0):
    torch.manual_seed(seed); np.random.seed(seed)
    pos = Y > 0
    c0 = inv_softplus(np.where(pos.sum(0) > 0, (Y * pos).sum(0) / np.maximum(pos.sum(0), 1), 0.5))
    net = BioHAN(X.shape[1], Y.shape[1], c0).to(DEV)
    opt = torch.optim.AdamW(net.parameters(), lr=LR, weight_decay=WD)
    T = lambda a: torch.as_tensor(a, device=DEV)
    Xt, LQt, Rt, Yt, Wt = T(X), T(LQ), T(R), T(Y), T(wts)
    probs = Wt / Wt.sum()
    n = len(X); steps = max(1, n // BATCH)
    best, best_ep, bad = -np.inf, 0, 0
    n_ep = epochs or EPOCHS_MAX
    for ep in range(1, n_ep + 1):
        net.train()
        for _ in range(steps):
            idx = torch.multinomial(probs, BATCH, replacement=True)
            qh, mh = net(Xt[idx], LQt[idx], Rt[idx])
            y = Yt[idx]; on = (y > 0).float()
            loss = nn.functional.binary_cross_entropy(qh.clamp(1e-5, 1 - 1e-5), on) \
                + ((mh - y) ** 2 * on).sum() / on.sum().clamp(min=1)
            opt.zero_grad(); loss.backward(); opt.step()
        if epochs is None:
            sc = float(np.nanmean(pcc_columns(Yv, predict(net, Xv, LQv, Rv))))
            if sc > best + 1e-4:
                best, best_ep, bad = sc, ep, 0
            else:
                bad += 1
                if bad >= PATIENCE:
                    break
    return net, best_ep, best


@torch.no_grad()
def predict(net, X, LQ, R, chunk=50_000):
    net.eval(); out = np.empty(R.shape, np.float32)
    for a in range(0, len(X), chunk):
        T = lambda v: torch.as_tensor(v[a:a + chunk], device=DEV)
        qh, mh = net(T(X), T(LQ), T(R))
        out[a:a + chunk] = (qh * mh).cpu().numpy()
    return out


def fit_biohan(X, LQ, R, Y, ct, groups, seed):
    freq = np.bincount(ct, minlength=ct.max() + 1).astype(np.float64)
    wts = (1.0 / np.sqrt(np.maximum(freq[ct], 1))).astype(np.float32)
    gs = np.unique(groups)
    val = groups == gs[seed % len(gs)]
    _, ep, _ = train(X[~val], LQ[~val], R[~val], Y[~val], wts[~val], None, X[val], LQ[val], R[val], Y[val], seed)
    net, _, _ = train(X, LQ, R, Y, wts, epochs=max(ep, 1), seed=seed)
    return net, max(ep, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slides", nargs="+")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    a = CFG["anchors"]
    print("device:", DEV, torch.cuda.get_device_name(0) if DEV.type == "cuda" else "", flush=True)
    for slide in args.slides:
        t0 = time.time()
        d = load(slide)
        h = load_hurdle(slide, d["ids"])
        q = h.pop("q"); R = np.maximum(h.pop("mu"), 0).astype(np.float32)
        S = (q * R).astype(np.float32)
        ctp = h["ctp"].astype(np.float32); ct = ctp.argmax(1)
        Y, P, xy = d["y"], d["pred"], d["xy_um"]
        _, groups = gene_groups(d, slide)
        lq = np.log(np.clip(q, 1e-4, 1 - 1e-4) / (1 - np.clip(q, 1e-4, 1 - 1e-4))).astype(np.float32)
        del q
        featS, featL = slide_stats(S), slide_stats(lq)
        NS = (neighbour_mean(xy, 30.0, ct) @ S).astype(np.float32)
        NC = (neighbour_mean(xy, 50.0) @ ctp).astype(np.float32)
        featN, featC = slide_stats(NS), slide_stats(NC)
        z = lambda A, f: ((A - f[0]) / f[1]).astype(np.float32)
        base = lambda idx: [z(S[idx], featS), z(lq[idx], featL), ctp[idx]]
        X_full = lambda idx: np.hstack(base(idx) + [z(NS[idx], featN), z(NC[idx], featC)])
        X_noni = lambda idx: np.hstack(base(idx))
        Xh = lambda idx: np.hstack([z(S[idx], featS), z(lq[idx], featL)])
        print(f"[{slide}] {len(Y)} cells x {Y.shape[1]} genes; features built ({time.time()-t0:.0f}s)", flush=True)
        tree = cKDTree(xy); pcs = prediction_pcs(P)
        out = HERE / ("outputs_v9_smoke" if args.smoke else "outputs_v9") / slide
        out.mkdir(parents=True, exist_ok=True)
        rows, genes_tab = [], []
        configs = [(1.0, 4, "diverse", 0)] if args.smoke else [(ar, n, pl, s) for ar, n, pl in CONFIGS for s in a["repeats"]]
        for ar, n, pl, s in configs:
            t1 = time.time()
            centres, rad, aidx, ag = place_cores(xy, ar, n, pl, s, P=P, pcs=pcs, tree=tree)
            eidx = np.flatnonzero(evaluation_mask(xy, centres, rad, a["evaluation_buffer_um"]))
            Ya = Y[aidx]
            rrs = adapt.fit_rrc_authors(S[aidx], Ya, *featS)
            hac = hac_fit(Xh(aidx), (1 / (1 + np.exp(-lq[aidx]))).astype(np.float32), R[aidx], Ya, ag)
            net_f, ep_f = fit_biohan(X_full(aidx), lq[aidx], R[aidx], Ya, ct[aidx], ag, s)
            net_n, ep_n = fit_biohan(X_noni(aidx), lq[aidx], R[aidx], Ya, ct[aidx], ag, s)
            fit_s = time.time() - t1
            qe = 1 / (1 + np.exp(-lq[eidx]))
            preds = {"GHIST+": P[eidx], "GHIST+-Soft": S[eidx], "RRC-soft": adapt.predict_rrc_authors(rrs, S[eidx]),
                     "HAC": _hac_predict(hac, Xh(eidx), qe, R[eidx]),
                     "BioHAN": predict(net_f, X_full(eidx), lq[eidx], R[eidx]),
                     "BioHAN-noniche": predict(net_n, X_noni(eidx), lq[eidx], R[eidx])}
            res, pg = evaluate(Y[eidx], preds, groups, xy[eidx], s, refs=("RRC-soft", "BioHAN-noniche"))
            key = f"a{ar}_n{n}_{pl}_s{s}"
            for meth, rw in res.items():
                rw.update(slide=slide, config=key, seed=s, method=meth, n_anchor_cells=int(len(aidx)),
                          n_eval_cells=int(len(eidx)), epochs_full=ep_f, epochs_noniche=ep_n, fit_seconds=round(fit_s, 1))
                rows.append(rw)
                genes_tab.append(pd.DataFrame({"config": key, "method": meth, "gene_idx": np.arange(len(pg[meth])), "pcc": pg[meth]}))
            f = lambda m: res[m]["median_pcc_all"]
            print(f"[{slide}] {key}: anchors={len(aidx)} ep={ep_f}/{ep_n} fit={fit_s:.0f}s | medPCC Soft {f('GHIST+-Soft'):.4f} "
                  f"RRC-soft {f('RRC-soft'):.4f} HAC {f('HAC'):.4f} BioHAN {f('BioHAN'):.4f} noniche {f('BioHAN-noniche'):.4f} "
                  f"(CI vs RRC-soft [{res['BioHAN']['ci_vs_RRC-soft_all_lo']:+.4f},{res['BioHAN']['ci_vs_RRC-soft_all_hi']:+.4f}]) | "
                  f"SSIM {res['RRC-soft']['median_ssim_svg50']:.3f}/{res['BioHAN']['median_ssim_svg50']:.3f} "
                  f"CMD {res['RRC-soft']['cmd_all']:.3f}/{res['BioHAN']['cmd_all']:.3f} ({time.time()-t0:.0f}s)", flush=True)
            pd.DataFrame(rows).to_csv(out / "runs.csv", index=False, float_format="%.6g")
            pd.concat(genes_tab).to_csv(out / "per_gene.csv.gz", index=False, float_format="%.5g")
            del preds, res, pg, net_f, net_n; gc.collect(); torch.cuda.empty_cache()
        print(f"[{slide}] done ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
