"""v8: biology-driven extensions on top of RRC-soft (protocol config_v8.json, frozen before any run).

Fix 2  Niche: tissue niches are spatially coherent within a cell type. Smooth the soft output S over neighbours of
       the SAME predicted cell type (GHIST+'s H&E cell-type head): S_n = (1-lam) S + lam * mean_{same-type nbrs < r} S.
       lam, r chosen on core cells only (no RNA outside cores). Niche-RRC = authors' RRC applied to S_n.
Fix 3  Programme: genes act in co-regulated programmes. Project the RRC-soft prediction onto the top-k principal
       components of the measured core RNA (standardised), k chosen by leave-one-core-out on core cells.
Both   Niche-RRC followed by the programme projection.

Usage:  python bio_run.py breast5k [atera breast5k_pan] [--smoke]
Writes  outputs_v8/<slide>/{runs.csv, per_gene.csv.gz}
"""
from __future__ import annotations

import argparse
import gc
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

import adapt
from anchors import evaluation_mask, place_cores, prediction_pcs
from hurdle import load_hurdle
from metrics import pcc_columns
from one_shot_rrc import fit_rrc as rrc_fit
from p1_run import CONFIGS, evaluate, slide_stats
from run_v3 import CFG, gene_groups
from slides import load

HERE = Path(__file__).resolve().parent
LAMS = [0.0, 0.25, 0.5, 0.75]
RADII = [15.0, 30.0]
KS = [8, 16, 32, 64]


def same_type_matrix(xy, ct, r):
    tree = cKDTree(xy)
    pairs = tree.query_pairs(r, output_type="ndarray")
    pairs = pairs[ct[pairs[:, 0]] == ct[pairs[:, 1]]]
    n = len(xy)
    i = np.concatenate([pairs[:, 0], pairs[:, 1]])
    j = np.concatenate([pairs[:, 1], pairs[:, 0]])
    W = sparse.csr_matrix((np.ones(len(i), np.float32), (i, j)), shape=(n, n))
    deg = np.asarray(W.sum(1)).ravel()
    return sparse.diags(np.where(deg > 0, 1.0 / np.maximum(deg, 1), 0).astype(np.float32)) @ W, deg > 0


def niche(S, WS, has, lam):
    out = S.copy()
    out[has] = (1 - lam) * S[has] + lam * WS[has]
    return out


def choose_niche(S, WS_by_r, has_by_r, Y, aidx):
    """lam, r maximising mean per-gene PCC on core cells (zero-shot; uses core RNA only)."""
    best = (-np.inf, 0.0, RADII[0])
    for r in RADII:
        for lam in LAMS:
            Sn = niche(S[aidx], WS_by_r[r][aidx], has_by_r[r][aidx], lam)
            sc = float(np.nanmean(pcc_columns(Y[aidx], Sn)))
            if sc > best[0]:
                best = (sc, lam, r)
    return best[1], best[2]


def rrc_oof(B, Y, groups, feat, r):
    x = ((B - feat[0]) / feat[1]).astype(np.float32)
    oof = np.zeros_like(Y, dtype=np.float32)
    for f in np.unique(groups):
        te, tr = groups == f, groups != f
        xm, ym, coef = rrc_fit(x[tr], (Y[tr] - B[tr]).astype(np.float32), r["alpha"])
        oof[te] = np.clip(B[te] + r["blend"] * ((x[te] - xm) @ coef + ym), 0, None)
    return oof


def prog_fit(Y, k):
    m = Y.mean(0); s = Y.std(0); s[s == 0] = 1.0
    _, _, Vt = np.linalg.svd((Y - m) / s, full_matrices=False)
    return m.astype(np.float32), s.astype(np.float32), Vt[:k].T.astype(np.float32)


def prog_apply(Q, p, chunk=200_000):
    m, s, V = p
    out = np.empty_like(Q)
    for a in range(0, len(Q), chunk):
        Z = (Q[a:a + chunk] - m) / s
        out[a:a + chunk] = np.clip((Z @ V) @ V.T * s + m, 0, None)
    return out


def choose_k(oof, Y, groups):
    """k maximising leave-one-core-out mean PCC of the projected out-of-fold RRC predictions."""
    best = (-np.inf, None)
    for k in KS + [None]:
        if k is None:
            sc = float(np.nanmean(pcc_columns(Y, oof)))
        else:
            proj = np.zeros_like(oof)
            for f in np.unique(groups):
                te = groups == f
                if (~te).sum() <= k:
                    continue
                proj[te] = prog_apply(oof[te], prog_fit(Y[~te], k))
            sc = float(np.nanmean(pcc_columns(Y, proj)))
        if sc > best[0]:
            best = (sc, k)
    return best[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slides", nargs="+")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    a = CFG["anchors"]
    for slide in args.slides:
        t0 = time.time()
        d = load(slide)
        h = load_hurdle(slide, d["ids"])
        S = (h["q"] * np.maximum(h.pop("mu"), 0)).astype(np.float32)
        del h["q"]
        ct = h["ctp"].argmax(1)
        Y, P, xy = d["y"], d["pred"], d["xy_um"]
        _, groups = gene_groups(d, slide)
        WS, HAS = {}, {}
        for r in RADII:
            W, HAS[r] = same_type_matrix(xy, ct, r)
            WS[r] = (W @ S).astype(np.float32)
            del W
        print(f"[{slide}] graphs built; cells with same-type neighbours: "
              + ", ".join(f"r={r:g}: {HAS[r].mean():.2f}" for r in RADII) + f" ({time.time()-t0:.0f}s)", flush=True)
        featS = slide_stats(S)
        tree = cKDTree(xy); pcs = prediction_pcs(P)
        out = HERE / ("outputs_v8_smoke" if args.smoke else "outputs_v8") / slide
        out.mkdir(parents=True, exist_ok=True)
        rows, genes_tab = [], []
        configs = [(1.0, 4, "diverse", 0)] if args.smoke else [(ar, n, pl, s) for ar, n, pl in CONFIGS for s in a["repeats"]]
        for ar, n, pl, s in configs:
            centres, rad, aidx, ag = place_cores(xy, ar, n, pl, s, P=P, pcs=pcs, tree=tree)
            eidx = np.flatnonzero(evaluation_mask(xy, centres, rad, a["evaluation_buffer_um"]))
            lam, r = choose_niche(S, WS, HAS, Y, aidx)
            Sn = niche(S, WS[r], HAS[r], lam)
            featN = slide_stats(Sn)
            rrs = adapt.fit_rrc_authors(S[aidx], Y[aidx], *featS)
            rrn = adapt.fit_rrc_authors(Sn[aidx], Y[aidx], *featN)
            k_s = choose_k(rrc_oof(S[aidx], Y[aidx], ag, featS, rrs), Y[aidx], ag)
            k_n = choose_k(rrc_oof(Sn[aidx], Y[aidx], ag, featN, rrn), Y[aidx], ag)
            Ye = Y[eidx]
            q_rs = adapt.predict_rrc_authors(rrs, S[eidx])
            q_rn = adapt.predict_rrc_authors(rrn, Sn[eidx])
            preds = {"GHIST+": P[eidx], "GHIST+-Soft": S[eidx], "Niche-Soft": Sn[eidx], "RRC-soft": q_rs,
                     "Niche-RRC": q_rn,
                     "Prog-RRC": q_rs if k_s is None else prog_apply(q_rs, prog_fit(Y[aidx], k_s)),
                     "Niche-Prog-RRC": q_rn if k_n is None else prog_apply(q_rn, prog_fit(Y[aidx], k_n))}
            res, pg = evaluate(Ye, preds, groups, xy[eidx], s, refs=("GHIST+-Soft", "RRC-soft"))
            key = f"a{ar}_n{n}_{pl}_s{s}"
            for meth, rw in res.items():
                rw.update(slide=slide, config=key, seed=s, method=meth, n_anchor_cells=int(len(aidx)),
                          n_eval_cells=int(len(eidx)), niche_lam=lam, niche_r=r, prog_k_soft=k_s or 0, prog_k_niche=k_n or 0)
                rows.append(rw)
                genes_tab.append(pd.DataFrame({"config": key, "method": meth, "gene_idx": np.arange(len(pg[meth])), "pcc": pg[meth]}))
            f = lambda m: res[m]["median_pcc_all"]
            print(f"[{slide}] {key}: lam={lam} r={r:g} k={k_s}/{k_n} | medPCC Soft {f('GHIST+-Soft'):.4f} Niche-Soft {f('Niche-Soft'):.4f} "
                  f"RRC-soft {f('RRC-soft'):.4f} Niche-RRC {f('Niche-RRC'):.4f} Prog-RRC {f('Prog-RRC'):.4f} Both {f('Niche-Prog-RRC'):.4f} | "
                  f"CMD RRC-soft {res['RRC-soft']['cmd_all']:.3f} Prog {res['Prog-RRC']['cmd_all']:.3f} | "
                  f"SSIM RRC-soft {res['RRC-soft']['median_ssim_svg50']:.3f} Niche-RRC {res['Niche-RRC']['median_ssim_svg50']:.3f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
            pd.DataFrame(rows).to_csv(out / "runs.csv", index=False, float_format="%.6g")
            pd.concat(genes_tab).to_csv(out / "per_gene.csv.gz", index=False, float_format="%.5g")
            del preds, res, pg, Ye, Sn, q_rs, q_rn; gc.collect()
        print(f"[{slide}] done ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
