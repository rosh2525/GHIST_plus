"""v9 ablations of BioHAN (protocol config_v9.json, ablation section): which biological component carries the gain?

Variants, all on identical cores and evaluation cells, each refit with the same training procedure:
  BioHAN              full model (reference)
  BioHAN-nobottleneck no 64-factor gene-programme bottleneck (trunk output width 512)
  BioHAN-noweight     no rare-cell-type up-weighting (uniform sampling)
  BioHAN-noniche      no niche inputs (re-reported for completeness)
Also re-reports RRC-soft. Seeds 0-2 of the primary core layout plus the paper-matched layout.

Usage:  python biohan_ablate.py breast5k [atera]
"""
from __future__ import annotations

import argparse
import gc
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.spatial import cKDTree

import adapt
import biohan_run as B
from anchors import evaluation_mask, place_cores, prediction_pcs
from hurdle import load_hurdle
from p1_run import evaluate, slide_stats
from run_v3 import CFG, gene_groups
from slides import load

HERE = Path(__file__).resolve().parent
LAYOUTS = [(1.0, 4, "diverse"), (3.5, 1, "diverse")]


def fit_variant(X, LQ, R, Y, ct, groups, seed, weighted=True, bottleneck=True):
    k_saved = B.K
    B.K = B.K if bottleneck else B.HID
    try:
        if weighted:
            return B.fit_biohan(X, LQ, R, Y, ct, groups, seed)
        wts = np.ones(len(X), np.float32)
        gs = np.unique(groups); val = groups == gs[seed % len(gs)]
        _, ep, _ = B.train(X[~val], LQ[~val], R[~val], Y[~val], wts[~val], None, X[val], LQ[val], R[val], Y[val], seed)
        net, _, _ = B.train(X, LQ, R, Y, wts, epochs=max(ep, 1), seed=seed)
        return net, max(ep, 1)
    finally:
        B.K = k_saved


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slides", nargs="+")
    args = ap.parse_args()
    a = CFG["anchors"]
    print("device:", B.DEV, flush=True)
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
        NS = (B.neighbour_mean(xy, 30.0, ct) @ S).astype(np.float32)
        NC = (B.neighbour_mean(xy, 50.0) @ ctp).astype(np.float32)
        featN, featC = slide_stats(NS), slide_stats(NC)
        z = lambda A, f: ((A - f[0]) / f[1]).astype(np.float32)
        base = lambda idx: [z(S[idx], featS), z(lq[idx], featL), ctp[idx]]
        X_full = lambda idx: np.hstack(base(idx) + [z(NS[idx], featN), z(NC[idx], featC)])
        X_noni = lambda idx: np.hstack(base(idx))
        tree = cKDTree(xy); pcs = prediction_pcs(P)
        out = HERE / "outputs_v9_ablate" / slide
        out.mkdir(parents=True, exist_ok=True)
        rows, genes_tab = [], []
        for ar, n, pl in LAYOUTS:
            for s in a["repeats"]:
                t1 = time.time()
                centres, rad, aidx, ag = place_cores(xy, ar, n, pl, s, P=P, pcs=pcs, tree=tree)
                eidx = np.flatnonzero(evaluation_mask(xy, centres, rad, a["evaluation_buffer_um"]))
                Ya = Y[aidx]
                rrs = adapt.fit_rrc_authors(S[aidx], Ya, *featS)
                Xa, Xe = X_full(aidx), X_full(eidx)
                nets = {"BioHAN": fit_variant(Xa, lq[aidx], R[aidx], Ya, ct[aidx], ag, s),
                        "BioHAN-nobottleneck": fit_variant(Xa, lq[aidx], R[aidx], Ya, ct[aidx], ag, s, bottleneck=False),
                        "BioHAN-noweight": fit_variant(Xa, lq[aidx], R[aidx], Ya, ct[aidx], ag, s, weighted=False)}
                net_n, ep_n = fit_variant(X_noni(aidx), lq[aidx], R[aidx], Ya, ct[aidx], ag, s)
                preds = {"RRC-soft": adapt.predict_rrc_authors(rrs, S[eidx])}
                for k, (net, _) in nets.items():
                    preds[k] = B.predict(net, Xe, lq[eidx], R[eidx])
                preds["BioHAN-noniche"] = B.predict(net_n, X_noni(eidx), lq[eidx], R[eidx])
                res, pg = evaluate(Y[eidx], preds, groups, xy[eidx], s, refs=("BioHAN", "RRC-soft"))
                key = f"a{ar}_n{n}_{pl}_s{s}"
                for meth, rw in res.items():
                    rw.update(slide=slide, config=key, seed=s, method=meth, n_anchor_cells=int(len(aidx)),
                              epochs={**{k: e for k, (_, e) in nets.items()}, "BioHAN-noniche": ep_n}.get(meth, 0))
                    rows.append(rw)
                    genes_tab.append(pd.DataFrame({"config": key, "method": meth, "gene_idx": np.arange(len(pg[meth])), "pcc": pg[meth]}))
                f = lambda m: res[m]["median_pcc_all"]
                print(f"[{slide}] {key}: RRC-soft {f('RRC-soft'):.4f} BioHAN {f('BioHAN'):.4f} nobottleneck {f('BioHAN-nobottleneck'):.4f} "
                      f"noweight {f('BioHAN-noweight'):.4f} noniche {f('BioHAN-noniche'):.4f} | CMD full {res['BioHAN']['cmd_all']:.3f} "
                      f"nobottleneck {res['BioHAN-nobottleneck']['cmd_all']:.3f} ({time.time()-t1:.0f}s)", flush=True)
                pd.DataFrame(rows).to_csv(out / "runs.csv", index=False, float_format="%.6g")
                pd.concat(genes_tab).to_csv(out / "per_gene.csv.gz", index=False, float_format="%.5g")
                del preds, res, pg, nets, net_n; gc.collect(); torch.cuda.empty_cache()
        print(f"[{slide}] done ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
