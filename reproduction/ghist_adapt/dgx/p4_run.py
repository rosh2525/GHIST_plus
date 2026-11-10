"""P4b: image-aware hurdle correction (HAC-img). HAC with UNI2-h image features added to the occurrence model.

Run only after (1) extract_uni2h.py produced embeddings for the slide and its QA overlay was checked, and
(2) config_v6.json was frozen. Same cores, evaluation and metrics as P1 (config_v4.json).

Usage:  python dgx/p4_run.py breast5k --emb /scratch/.../breast5k_uni2h      (expects <emb>.npy and <emb>_ids.npy)
Writes  outputs_v6/<slide>/runs.csv, per_gene.csv.gz
"""
from __future__ import annotations

import argparse
import gc
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import adapt  # noqa: E402
from anchors import evaluation_mask, place_cores, prediction_pcs  # noqa: E402
from hurdle import load_hurdle  # noqa: E402
from p1_run import CONFIGS, _hac_predict, evaluate, hac_fit, slide_stats  # noqa: E402
from run_v3 import CFG, gene_groups  # noqa: E402
from slides import load  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slide")
    ap.add_argument("--emb", required=True)
    ap.add_argument("--n-pc", type=int, default=None)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    cfg6 = json.loads((HERE / "config_v6.json").read_text())
    if not cfg6.get("frozen"):
        raise SystemExit("config_v6.json is not frozen: set 'frozen' (date, before any P4 metric) first")
    n_pc = args.n_pc or cfg6["image_features"]["n_pc"]
    t0 = time.time()
    d = load(args.slide); h = load_hurdle(args.slide, d["ids"])
    Y, P, q = d["y"], d["pred"], h["q"]
    R = np.maximum(h["mu"], 0); S = (q * R).astype(np.float32)
    E = np.load(args.emb + ".npy", mmap_mode="r"); eid = np.load(args.emb + "_ids.npy")
    pos = pd.Series(np.arange(len(eid)), index=eid).reindex(d["ids"])
    if pos.isna().any():
        raise SystemExit(f"{int(pos.isna().sum())} cohort cells have no embedding")
    pos = pos.to_numpy(np.int64)
    rng = np.random.default_rng(0)
    sub = np.asarray(E[np.sort(pos[rng.choice(len(pos), min(50_000, len(pos)), replace=False)])], dtype=np.float64)
    mu_e, sd_e = sub.mean(0), sub.std(0); sd_e[sd_e == 0] = 1
    _, _, Vt = np.linalg.svd((sub - mu_e) / sd_e, full_matrices=False)
    V = Vt[:n_pc].T.astype(np.float32); del sub
    F = np.empty((len(pos), n_pc), dtype=np.float32)                   # RNA-free image PCs, computed in chunks
    for a in range(0, len(pos), 50_000):
        blk = np.asarray(E[pos[a:a + 50_000]], dtype=np.float32)
        F[a:a + 50_000] = ((blk - mu_e.astype(np.float32)) / sd_e.astype(np.float32)) @ V
    del E
    featP, featS, featF = slide_stats(P), slide_stats(S), slide_stats(F)
    lq = np.log(np.clip(q, 1e-4, 1 - 1e-4) / (1 - np.clip(q, 1e-4, 1 - 1e-4))).astype(np.float32); featL = slide_stats(lq)
    _, groups = gene_groups(d, args.slide)
    tree = cKDTree(d["xy_um"]); pcs = prediction_pcs(P)
    out = HERE / ("outputs_v6_smoke" if args.smoke else "outputs_v6") / args.slide
    out.mkdir(parents=True, exist_ok=True)
    rows, gt = [], []
    a = CFG["anchors"]
    configs = [(1.0, 4, "diverse", 0)] if args.smoke else [(ar, n, pl, s) for ar, n, pl in CONFIGS for s in a["repeats"]]

    def X(idx, with_img):
        parts = [(S[idx] - featS[0]) / featS[1], (lq[idx] - featL[0]) / featL[1]]
        if with_img:
            parts.append((F[idx] - featF[0]) / featF[1])
        return np.hstack(parts).astype(np.float32)

    for ar, n, pl, s in configs:
        centres, r, aidx, ag = place_cores(d["xy_um"], ar, n, pl, s, P=P, pcs=pcs, tree=tree)
        eidx = np.flatnonzero(evaluation_mask(d["xy_um"], centres, r, a["evaluation_buffer_um"]))
        rrc = adapt.fit_rrc_authors(P[aidx], Y[aidx], *featP)
        hac = hac_fit(X(aidx, False), q[aidx], R[aidx], Y[aidx], ag)
        hci = hac_fit(X(aidx, True), q[aidx], R[aidx], Y[aidx], ag)
        preds = {"GHIST+": P[eidx], "GHIST+-Soft": S[eidx], "RRC-authors": adapt.predict_rrc_authors(rrc, P[eidx]),
                 "HAC": _hac_predict(hac, X(eidx, False), q[eidx], R[eidx]),
                 "HAC-img": _hac_predict(hci, X(eidx, True), q[eidx], R[eidx])}
        res, pg = evaluate(Y[eidx], preds, groups, d["xy_um"][eidx], s, refs=("GHIST+", "RRC-authors", "HAC"))
        key = f"a{ar}_n{n}_{pl}_s{s}"
        for m, rw in res.items():
            rw.update(slide=args.slide, config=key, method=m, n_anchor_cells=int(len(aidx)), hac_alpha=hac["alpha"], hac_img_alpha=hci["alpha"], n_pc=n_pc)
            rows.append(rw)
            gt.append(pd.DataFrame({"config": key, "method": m, "gene_idx": np.arange(len(pg[m])), "pcc": pg[m]}))
        print(f"[{args.slide}] {key}: HAC {res['HAC']['median_pcc_all']:.4f} HAC-img {res['HAC-img']['median_pcc_all']:.4f} RRC {res['RRC-authors']['median_pcc_all']:.4f} "
              f"CI(img-HAC) [{res['HAC-img']['ci_vs_HAC_all_lo']:+.4f},{res['HAC-img']['ci_vs_HAC_all_hi']:+.4f}] ({time.time()-t0:.0f}s)", flush=True)
        pd.DataFrame(rows).to_csv(out / "runs.csv", index=False, float_format="%.6g")
        pd.concat(gt).to_csv(out / "per_gene.csv.gz", index=False, float_format="%.5g")
        del preds, res, pg; gc.collect()


if __name__ == "__main__":
    main()
