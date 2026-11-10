"""Match Breast5K cells/genes, verify prediction scale, build spatial splits.

Usage:  python prepare_data.py [--skip-cache]

Writes
  local_artifacts/adapter_cache/breast5k_matched.npz  (large, git-ignored)
  reproduction/residual_adapter/splits/split_seed{s}.csv.gz + split_summary.json
  reproduction/residual_adapter/outputs/scale_check.json
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pandas as pd

from common import (AUTHOR_PRED_SCALE, CACHE, CACHE_DIR, EXPECTED, OUT_DIR, PATHS, REPO,
                    SPLIT_DIR, load_cache, load_config, pcc_per_gene, rmse_per_gene)


def build_cache() -> None:
    t0 = time.time()
    with np.load(PATHS["prediction"], allow_pickle=True) as z:  # gene_names is an object array
        pred_ids = z["cell_ids"].astype(np.int64)
        pred_genes = z["gene_names"].astype(str)
        pred_all = z["pred_expr_scaled"]
    header = pd.read_csv(PATHS["gt"], nrows=0).columns.astype(str).str.strip().tolist()
    id_col = header[0]
    gt_genes = set(header[1:])
    genes = [g for g in pred_genes if g in gt_genes]

    coords = pd.read_csv(PATHS["coordinates"]).rename(columns={"cell_id": "id", "x_coord": "x", "y_coord": "y"})
    coords["id"] = coords["id"].astype(np.int64)
    gt_ids = pd.read_csv(PATHS["gt"], usecols=[id_col])[id_col].astype(np.int64).to_numpy()
    ids = np.intersect1d(np.intersect1d(pred_ids, gt_ids, assume_unique=True), coords["id"].to_numpy(), assume_unique=True)
    if (len(ids), len(genes)) != EXPECTED:
        raise AssertionError(f"cohort {(len(ids), len(genes))} != expected {EXPECTED}")
    print(f"cohort {len(ids)} cells x {len(genes)} genes ({time.time()-t0:.0f}s)")

    # measured counts, chunked to bound memory
    gt_row = {v: i for i, v in enumerate(ids)}
    raw = np.zeros((len(ids), len(genes)), dtype=np.float32)
    seen = 0
    for chunk in pd.read_csv(PATHS["gt"], usecols=[id_col, *genes], chunksize=50000, dtype={g: np.float32 for g in genes}):
        cid = chunk[id_col].astype(np.int64).to_numpy()
        keep = np.isin(cid, ids)
        rows = np.fromiter((gt_row[v] for v in cid[keep]), dtype=np.int64)
        raw[rows] = chunk.loc[keep, genes].to_numpy(np.float32)
        seen += keep.sum()
        print(f"  read {seen} matched rows ({time.time()-t0:.0f}s)", flush=True)
    if seen != len(ids):
        raise AssertionError("missing measured rows")

    prow = np.searchsorted(pred_ids[np.argsort(pred_ids)], ids)
    order = np.argsort(pred_ids)[prow]
    gcol = [int(np.flatnonzero(pred_genes == g)[0]) for g in genes]
    pred_unscaled = pred_all[np.ix_(order, gcol)].astype(np.float32)
    c = coords.set_index("id").loc[ids]
    if not (np.isfinite(raw).all() and (raw >= 0).all() and np.isfinite(pred_unscaled).all()):
        raise ValueError("non-finite or negative values")
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(CACHE, ids=ids, genes=np.asarray(genes, dtype="U32"), raw=raw,
             pred_unscaled=pred_unscaled,
             x=c["x"].to_numpy(np.float32), y=c["y"].to_numpy(np.float32))
    print(f"cache written: {CACHE} ({time.time()-t0:.0f}s)")


def scale_check(d: dict) -> dict:
    """Verify our loader reproduces the authors' released Breast5K per-gene PCC table (3i)."""
    pcc = pcc_per_gene(d["gt"], d["pred"])
    report = {"n_cells": int(len(d["ids"])), "n_genes": int(len(d["genes"])),
              "author_scale": AUTHOR_PRED_SCALE, "median_pcc_full_slide": float(np.nanmedian(pcc))}
    ref_path = REPO / "results/reproduction/Figure3/3i.csv"
    if ref_path.exists():
        ref = pd.read_csv(ref_path).set_index("gene")["PCC"]
        common = [g for g in d["genes"] if g in ref.index]
        ours = pd.Series(pcc, index=d["genes"]).loc[common]
        report["reference_3i_genes_matched"] = len(common)
        report["max_abs_pcc_diff_vs_3i"] = float(np.abs(ours - ref.loc[common]).max())
    # unit sanity: magnitudes in log1p space. Reported only; the 0.5 scale is the authors', not tuned here.
    rng = np.random.default_rng(0)
    sub = rng.choice(len(d["ids"]), 50000, replace=False)
    y = d["gt"][sub]
    report["mean_log1p_measured"] = float(y.mean())
    for s in (0.5, 1.0):
        p = d["pred_unscaled"][sub] * s
        report[f"mean_pred_scale_{s}"] = float(p.mean())
        report[f"median_rmse_scale_{s}_random50k"] = float(np.median(rmse_per_gene(y, p)))
    return report


def make_splits(d: dict, cfg: dict) -> dict:
    s = cfg["split"]
    x, y, ids = d["x"].astype(np.float64), d["y"].astype(np.float64), d["ids"]
    bx = np.floor((x - x.min()) / s["block_um"]).astype(int)
    by = np.floor((y - y.min()) / s["block_um"]).astype(int)
    key = bx * 100000 + by
    blocks, inverse, counts = np.unique(key, return_inverse=True, return_counts=True)
    bxy = np.stack([blocks // 100000, blocks % 100000], 1)
    eligible = np.flatnonzero(counts >= s["min_cells_per_block"])
    gap = s["buffer_blocks"] + 1  # Chebyshev distance required between different splits
    summary = {"n_blocks": int(len(blocks)), "n_eligible_blocks": int(len(eligible)),
               "median_cells_per_eligible_block": float(np.median(counts[eligible])), "seeds": {}}
    SPLIT_DIR.mkdir(parents=True, exist_ok=True)
    for seed in s["seeds"]:
        rng = np.random.default_rng(seed)
        order = rng.permutation(eligible)
        assign: dict[str, list[int]] = {"train": [], "val": [], "test": []}
        target = {"train": s["n_train"], "val": s["n_val"], "test": s["n_test"]}

        def far(b, others):
            if not others:
                return True
            return np.abs(bxy[others] - bxy[b]).max(1).min() >= gap

        for name in ("train", "val", "test"):
            others = [b for k, v in assign.items() if k != name for b in v]
            for b in order:
                if sum(counts[assign[name]]) >= target[name]:
                    break
                if b in others or b in assign[name] or not far(b, others):
                    continue
                assign[name].append(int(b))
            if sum(counts[assign[name]]) < target[name]:
                raise RuntimeError(f"seed {seed}: not enough blocks for {name}")
        rows = []
        for name, bl in assign.items():
            cand = np.flatnonzero(np.isin(inverse, bl))
            pick = np.sort(rng.choice(cand, target[name], replace=False))
            rows.append(pd.DataFrame({"cell_id": ids[pick], "split": name, "block": blocks[inverse[pick]]}))
        dev = assign["train"] + assign["val"]
        dist = np.abs(bxy[:, None, :] - bxy[None, dev, :]).max(2).min(1)
        rest_rows = np.flatnonzero(dist[inverse] >= gap)
        rows.append(pd.DataFrame({"cell_id": ids[rest_rows], "split": "rest", "block": blocks[inverse[rest_rows]]}))
        df = pd.concat(rows, ignore_index=True)
        df.to_csv(SPLIT_DIR / f"split_seed{seed}.csv.gz", index=False)
        summary["seeds"][str(seed)] = {k: {"n_blocks": len(v), "n_cells_in_blocks": int(counts[v].sum()),
                                            "n_cells_used": int(target[k]), "blocks": [int(b) for b in blocks[v]]}
                                       for k, v in assign.items()}
        summary["seeds"][str(seed)]["rest_cells"] = int(len(rest_rows))
        print(f"seed {seed}: " + ", ".join(f"{k}={len(v)} blocks" for k, v in assign.items()) + f", rest={len(rest_rows)} cells")
    (SPLIT_DIR / "split_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-cache", action="store_true", help="reuse an existing cache")
    args = ap.parse_args()
    if not (args.skip_cache and CACHE.exists()):
        build_cache()
    d = load_cache()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rep = scale_check(d)
    (OUT_DIR / "scale_check.json").write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))
    if "max_abs_pcc_diff_vs_3i" in rep and rep["max_abs_pcc_diff_vs_3i"] > 1e-4:
        raise AssertionError("loader does not reproduce the authors' released Breast5K PCC table")
    make_splits(d, load_config())


if __name__ == "__main__":
    main()
