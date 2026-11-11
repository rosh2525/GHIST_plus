"""Score a retrained GHIST (or any) within-slide prediction the way the paper's Figure3 does.

Cohort (Figure3.ipynb load_scope, "within"): cells and genes common to measured RNA, coordinates and every
released within-slide prediction -> expected 12,719 cells x 280 genes. Target = log1p(raw counts).
Gene groups: Giotto spatial score (authors' giotto_scores, k = 8) on the measured cohort; headline metric =
median per-gene PCC over the top-20 spatially variable genes. Paper (manuscript): GHIST+ 0.744, GHIST 0.714.

Usage: python eval_within.py --bundle BUNDLE --ours NAME=path.csv [NAME=path.csv ...] --out results.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from metrics import giotto_scores, pcc_columns  # noqa: E402

ID_COLUMNS = ["c_id", "cell_id", "cellID", "Cell ID", "id", "V1"]
NON_GENES = {"x", "y", "x_coord", "y_coord", "x coordinate", "y coordinate"}
PAPER = {"GHIST+": 0.744, "GHIST": 0.714, "SpatialEx": 0.521}


def read(path):
    df = pd.read_csv(path)
    key = next((c for c in ID_COLUMNS if c in df.columns), df.columns[0])
    df[key] = pd.to_numeric(df[key], errors="raise").astype(np.int64)
    df = df.set_index(key)
    return df[[c for c in df.columns if str(c).strip() not in NON_GENES]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--ours", nargs="*", default=[])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    B = Path(a.bundle)
    released = {"GHIST+": (B / "GHIST_plus/predictions/breast_single/trainval_slide3_pred_expr_scaled.csv", 0.5),
                "GHIST": (B / "other_models/within_slide/GHIST_pred_expr.csv", 1.0),
                "SpatialEx": (B / "other_models/within_slide/SpatialEx_UNI2H_pred_expr.csv", 1.0)}
    gt = read(B / "evaluation_data/breast2/cell_gene_matrix_filtered.csv")
    co = pd.read_csv(B / "evaluation_data/breast2/cell_coords.csv")
    co = co.rename(columns={co.columns[0]: "id", "x_coord": "x", "y_coord": "y"}).set_index("id")
    preds = {k: (read(p), s) for k, (p, s) in released.items()}
    ids = set(gt.index) & set(co.index)
    genes = set(gt.columns)
    for df, _ in preds.values():
        ids &= set(df.index); genes &= set(df.columns)
    ids = np.array(sorted(ids)); genes = sorted(genes)
    print(f"released cohort: {len(ids)} cells x {len(genes)} genes (paper notebook expects 12719 x 280)")
    Y = np.log1p(gt.loc[ids, genes].to_numpy(np.float64))
    xy = co.loc[ids, ["x", "y"]].to_numpy(np.float32)
    score = giotto_scores(Y.astype(np.float32), xy)
    top20 = np.argsort(-score, kind="stable")[:20]
    res = {"cohort_cells": int(len(ids)), "cohort_genes": len(genes), "top20_svg": [genes[i] for i in top20], "methods": {}}
    for name, (df, s) in preds.items():
        pcc = pcc_columns(Y, df.loc[ids, genes].to_numpy(np.float64) * s)
        res["methods"][name] = {"median_pcc_top20_svg": float(np.nanmedian(pcc[top20])),
                                "median_pcc_all": float(np.nanmedian(pcc)), "paper_top20": PAPER.get(name), "cells": int(len(ids))}
    for spec in a.ours:
        name, path = spec.split("=", 1)
        df = read(path)
        keep = np.isin(ids, df.index.to_numpy())
        miss = [g for g in genes if g not in df.columns]
        if miss:
            raise ValueError(f"{name}: missing genes {miss[:5]}")
        pcc = pcc_columns(Y[keep], df.loc[ids[keep], genes].to_numpy(np.float64))
        ref = pcc_columns(Y[keep], preds["GHIST"][0].loc[ids[keep], genes].to_numpy(np.float64))
        res["methods"][name] = {"median_pcc_top20_svg": float(np.nanmedian(pcc[top20])), "median_pcc_all": float(np.nanmedian(pcc)),
                                "cells": int(keep.sum()), "released_GHIST_on_same_cells_top20": float(np.nanmedian(ref[top20])),
                                "paper_top20": PAPER["GHIST"]}
    for k, v in res["methods"].items():
        print(f"{k:28s} top-20 SVG median PCC {v['median_pcc_top20_svg']:.4f} | all genes {v['median_pcc_all']:.4f} | "
              f"cells {v['cells']} | paper {v.get('paper_top20')}")
    Path(a.out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
