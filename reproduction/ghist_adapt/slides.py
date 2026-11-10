"""Slide registry and matched-cohort caches for GHIST+-Adapt v3.

Cohort rule (authors' Figure3.ipynb): genes = GHIST+ prediction genes that are also measured,
cells = intersection of prediction / measured / coordinate IDs. Target = log1p(raw counts).
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BUNDLE = Path(os.getenv("GHIST_BUNDLE_ROOT") or REPO / "local_artifacts/bundle")
FETCH = Path(os.getenv("GHIST_V3_FETCH") or REPO / "local_artifacts/bundle_v3_fetch")  # re-downloaded copies
CACHE = Path(os.getenv("GHIST_ADAPT_CACHE") or REPO / "local_artifacts/ghist_adapt_cache")

SLIDES = {
    "breast5k": dict(
        pred=BUNDLE / "GHIST_plus/predictions/breast_multi/breast5k/test_slide6_pred_expr_scaled.npz",
        gt=BUNDLE / "evaluation_data/breast5k/cell_gene_matrix_filtered.csv",
        coords=BUNDLE / "evaluation_data/breast5k/cell_coords_histology_filtered.csv",
        pred_scale=0.5, um_per_unit=1.0,
        unit_note="micrometres: median nearest-neighbour distance 7.0 units"),
    "atera": dict(
        pred=BUNDLE / "GHIST_plus/predictions/breast_multi/atera/test_slide9101_pred_expr_scaled.npz",
        gt=BUNDLE / "evaluation_data/atera/cell_gene_matrix_filtered.csv",
        coords=BUNDLE / "evaluation_data/atera/cell_coords.csv",
        pred_scale=0.5, um_per_unit=0.2125,
        unit_note="ASSUMED H&E pixels at 0.2125 um: nearest-neighbour distance 28.3 units matches BreastCancer2"),
    "breast5k_pan": dict(  # same slide, predicted by the released PAN-CANCER checkpoint (Breast5K = its held-out slide 11)
        pred=BUNDLE / "GHIST_plus/predictions/pancancer/breast_5k/test_slide11_pred_expr_scaled.npz",
        gt=BUNDLE / "evaluation_data/breast5k/cell_gene_matrix_filtered.csv",
        coords=BUNDLE / "evaluation_data/breast5k/cell_coords_histology_filtered.csv",
        pred_scale=0.5, um_per_unit=1.0,   # pancancer config: data.expr_scale = 2.0
        unit_note="micrometres (same coordinates as breast5k)", expected=None),
    "breast2": dict(
        pred=BUNDLE / "GHIST_plus/predictions/breast_multi/test_slide3_pred_expr_scaled.csv",
        gt=FETCH / "evaluation_data/breast2/cell_gene_matrix_filtered.csv",
        coords=FETCH / "evaluation_data/breast2/cell_coords.csv",
        pred_scale=1.0, um_per_unit=0.2125,
        unit_note="H&E pixels at 0.2125 um (paper: 500 pixels ~ 106 um)"),
}


def _csv_ids_genes(path):
    cols = pd.read_csv(path, nrows=0).columns.astype(str).str.strip().tolist()
    return cols[0], cols[1:]


def _read_matrix_csv(path, id_col, genes, ids):
    row = {v: i for i, v in enumerate(ids)}
    out = np.zeros((len(ids), len(genes)), dtype=np.float32)
    seen = 0
    for chunk in pd.read_csv(path, usecols=[id_col, *genes], chunksize=50000, dtype={g: np.float32 for g in genes}):
        cid = pd.to_numeric(chunk[id_col], errors="raise").astype(np.int64).to_numpy()
        keep = np.isin(cid, ids)
        rows = np.fromiter((row[v] for v in cid[keep]), dtype=np.int64)
        out[rows] = chunk.loc[keep, genes].to_numpy(np.float32)
        seen += int(keep.sum())
    if seen != len(ids):
        raise AssertionError(f"{path}: matched {seen} rows, expected {len(ids)}")
    return out


def build(name: str) -> Path:
    s = SLIDES[name]
    t0 = time.time()
    if s["pred"].suffix == ".npz":
        with np.load(s["pred"], allow_pickle=True) as z:  # gene_names is an object array
            pid = z["cell_ids"].astype(np.int64)
            pgenes = z["gene_names"].astype(str).tolist()
            pvals = z["pred_expr_scaled"]
    else:
        pid_col, pgenes = _csv_ids_genes(s["pred"])
        pid = pd.read_csv(s["pred"], usecols=[pid_col])[pid_col].astype(np.int64).to_numpy()
        pvals = None
    gid_col, ggenes = _csv_ids_genes(s["gt"])
    gset = set(ggenes)
    genes = [g for g in pgenes if g in gset]
    gids = pd.read_csv(s["gt"], usecols=[gid_col])[gid_col].astype(np.int64).to_numpy()
    c = pd.read_csv(s["coords"])
    c = c.rename(columns={c.columns[0]: "id", "x_coord": "x", "y_coord": "y"})
    c["id"] = c["id"].astype(np.int64)
    for arr, what in ((pid, "prediction"), (gids, "measured"), (c["id"].to_numpy(), "coordinate")):
        if len(arr) != len(np.unique(arr)):
            raise ValueError(f"{name}: duplicate {what} cell IDs")
    ids = np.intersect1d(np.intersect1d(pid, gids), c["id"].to_numpy())
    raw = _read_matrix_csv(s["gt"], gid_col, genes, ids)
    if pvals is not None:
        order = np.argsort(pid)
        rows = order[np.searchsorted(pid[order], ids)]
        cols = [pgenes.index(g) for g in genes]
        pred = pvals[np.ix_(rows, cols)].astype(np.float32)
    else:
        pred = _read_matrix_csv(s["pred"], pid_col, genes, ids)
    pred *= np.float32(s["pred_scale"])
    xy = c.set_index("id").loc[ids, ["x", "y"]].to_numpy(np.float64) * s["um_per_unit"]
    if not (np.isfinite(raw).all() and (raw >= 0).all() and np.isfinite(pred).all() and np.isfinite(xy).all()):
        raise ValueError(f"{name}: invalid values")
    CACHE.mkdir(parents=True, exist_ok=True)
    out = CACHE / f"{name}.npz"
    np.savez(out, ids=ids, genes=np.asarray(genes, dtype="U32"), raw=raw, pred=pred, xy_um=xy.astype(np.float32))
    print(f"{name}: {len(ids)} cells x {len(genes)} genes cached ({time.time() - t0:.0f}s) -> {out}")
    return out


def load(name: str) -> dict:
    path = CACHE / f"{name}.npz"
    if not path.exists():
        build(name)
    with np.load(path) as z:
        d = {k: z[k] for k in z.files}
    d["genes"] = d["genes"].astype(str)
    d["xy_um"] = d["xy_um"].astype(np.float64)
    d["y"] = np.log1p(d.pop("raw"))
    return d


if __name__ == "__main__":
    import sys
    for n in sys.argv[1:] or list(SLIDES):
        build(n)
