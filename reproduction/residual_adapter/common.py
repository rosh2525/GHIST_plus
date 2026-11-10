"""Shared paths, loaders and metrics for the GHIST+ residual-adapter experiment.

Loading mirrors the authors' Figure3.ipynb (`load_external_scope`) for Breast5K:
genes = prediction genes that also appear in the measured matrix, cells = the
intersection of prediction / measured / coordinate IDs, target = log1p(raw counts),
prediction = released `pred_expr_scaled` * 0.5.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BUNDLE = Path(os.getenv("GHIST_BUNDLE_ROOT") or REPO / "local_artifacts" / "bundle").expanduser()
CACHE_DIR = Path(os.getenv("ADAPTER_CACHE") or REPO / "local_artifacts" / "adapter_cache").expanduser()
CACHE = CACHE_DIR / "breast5k_matched.npz"
SPLIT_DIR = HERE / "splits"
OUT_DIR = HERE / "outputs"
CONFIG_PATH = HERE / "config.json"

PATHS = {
    "prediction": BUNDLE / "GHIST_plus/predictions/breast_multi/breast5k/test_slide6_pred_expr_scaled.npz",
    "gt": BUNDLE / "evaluation_data/breast5k/cell_gene_matrix_filtered.csv",
    "coordinates": BUNDLE / "evaluation_data/breast5k/cell_coords_histology_filtered.csv",
}
EXPECTED = (564723, 258)
AUTHOR_PRED_SCALE = 0.5  # Figure3.ipynb load_external_scope


def load_config() -> dict:
    return json.loads(CONFIG_PATH.read_text())


def load_cache():
    with np.load(CACHE, allow_pickle=False) as z:
        data = {k: z[k] for k in z.files}
    data["genes"] = data["genes"].astype(str)
    data["gt"] = np.log1p(data["raw"].astype(np.float32))
    data["pred"] = data["pred_unscaled"] * np.float32(AUTHOR_PRED_SCALE)
    return data


def load_split(seed: int) -> pd.DataFrame:
    return pd.read_csv(SPLIT_DIR / f"split_seed{seed}.csv.gz")


def rows_for(ids_all: np.ndarray, ids_subset: np.ndarray) -> np.ndarray:
    """Row indices of ids_subset inside the sorted id vector ids_all."""
    idx = np.searchsorted(ids_all, ids_subset)
    if not np.array_equal(ids_all[idx], ids_subset):
        raise KeyError("split contains cell IDs absent from the cache")
    return idx


# ---------------------------------------------------------------- metrics
def pcc_per_gene(y: np.ndarray, p: np.ndarray) -> np.ndarray:
    y = y.astype(np.float64) - y.mean(0)
    p = p.astype(np.float64) - p.mean(0)
    den = np.sqrt((y * y).sum(0) * (p * p).sum(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (y * p).sum(0) / den
    out[den == 0] = np.nan
    return out


def rmse_per_gene(y: np.ndarray, p: np.ndarray) -> np.ndarray:
    return np.sqrt(((y.astype(np.float64) - p) ** 2).mean(0))


def positive_rmse_per_gene(y: np.ndarray, p: np.ndarray) -> np.ndarray:
    pos = y > 0
    sq = np.where(pos, (y.astype(np.float64) - p) ** 2, 0.0)
    n = pos.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, np.sqrt(sq.sum(0) / np.maximum(n, 1)), np.nan)


def std_ratio_per_gene(y: np.ndarray, p: np.ndarray) -> np.ndarray:
    sy = y.std(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(sy > 0, p.std(0) / sy, np.nan)
