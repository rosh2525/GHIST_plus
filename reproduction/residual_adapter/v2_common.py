"""v2 helpers: extra frozen GHIST+ outputs, train pools, features, MLP definition."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from common import CACHE_DIR, HERE, PATHS, load_split, rows_for

CONFIG_V2 = HERE / "config_v2.json"
EXTRAS = CACHE_DIR / "breast5k_extras.npz"
SPLIT_V2 = HERE / "splits_v2"
OUT_V2 = HERE / "outputs_v2"
MODEL_DIR_V2 = CACHE_DIR.parent / "adapter_v2_models"  # *.pt is git-ignored; weights kept locally


def load_config_v2() -> dict:
    return json.loads(CONFIG_V2.read_text())


def load_extras(ids: np.ndarray, genes: np.ndarray) -> dict:
    if not EXTRAS.exists():
        with np.load(PATHS["prediction"], allow_pickle=True) as z:
            pid = z["cell_ids"].astype(np.int64)
            pg = z["gene_names"].astype(str)
            order = np.argsort(pid)
            order = order[np.searchsorted(pid[order], ids)]
            if not np.array_equal(pid[order], ids):
                raise KeyError("prediction IDs do not cover the cohort")
            gcol = [int(np.flatnonzero(pg == g)[0]) for g in genes]
            occ = z["occurrence_probability"][order][:, gcol].astype(np.float32)
            ctp = z["ct_probs"][order].astype(np.float32)
        np.savez(EXTRAS, ids=ids, occ=occ, ctp=ctp)
    with np.load(EXTRAS) as z:
        if not np.array_equal(z["ids"], ids):
            raise ValueError("extras cache is out of sync with the main cache")
        return {"occ": z["occ"], "ctp": z["ctp"]}


def block_index(d: dict, block_um: float) -> np.ndarray:
    x, y = d["x"].astype(np.float64), d["y"].astype(np.float64)
    return np.stack([np.floor((x - x.min()) / block_um), np.floor((y - y.min()) / block_um)], 1).astype(int)


def build_pool(d: dict, seed: int, block_um: float) -> pd.DataFrame:
    """Permuted train pool for one split seed; cached as a manifest so it is fixed and inspectable."""
    path = SPLIT_V2 / f"train_pool_seed{seed}.csv.gz"
    if path.exists():
        return pd.read_csv(path)
    sp = load_split(seed)
    ids = d["ids"]
    held_rows = rows_for(ids, np.sort(sp.loc[sp.split.isin(["val", "test"]), "cell_id"].to_numpy()))
    blk = block_index(d, block_um)
    held = np.unique(blk[held_rows], axis=0)
    far = np.ones(len(ids), bool)
    for h in held:
        far &= np.abs(blk - h).max(1) >= 2
    pool = np.flatnonzero(far)
    perm = np.random.default_rng(1000 + seed).permutation(pool)
    df = pd.DataFrame({"cell_id": ids[perm], "rank": np.arange(len(perm))})
    SPLIT_V2.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return df


def features(d: dict, ex: dict, rows: np.ndarray):
    P = d["pred"][rows].astype(np.float32)
    X = np.hstack([P, ex["occ"][rows], ex["ctp"][rows]]).astype(np.float32)
    return P, X


def make_mlp(n_in: int, n_out: int, hidden, dropout: float):
    import torch.nn as nn

    layers, prev = [], n_in
    for h in hidden:
        layers += [nn.Linear(prev, h), nn.GELU(), nn.Dropout(dropout)]
        prev = h
    last = nn.Linear(prev, n_out)
    nn.init.zeros_(last.weight)
    nn.init.zeros_(last.bias)
    return nn.Sequential(*layers, last)
