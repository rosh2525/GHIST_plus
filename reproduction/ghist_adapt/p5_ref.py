"""Wu et al. 2021 breast cancer scRNA-seq atlas (GSE176078), restricted to the Xenium panel genes.

Builds ghist_adapt_cache/wu2021_panel.npz: log1p(CP10k) expression (cells x panel genes, float32), cell types,
patient IDs. Usage: python p5_ref.py   (reads work/reference/count_matrix_sparse.mtx, ~2.4 GB, streamed once)
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

from slides import CACHE, load

REF = Path(os.getenv("WU2021_DIR") or Path(__file__).resolve().parents[3] / "work" / "reference")
OUT = CACHE / "wu2021_panel.npz"


def build():
    t0 = time.time()
    panel = sorted(set(load("breast5k")["genes"]) | set(load("atera")["genes"]))
    genes = pd.read_csv(REF / "count_matrix_genes.tsv", header=None)[0].astype(str).tolist()
    gidx = {g: i for i, g in enumerate(genes)}
    keep = [g for g in panel if g in gidx]
    row_to_col = np.full(len(genes) + 1, -1, dtype=np.int64)          # mtx rows are 1-based gene indices
    for j, g in enumerate(keep):
        row_to_col[gidx[g] + 1] = j
    meta = pd.read_csv(REF / "metadata.csv", index_col=0)
    n_cells = len(meta)
    X = np.zeros((n_cells, len(keep)), dtype=np.float32)
    for chunk in pd.read_csv(REF / "count_matrix_sparse.mtx", sep=" ", skiprows=2, header=None, names=["g", "c", "v"],
                             dtype={"g": np.int32, "c": np.int32, "v": np.float32}, chunksize=20_000_000):
        col = row_to_col[chunk["g"].to_numpy()]
        m = col >= 0
        X[chunk["c"].to_numpy()[m] - 1, col[m]] += chunk["v"].to_numpy()[m]
        print(f"  {time.time()-t0:.0f}s", flush=True)
    lib = meta["nCount_RNA"].to_numpy(np.float32)
    X = np.log1p(X / lib[:, None] * 1e4).astype(np.float32)
    np.savez(OUT, X=X, genes=np.asarray(keep, dtype="U32"), celltype_major=meta["celltype_major"].to_numpy().astype("U64"),
             celltype_minor=meta["celltype_minor"].to_numpy().astype("U64"), patient=meta["orig.ident"].to_numpy().astype("U32"))
    print(f"wu2021: {n_cells} cells x {len(keep)}/{len(panel)} panel genes -> {OUT} ({time.time()-t0:.0f}s)")


def load_ref():
    if not OUT.exists():
        build()
    with np.load(OUT) as z:
        return {k: z[k] for k in z.files}


if __name__ == "__main__":
    build()
