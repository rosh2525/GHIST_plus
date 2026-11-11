"""Spillover test: does Xenium's measured RNA contain transcripts leaking from touching neighbours,
and does modelling that leak improve GHIST+ predictions?

Decision rule (written before the first run, 2026-11-11):
  D1 (measured data only). Immune cells are defined by measured RNA (CD3E, MS4A1 or MZB1 > 0) and tumour cells by
     GHIST+'s H&E cell-type head (Malignant + Epithelial prob >= 0.6). Spillover is SUPPORTED if detection of the
     epithelial markers in immune cells 0-10 um from tumour is >= 2x the detection at 30-50 um, while immune markers
     in the same cells change by < 1.5x.
  D2 (prediction). Base = soft output S = q * relu(mu). Per gene, fit y ~ S + contact-neighbour mean of S
     (+ far-ring mean as control) on half of 400 um spatial blocks, test on the other half.
     Spillover modelling is WORTH BUILDING if the contact term raises median test PCC (all genes) by >= 0.005
     and by more than the far-ring term does.
Development slide: breast5k. ATERA is run afterwards with this file unchanged.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from slides import load  # noqa: E402
from hurdle import load_hurdle  # noqa: E402

EPI = ["EPCAM", "CDH1", "ERBB2", "GATA3", "FOXA1", "ESR1", "AGR3", "KRT8", "KRT18", "KRT19", "KRT7", "TFF1", "ELF3"]
IMM_DEF = ["CD3E", "MS4A1", "MZB1"]
IMM_CTRL = ["CD3E", "CD8A", "MS4A1", "MZB1", "CD68", "CD163", "PTPRC", "CD2", "CD3D"]
TUMOUR_CT = (4, 5)  # Malignant, Epithelial in the ct_probs order
BINS = [(0, 10), (10, 15), (15, 20), (20, 30), (30, 50), (50, 100)]
CONTACT = (0.0, 15.0)   # um, centre-to-centre: touching cells (median nearest-neighbour ~7 um on Breast5K)
FAR = (25.0, 40.0)      # um, same tissue neighbourhood but not touching


def ring_matrix(xy, lo, hi):
    tree = cKDTree(xy)
    pairs = tree.query_pairs(hi, output_type="ndarray")
    d = np.linalg.norm(xy[pairs[:, 0]] - xy[pairs[:, 1]], axis=1)
    pairs = pairs[(d > lo) & (d <= hi)]
    n = len(xy)
    r = np.concatenate([pairs[:, 0], pairs[:, 1]])
    c = np.concatenate([pairs[:, 1], pairs[:, 0]])
    W = sparse.csr_matrix((np.ones(len(r), np.float32), (r, c)), shape=(n, n))
    deg = np.asarray(W.sum(1)).ravel()
    inv = np.where(deg > 0, 1.0 / np.maximum(deg, 1), 0).astype(np.float32)
    return sparse.diags(inv) @ W, deg


def d1(d, h, out):
    genes = list(d["genes"])
    y, xy = d["y"], d["xy_um"]
    tumour = h["ctp"][:, list(TUMOUR_CT)].sum(1) >= 0.6
    imm_cols = [genes.index(g) for g in IMM_DEF if g in genes]
    immune = (y[:, imm_cols] > 0).any(1) & ~tumour
    dist, _ = cKDTree(xy[tumour]).query(xy[immune], k=1)
    epi = [g for g in EPI if g in genes]
    ctrl = [g for g in IMM_CTRL if g in genes]
    rows = []
    for lo, hi in BINS:
        m = (dist > lo) & (dist <= hi)
        Yi = y[immune][m]
        rows.append({"bin_um": f"{lo}-{hi}", "n_cells": int(m.sum()),
                     **{f"det_{g}": float((Yi[:, genes.index(g)] > 0).mean()) if m.any() else None for g in epi + ctrl}})
    near, far = rows[0], rows[4]
    ratio = lambda g: (near[f"det_{g}"] / far[f"det_{g}"]) if far[f"det_{g}"] else float("inf")
    epi_ratio = float(np.median([ratio(g) for g in epi]))
    ctrl_ratio = float(np.median([ratio(g) for g in ctrl]))
    res = {"n_tumour": int(tumour.sum()), "n_immune": int(immune.sum()), "epi_genes": epi, "ctrl_genes": ctrl,
           "bins": rows, "median_epi_ratio_0_10_vs_30_50": epi_ratio,
           "median_ctrl_ratio_0_10_vs_30_50": ctrl_ratio,
           "D1_supported": bool(epi_ratio >= 2.0 and ctrl_ratio < 1.5)}
    print(f"D1: immune {res['n_immune']}, tumour {res['n_tumour']}; epithelial-marker detection ratio "
          f"(0-10 vs 30-50 um) = {epi_ratio:.2f}, immune-marker ratio = {ctrl_ratio:.2f} -> "
          f"{'SUPPORTED' if res['D1_supported'] else 'not supported'}")
    for r in rows:
        print("  ", r["bin_um"], r["n_cells"], " ".join(f"{g}={r[f'det_{g}']:.3f}" for g in epi[:6] if r[f"det_{g}"] is not None))
    return res


def pcc_cols(Y, P):
    Y = Y - Y.mean(0)
    P = P - P.mean(0)
    den = np.sqrt((Y * Y).sum(0) * (P * P).sum(0))
    return np.where(den > 0, (Y * P).sum(0) / np.maximum(den, 1e-12), np.nan)


def d2(d, h, out, seed=0):
    y, xy = d["y"], d["xy_um"]
    S = (h["q"] * np.maximum(h["mu"], 0)).astype(np.float32)
    Wc, degc = ring_matrix(xy, *CONTACT)
    Wf, degf = ring_matrix(xy, *FAR)
    C = (Wc @ S).astype(np.float32)
    F = (Wf @ S).astype(np.float32)
    blk = np.floor(xy / 400.0).astype(np.int64)
    key = blk[:, 0] * 100003 + blk[:, 1]
    ub = np.unique(key)
    rng = np.random.default_rng(seed)
    train_blocks = rng.permutation(ub)[: len(ub) // 2]
    tr = np.isin(key, train_blocks)
    te = ~tr
    G = y.shape[1]
    preds = {k: np.zeros((te.sum(), G), np.float32) for k in ("S", "S+contact", "S+far", "S+contact+far")}
    coef_c = np.zeros(G)
    for g in range(G):
        feats = {"S": [S[:, g]], "S+contact": [S[:, g], C[:, g]], "S+far": [S[:, g], F[:, g]],
                 "S+contact+far": [S[:, g], C[:, g], F[:, g]]}
        for k, cols in feats.items():
            X = np.column_stack([np.ones(len(y))] + cols)
            beta, *_ = np.linalg.lstsq(X[tr], y[tr, g], rcond=None)
            preds[k][:, g] = X[te] @ beta
            if k == "S+contact+far":
                coef_c[g] = beta[2]
    Yte = y[te]
    pcc = {k: pcc_cols(Yte, v) for k, v in preds.items()}
    med = {k: float(np.nanmedian(v)) for k, v in pcc.items()}
    gain_c = med["S+contact"] - med["S"]
    gain_f = med["S+far"] - med["S"]
    res = {"median_pcc": med, "gain_contact": gain_c, "gain_far": gain_f,
           "gain_contact_given_far": med["S+contact+far"] - med["S+far"],
           "genes_with_positive_contact_coef": int((coef_c > 0).sum()), "n_genes": G,
           "median_contact_degree": float(np.median(degc)), "median_far_degree": float(np.median(degf)),
           "n_train": int(tr.sum()), "n_test": int(te.sum()),
           "D2_worth_building": bool(gain_c >= 0.005 and gain_c > gain_f)}
    print(f"D2: median test PCC {json.dumps({k: round(v, 4) for k, v in med.items()})}")
    print(f"    contact gain {gain_c:+.4f}, far gain {gain_f:+.4f}, contact beyond far {res['gain_contact_given_far']:+.4f};"
          f" contact coef > 0 in {res['genes_with_positive_contact_coef']}/{G} genes -> "
          f"{'WORTH BUILDING' if res['D2_worth_building'] else 'not worth building'}")
    np.save(out / "d2_pcc_per_gene.npy", np.vstack([pcc[k] for k in preds]))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slide", choices=["breast5k", "atera", "breast5k_pan"])
    ap.add_argument("--out", default=str(HERE / "outputs"))
    a = ap.parse_args()
    out = Path(a.out) / a.slide
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    d = load(a.slide)
    h = load_hurdle(a.slide, d["ids"])
    print(f"{a.slide}: {len(d['ids'])} cells x {len(d['genes'])} genes loaded ({time.time() - t0:.0f}s)")
    res = {"slide": a.slide, "D1": d1(d, h, out), "D2": d2(d, h, out)}
    (out / "spill_test.json").write_text(json.dumps(res, indent=1))
    print(f"done in {time.time() - t0:.0f}s -> {out / 'spill_test.json'}")


if __name__ == "__main__":
    main()
