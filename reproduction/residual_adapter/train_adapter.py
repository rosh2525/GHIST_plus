"""Fit per-gene calibration and the multigene residual ridge adapter for every split.

Only training cells are used for fitting; only validation cells for choosing alpha.
Test and 'rest' cells are not loaded here.

Usage:  python train_adapter.py
Writes  outputs/seed{s}/{calibration.csv, adapter_weights.csv.gz, adapter_scaler.csv, val_curve.csv, fit_summary.json}
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from common import OUT_DIR, load_cache, load_config, load_split, pcc_per_gene, rmse_per_gene, rows_for


# ------------------------------------------------------------- models
def fit_calibration(P: np.ndarray, Y: np.ndarray):
    """Independent per-gene OLS: y_g ~ a_g + b_g * p_g."""
    pm, ym = P.mean(0), Y.mean(0)
    var = ((P - pm) ** 2).mean(0)
    cov = ((P - pm) * (Y - ym)).mean(0)
    b = np.where(var > 0, cov / np.where(var > 0, var, 1), 0.0)
    a = ym - b * pm
    return a, b


def apply_calibration(P, a, b):
    return np.clip(a + b * P, 0, None)


def apply_adapter(P, mu, sd, W, c):
    """Corrected = clip(P + c + ((P - mu) / sd) @ W, 0)."""
    return np.clip(P + c + ((P - mu) / sd) @ W, 0, None)


def ridge_path(Z: np.ndarray, R: np.ndarray, alphas):
    """Closed-form multi-output ridge for many alphas via one SVD (Z already train-standardised)."""
    zm, rm = Z.mean(0), R.mean(0)
    U, S, Vt = np.linalg.svd(Z - zm, full_matrices=False)
    UtR = U.T @ (R - rm)
    for alpha in alphas:
        W = Vt.T @ ((S / (S ** 2 + alpha))[:, None] * UtR)
        yield alpha, W, rm - zm @ W


def summarise(Y, Q):
    pcc, rmse = pcc_per_gene(Y, Q), rmse_per_gene(Y, Q)
    return {"median_pcc": float(np.nanmedian(pcc)), "mean_pcc": float(np.nanmean(pcc)),
            "mean_rmse": float(rmse.mean()), "median_rmse": float(np.median(rmse))}


def main() -> None:
    cfg = load_config()
    d = load_cache()
    ids, genes = d["ids"], d["genes"]
    for seed in cfg["split"]["seeds"]:
        sp = load_split(seed)
        tr = rows_for(ids, np.sort(sp.loc[sp.split == "train", "cell_id"].to_numpy()))
        va = rows_for(ids, np.sort(sp.loc[sp.split == "val", "cell_id"].to_numpy()))
        Ptr, Ytr = d["pred"][tr].astype(np.float64), d["gt"][tr].astype(np.float64)
        Pva, Yva = d["pred"][va].astype(np.float64), d["gt"][va].astype(np.float64)

        a, b = fit_calibration(Ptr, Ytr)
        mu = Ptr.mean(0)
        sd = Ptr.std(0)
        sd[sd == 0] = 1.0
        Ztr = (Ptr - mu) / sd

        curve = []
        best = None
        for alpha, W, c in ridge_path(Ztr, Ytr - Ptr, cfg["adapter"]["alphas"]):
            m = summarise(Yva, apply_adapter(Pva, mu, sd, W, c))
            m["alpha"] = alpha
            curve.append(m)
            score = (m["median_pcc"], -m["mean_rmse"])
            if best is None or score > best[0]:
                best = (score, alpha, W, c)
        _, alpha, W, c = best

        val = {"GHIST+": summarise(Yva, Pva),
               "calibration": summarise(Yva, apply_calibration(Pva, a, b)),
               "adapter": next(m for m in curve if m["alpha"] == alpha)}
        gate = val["adapter"]["median_pcc"] > val["calibration"]["median_pcc"]

        out = OUT_DIR / f"seed{seed}"
        out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"gene": genes, "a": a, "b": b}).to_csv(out / "calibration.csv", index=False)
        pd.DataFrame({"gene": genes, "mu": mu, "sd": sd}).to_csv(out / "adapter_scaler.csv", index=False)
        wdf = pd.DataFrame(np.vstack([c[None, :], W]), columns=genes, index=["(intercept)", *genes])
        wdf.index.name = "input \\ output"
        wdf.to_csv(out / "adapter_weights.csv.gz", float_format="%.6g")
        pd.DataFrame(curve).to_csv(out / "val_curve.csv", index=False)
        summary = {"seed": seed, "n_train": int(len(tr)), "n_val": int(len(va)), "selected_alpha": alpha,
                   "selected_on_alpha_grid_edge": alpha in (cfg["adapter"]["alphas"][0], cfg["adapter"]["alphas"][-1]),
                   "validation": val, "adapter_beats_calibration_on_validation": bool(gate),
                   "n_calibration_negative_slopes": int((b < 0).sum())}
        (out / "fit_summary.json").write_text(json.dumps(summary, indent=2))
        print(f"seed {seed}: alpha={alpha}  val median PCC  GHIST+={val['GHIST+']['median_pcc']:.4f}  "
              f"calib={val['calibration']['median_pcc']:.4f}  adapter={val['adapter']['median_pcc']:.4f}  "
              f"| mean RMSE {val['GHIST+']['mean_rmse']:.4f}/{val['calibration']['mean_rmse']:.4f}/{val['adapter']['mean_rmse']:.4f}"
              f"  gate={'PASS' if gate else 'FAIL'}")


if __name__ == "__main__":
    main()
