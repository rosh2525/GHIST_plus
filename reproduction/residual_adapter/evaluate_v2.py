"""v2 locked evaluation on the 20,000 test cells of each split (same test cells as v1).

Usage:  python evaluate_v2.py          (needs torch for the MLP)
Writes  outputs_v2/{summary.csv, summary_across_seeds.csv, summary.md}, outputs_v2/seed{s}/per_gene.csv.gz,
        outputs_v2/figures/*.png
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

from common import (load_cache, load_split, pcc_per_gene, positive_rmse_per_gene, rmse_per_gene, rows_for,
                    std_ratio_per_gene)
from train_adapter import apply_calibration
from v2_common import MODEL_DIR_V2, OUT_V2, block_index, features, load_config_v2, load_extras, make_mlp

METHODS = ["GHIST+", "calibration", "ridge", "mlp"]
LABEL = {"GHIST+": "GHIST+ (frozen)", "calibration": "Per-gene calibration", "ridge": "Ridge residual adapter",
         "mlp": "MLP residual adapter (MPS)"}
COLOUR = {"GHIST+": "#D1495B", "calibration": "#EDAE49", "ridge": "#4C8DB8", "mlp": "#2A9D8F"}


def predictors(cfg, seed, N, n_in, n_out):
    out = OUT_V2 / f"seed{seed}" / f"n{N}"
    cal = pd.read_csv(out / "calibration.csv")
    sc = pd.read_csv(out / "scaler.csv")
    w = pd.read_csv(out / "ridge_weights.csv.gz").to_numpy(np.float64)
    c, W = w[0], w[1:]
    mu, sd = sc.mu.to_numpy(np.float32), sc.sd.to_numpy(np.float32)
    m = cfg["methods"]["mlp"]
    net = make_mlp(n_in, n_out, m["hidden"], m["dropout"])
    net.load_state_dict(torch.load(MODEL_DIR_V2 / f"seed{seed}_n{N}_mlp.pt", map_location="cpu"))
    net.eval()

    def mlp(P, X):
        with torch.no_grad():
            return np.clip(P + net(torch.from_numpy((X - mu) / sd)).numpy(), 0, None)

    return {"GHIST+": lambda P, X: P,
            "calibration": lambda P, X: apply_calibration(P, cal.a.to_numpy(), cal.b.to_numpy()),
            "ridge": lambda P, X: np.clip(P + c + ((X - mu) / sd) @ W, 0, None),
            "mlp": mlp}


def per_gene(Y, Q, genes):
    return pd.DataFrame({"gene": genes, "pcc": pcc_per_gene(Y, Q), "rmse": rmse_per_gene(Y, Q),
                         "rmse_positive": positive_rmse_per_gene(Y, Q), "std_ratio": std_ratio_per_gene(Y, Q)})


def main() -> None:
    cfg = load_config_v2()
    d = load_cache()
    ex = load_extras(d["ids"], d["genes"])
    ids, genes = d["ids"], d["genes"]
    rows, tables = [], []
    for seed in (0, 1, 2):
        sp = load_split(seed)
        te = rows_for(ids, np.sort(sp.loc[sp.split == "test", "cell_id"].to_numpy()))
        P, X = features(d, ex, te)
        Y = d["gt"][te].astype(np.float32)
        base = per_gene(Y, P, genes)
        seed_tables = []
        for N in cfg["train_sizes"]:
            preds = predictors(cfg, seed, N, X.shape[1], len(genes))
            res = {mth: (base if mth == "GHIST+" else per_gene(Y, f(P, X), genes)) for mth, f in preds.items()}
            for mth, t in res.items():
                dp, dc = t.pcc - base.pcc, t.pcc - res["calibration"].pcc
                rows.append({"seed": seed, "n_train": N, "method": mth, "median_pcc": t.pcc.median(), "mean_pcc": t.pcc.mean(),
                             "mean_rmse": t.rmse.mean(), "median_rmse_positive": t.rmse_positive.median(),
                             "median_std_ratio": t.std_ratio.median(),
                             "genes_pcc_up_vs_ghist": int((dp > 0.005).sum()), "genes_pcc_down_vs_ghist": int((dp < -0.005).sum()),
                             "genes_pcc_up_vs_calib": int((dc > 0.005).sum()), "genes_pcc_down_vs_calib": int((dc < -0.005).sum())})
                seed_tables.append(t.assign(seed=seed, n_train=N, method=mth))
        pd.concat(seed_tables).to_csv(OUT_V2 / f"seed{seed}" / "per_gene.csv.gz", index=False, float_format="%.6g")
        tables += seed_tables
        print(f"seed {seed} evaluated", flush=True)

    summ = pd.DataFrame(rows)
    summ.to_csv(OUT_V2 / "summary.csv", index=False, float_format="%.5f")
    agg = (summ.groupby(["n_train", "method"], sort=False)
           .agg(median_pcc=("median_pcc", "mean"), median_pcc_sd=("median_pcc", "std"), mean_pcc=("mean_pcc", "mean"),
                mean_rmse=("mean_rmse", "mean"), rmse_pos=("median_rmse_positive", "mean"), sd_ratio=("median_std_ratio", "mean"),
                up=("genes_pcc_up_vs_ghist", "mean"), down=("genes_pcc_down_vs_ghist", "mean"))
           .reset_index())
    agg.to_csv(OUT_V2 / "summary_across_seeds.csv", index=False, float_format="%.5f")

    figs = OUT_V2 / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    # learning curves
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    for ax, (col, name) in zip(axes, [("median_pcc", "test median per-gene PCC ↑"), ("mean_rmse", "test mean RMSE ↓"),
                                      ("median_rmse_positive", "test median RMSE, measured>0 ↓")]):
        for mth in METHODS:
            g = summ[summ.method == mth].groupby("n_train")[col]
            ax.errorbar(g.mean().index, g.mean().values, yerr=g.std().values, marker="o", capsize=3, color=COLOUR[mth], label=LABEL[mth])
        ax.set_xscale("log"); ax.set_xlabel("local training cells (log scale)"); ax.set_ylabel(name)
    axes[0].legend(fontsize=8)
    fig.suptitle("Breast5K within-slide adaptation; mean ± SD over 3 spatial splits (20k test cells each)", fontsize=10)
    fig.tight_layout(); fig.savefig(figs / "learning_curve.png", dpi=150); plt.close(fig)

    # per-gene scatter at the pre-specified training size
    Nf = cfg["figure_train_size"]
    allg = pd.concat(tables)
    t0 = allg[(allg.seed == 0) & (allg.n_train == Nf)].pivot(index="gene", columns="method", values="pcc")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
    for ax, mth in zip(axes, ["ridge", "mlp"]):
        ax.scatter(t0["GHIST+"], t0[mth], s=8, color=COLOUR[mth])
        lim = [np.nanmin(t0.values), np.nanmax(t0.values)]
        ax.plot(lim, lim, "k--", lw=0.8)
        up, down = int(((t0[mth] - t0["GHIST+"]) > 0.005).sum()), int(((t0[mth] - t0["GHIST+"]) < -0.005).sum())
        ax.set_xlabel("GHIST+ PCC"); ax.set_ylabel(f"{LABEL[mth]} PCC"); ax.set_title(f"seed 0, N={Nf}: {up} genes ↑, {down} ↓", fontsize=10)
    fig.tight_layout(); fig.savefig(figs / f"pcc_scatter_seed0_n{Nf}.png", dpi=150); plt.close(fig)

    # spatial maps: all cells of the largest seed-0 test block
    sp = load_split(0)
    block_id = sp[sp.split == "test"].block.value_counts().idxmax()
    blk = block_index(d, 400.0)
    r = np.flatnonzero(blk[:, 0] * 100000 + blk[:, 1] == block_id)
    P, X = features(d, ex, r)
    Y = d["gt"][r]
    preds = predictors(cfg, 0, Nf, X.shape[1], len(genes))
    Q = {m: preds[m](P, X) for m in ["GHIST+", "ridge", "mlp"]}
    chosen = [g for g in cfg["map_genes"] if g in list(genes)]
    fig, axes = plt.subplots(len(chosen), 4, figsize=(13, 3.1 * len(chosen)), squeeze=False)
    for i, g in enumerate(chosen):
        j = list(genes).index(g)
        for k, col in enumerate(["measured", "GHIST+", "ridge", "mlp"]):
            v = Y[:, j] if col == "measured" else Q[col][:, j]
            vmax = max(np.percentile(v, 99), 1e-3)  # own scale per panel: compares spatial pattern, not magnitude
            ax = axes[i, k]
            ax.scatter(d["x"][r], d["y"][r], c=v, s=3, cmap="viridis", vmin=0, vmax=vmax, rasterized=True)
            ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
            ttl = "measured" if col == "measured" else f"{LABEL[col]}\nPCC={pcc_per_gene(Y[:, [j]], v[:, None])[0]:.2f}"
            ax.set_title(ttl, fontsize=9)
        axes[i, 0].set_ylabel(g, fontsize=11, fontweight="bold")
    fig.suptitle(f"log1p expression, largest seed-0 test block ({len(r)} cells), adapters trained on N={Nf}; colour scale per panel (99th pct)", fontsize=10)
    fig.tight_layout(); fig.savefig(figs / f"spatial_maps_seed0_n{Nf}.png", dpi=130); plt.close(fig)

    fits = [json.loads((OUT_V2 / f"seed{s}" / f"n{N}" / "fit.json").read_text()) for s in (0, 1, 2) for N in cfg["train_sizes"]]
    lines = ["# v2 results (auto-generated by evaluate_v2.py)", "",
             "Test = the same 20,000 buffered test cells per split as v1. Mean over 3 spatial splits (± SD for median PCC). log1p units.", "",
             "| local train cells | method | median PCC | mean PCC | mean RMSE | median RMSE (measured>0) | SD ratio | genes PCC ↑ / ↓ vs GHIST+ |",
             "|---|---|---|---|---|---|---|---|"]
    for _, a in agg.iterrows():
        lines.append(f"| {a.n_train:,} | {LABEL[a.method]} | {a.median_pcc:.4f} ± {a.median_pcc_sd:.4f} | {a.mean_pcc:.4f} | {a.mean_rmse:.4f} | "
                     f"{a.rmse_pos:.4f} | {a.sd_ratio:.3f} | {a.up:.0f} / {a.down:.0f} |")
    lines += ["", "Gene counts use a ±0.005 PCC threshold, averaged over seeds.", "",
              "| seed | N | ridge alpha | MLP best step / run | val median PCC GHIST+ / calib / ridge / MLP | fit seconds |", "|---|---|---|---|---|---|"]
    for f in fits:
        v = f["validation"]
        lines.append(f"| {f['seed']} | {f['n_train']:,} | {f['ridge_alpha']:g} | {f['mlp_best_step']} / {f['mlp_steps_run']} | "
                     f"{v['GHIST+']['median_pcc']:.4f} / {v['calibration']['median_pcc']:.4f} / {v['ridge']['median_pcc']:.4f} / {v['mlp']['median_pcc']:.4f} | {f['seconds']} |")
    (OUT_V2 / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
