"""Locked evaluation of GHIST+ / per-gene calibration / multigene residual adapter.

Primary set: 'test' (20,000 cells in buffered held-out blocks). Secondary: 'rest'
(every cell in blocks at least one block away from train/validation).

Usage:  python evaluate.py
Writes  outputs/summary.csv, outputs/summary.md, outputs/seed{s}/per_gene.csv, outputs/figures/*.png
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (OUT_DIR, load_cache, load_config, load_split, pcc_per_gene, positive_rmse_per_gene,
                    rmse_per_gene, rows_for, std_ratio_per_gene)
from train_adapter import apply_adapter, apply_calibration

METHODS = ["GHIST+", "calibration", "adapter"]
LABEL = {"GHIST+": "GHIST+ (frozen)", "calibration": "Per-gene calibration", "adapter": "Multigene residual adapter"}
COLOUR = {"GHIST+": "#D1495B", "calibration": "#EDAE49", "adapter": "#2A9D8F"}


def load_models(seed, genes):
    out = OUT_DIR / f"seed{seed}"
    cal = pd.read_csv(out / "calibration.csv").set_index("gene").loc[genes]
    sc = pd.read_csv(out / "adapter_scaler.csv").set_index("gene").loc[genes]
    w = pd.read_csv(out / "adapter_weights.csv.gz", index_col=0)
    assert list(w.columns) == list(genes) and list(w.index[1:]) == list(genes)
    W, c = w.iloc[1:].to_numpy(np.float64), w.iloc[0].to_numpy(np.float64)
    return {"GHIST+": lambda P: P,
            "calibration": lambda P: apply_calibration(P, cal.a.to_numpy(), cal.b.to_numpy()),
            "adapter": lambda P: apply_adapter(P, sc.mu.to_numpy(), sc.sd.to_numpy(), W, c)}


def evaluate_rows(d, rows, models, chunk=100000):
    """Per-gene metrics, computed via sufficient statistics so 'rest' (~500k cells) fits in memory."""
    G = len(d["genes"])
    acc = {m: dict(n=0, sy=np.zeros(G), sq=np.zeros(G), syy=np.zeros(G), sqq=np.zeros(G), syq=np.zeros(G),
                   se=np.zeros(G), npos=np.zeros(G), sepos=np.zeros(G)) for m in models}
    for s in range(0, len(rows), chunk):
        r = rows[s:s + chunk]
        Y, P = d["gt"][r].astype(np.float64), d["pred"][r].astype(np.float64)
        pos = Y > 0
        for m, f in models.items():
            Q, a = f(P), acc[m]
            a["n"] += len(r)
            a["sy"] += Y.sum(0); a["sq"] += Q.sum(0)
            a["syy"] += (Y * Y).sum(0); a["sqq"] += (Q * Q).sum(0); a["syq"] += (Y * Q).sum(0)
            e2 = (Y - Q) ** 2
            a["se"] += e2.sum(0); a["npos"] += pos.sum(0); a["sepos"] += np.where(pos, e2, 0).sum(0)
    res = {}
    for m, a in acc.items():
        n = a["n"]
        vy = a["syy"] / n - (a["sy"] / n) ** 2
        vq = a["sqq"] / n - (a["sq"] / n) ** 2
        cov = a["syq"] / n - a["sy"] * a["sq"] / n ** 2
        with np.errstate(invalid="ignore", divide="ignore"):
            pcc = np.where((vy > 1e-12) & (vq > 1e-12), cov / np.sqrt(np.clip(vy, 0, None) * np.clip(vq, 0, None)), np.nan)
            res[m] = pd.DataFrame({"gene": d["genes"], "pcc": pcc, "rmse": np.sqrt(a["se"] / n),
                                   "rmse_positive": np.where(a["npos"] > 0, np.sqrt(a["sepos"] / np.maximum(a["npos"], 1)), np.nan),
                                   "std_ratio": np.where(vy > 1e-12, np.sqrt(np.clip(vq, 0, None) / vy), np.nan),
                                   "mean_pred": a["sq"] / n, "mean_measured": a["sy"] / n})
    return res


def summary_rows(seed, split, res):
    base, cal = res["GHIST+"], res["calibration"]
    out = []
    for m, t in res.items():
        dp, dc = t.pcc - base.pcc, t.pcc - cal.pcc
        out.append({"seed": seed, "set": split, "method": m, "n_genes_defined": int(t.pcc.notna().sum()),
                    "median_pcc": t.pcc.median(), "mean_pcc": t.pcc.mean(), "median_rmse": t.rmse.median(),
                    "mean_rmse": t.rmse.mean(), "median_rmse_positive": t.rmse_positive.median(),
                    "median_std_ratio": t.std_ratio.median(),
                    "genes_pcc_up_vs_ghist": int((dp > 0.005).sum()), "genes_pcc_down_vs_ghist": int((dp < -0.005).sum()),
                    "genes_pcc_up_vs_calib": int((dc > 0.005).sum()), "genes_pcc_down_vs_calib": int((dc < -0.005).sum()),
                    "genes_rmse_down_vs_ghist": int((t.rmse < base.rmse - 1e-3).sum()),
                    "genes_rmse_up_vs_ghist": int((t.rmse > base.rmse + 1e-3).sum())})
    return out


def plot_layout(d, cfg, path):
    seeds = cfg["split"]["seeds"]
    fig, axes = plt.subplots(1, len(seeds), figsize=(5 * len(seeds), 5), squeeze=False)
    sub = np.random.default_rng(0).choice(len(d["ids"]), 60000, replace=False)
    for ax, seed in zip(axes[0], seeds):
        sp = load_split(seed)
        ax.scatter(d["x"][sub], d["y"][sub], s=0.2, c="#dddddd", rasterized=True)
        for name, col in (("rest", "#9ecae1"), ("train", "#d62728"), ("val", "#ff7f0e"), ("test", "#2ca02c")):
            r = rows_for(d["ids"], np.sort(sp.loc[sp.split == name, "cell_id"].to_numpy()))
            if name == "rest":
                r = r[::10]
            ax.scatter(d["x"][r], d["y"][r], s=0.3, c=col, label=name, rasterized=True)
        ax.set_title(f"split seed {seed}"); ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
    axes[0, -1].legend(markerscale=15, loc="upper right", fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def plot_pcc(per_gene, path):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    t = per_gene[per_gene.set == "test"]
    data = [t[t.method == m].pcc.dropna() for m in METHODS]
    bp = axes[0].boxplot(data, patch_artist=True, showfliers=False)
    for patch, m in zip(bp["boxes"], METHODS):
        patch.set_facecolor(COLOUR[m])
    axes[0].set_xticks(range(1, 4), [LABEL[m] for m in METHODS], rotation=15)
    axes[0].set_ylabel("per-gene PCC (test, seed 0)")
    piv = t.pivot(index="gene", columns="method", values="pcc")
    axes[1].scatter(piv["GHIST+"], piv["adapter"], s=8, c=COLOUR["adapter"])
    lim = [np.nanmin(piv.values), np.nanmax(piv.values)]
    axes[1].plot(lim, lim, "k--", lw=0.8)
    axes[1].set_xlabel("GHIST+ PCC"); axes[1].set_ylabel("adapter PCC")
    axes[2].hist((piv["adapter"] - piv["GHIST+"]).dropna(), bins=40, color=COLOUR["adapter"])
    axes[2].axvline(0, c="k", lw=0.8); axes[2].set_xlabel("ΔPCC (adapter − GHIST+)"); axes[2].set_ylabel("genes")
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def plot_maps(d, cfg, models, path):
    genes = list(d["genes"])
    chosen = [g for g in cfg["map_genes"]["candidates"] if g in genes][:6]
    sp = load_split(cfg["split"]["seeds"][0])
    test = sp[sp.split == "test"]
    block = test.block.value_counts().idxmax()
    # all cells in the largest seed-0 test block
    bs = cfg["split"]["block_um"]
    bx = np.floor((d["x"] - d["x"].min()) / bs).astype(int); by = np.floor((d["y"] - d["y"].min()) / bs).astype(int)
    r = np.flatnonzero(bx * 100000 + by == block)
    Y, P = d["gt"][r], d["pred"][r].astype(np.float64)
    preds = {m: f(P) for m, f in models.items()}
    cols = ["measured", *METHODS]
    fig, axes = plt.subplots(len(chosen), 4, figsize=(13, 3.1 * len(chosen)), squeeze=False)
    for i, g in enumerate(chosen):
        j = genes.index(g)
        vmax = max(np.percentile(Y[:, j], 99), 1e-3)
        for k, col in enumerate(cols):
            v = Y[:, j] if col == "measured" else preds[col][:, j]
            ax = axes[i, k]
            ax.scatter(d["x"][r], d["y"][r], c=v, s=3, cmap="viridis", vmin=0, vmax=vmax, rasterized=True)
            ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
            title = col if col == "measured" else f"{LABEL[col]}\nPCC={pcc_per_gene(Y[:, [j]], preds[col][:, [j]])[0]:.2f}"
            ax.set_title(title, fontsize=9)
        axes[i, 0].set_ylabel(g, fontsize=11, fontweight="bold")
    fig.suptitle(f"log1p expression, largest seed-0 test block ({len(r)} cells); shared colour scale per gene", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)
    return chosen, int(len(r))


def main() -> None:
    cfg = load_config()
    d = load_cache()
    ids, genes = d["ids"], d["genes"]
    rows, all_gene = [], []
    for seed in cfg["split"]["seeds"]:
        sp = load_split(seed)
        models = load_models(seed, genes)
        for split in ("test", "rest"):
            r = rows_for(ids, np.sort(sp.loc[sp.split == split, "cell_id"].to_numpy()))
            res = evaluate_rows(d, r, models)
            rows += summary_rows(seed, split, res)
            for m, t in res.items():
                all_gene.append(t.assign(seed=seed, set=split, method=m, n_cells=len(r)))
        pg = pd.concat([t for t in all_gene if t.seed.iloc[0] == seed])
        pg.to_csv(OUT_DIR / f"seed{seed}" / "per_gene.csv.gz", index=False, float_format="%.6g")
        print(f"seed {seed} evaluated")
    summ = pd.DataFrame(rows)
    summ.to_csv(OUT_DIR / "summary.csv", index=False, float_format="%.5f")

    agg = (summ.groupby(["set", "method"], sort=False)
           .agg(median_pcc_mean=("median_pcc", "mean"), median_pcc_sd=("median_pcc", "std"),
                mean_rmse_mean=("mean_rmse", "mean"), mean_rmse_sd=("mean_rmse", "std"),
                median_rmse_positive=("median_rmse_positive", "mean"), median_std_ratio=("median_std_ratio", "mean"),
                genes_up_vs_ghist=("genes_pcc_up_vs_ghist", "mean"), genes_down_vs_ghist=("genes_pcc_down_vs_ghist", "mean"),
                genes_up_vs_calib=("genes_pcc_up_vs_calib", "mean"), genes_down_vs_calib=("genes_pcc_down_vs_calib", "mean"))
           .reset_index())
    agg.to_csv(OUT_DIR / "summary_across_seeds.csv", index=False, float_format="%.5f")

    figs = OUT_DIR / "figures"
    figs.mkdir(exist_ok=True)
    seed0 = cfg["split"]["seeds"][0]
    plot_layout(d, cfg, figs / "split_layout.png")
    plot_pcc(pd.concat([t for t in all_gene if t.seed.iloc[0] == seed0]), figs / "pcc_test_seed0.png")
    chosen, nblock = plot_maps(d, cfg, load_models(seed0, genes), figs / "spatial_maps_seed0.png")

    fits = {s: json.loads((OUT_DIR / f"seed{s}" / "fit_summary.json").read_text()) for s in cfg["split"]["seeds"]}
    lines = ["# Results (auto-generated by evaluate.py)", "",
             "Per-gene metrics over 258 shared genes; log1p units. `±` = SD over the three spatial splits.", "",
             "| set | method | median PCC | mean RMSE | median RMSE (measured>0) | median SD ratio (pred/measured) | genes PCC ↑ / ↓ vs GHIST+ | genes PCC ↑ / ↓ vs calibration |",
             "|---|---|---|---|---|---|---|---|"]
    for _, a in agg.iterrows():
        lines.append(f"| {a.set} | {LABEL[a.method]} | {a.median_pcc_mean:.4f} ± {a.median_pcc_sd:.4f} | "
                     f"{a.mean_rmse_mean:.4f} ± {a.mean_rmse_sd:.4f} | {a.median_rmse_positive:.4f} | {a.median_std_ratio:.3f} | "
                     f"{a.genes_up_vs_ghist:.0f} / {a.genes_down_vs_ghist:.0f} | {a.genes_up_vs_calib:.0f} / {a.genes_down_vs_calib:.0f} |")
    lines += ["", "Gene counts use a ±0.005 PCC threshold, averaged over seeds.", "",
              "| seed | selected alpha | val median PCC GHIST+ / calib / adapter | validation gate (adapter > calib) |", "|---|---|---|---|"]
    for s, f in fits.items():
        v = f["validation"]
        lines.append(f"| {s} | {f['selected_alpha']:g} | {v['GHIST+']['median_pcc']:.4f} / {v['calibration']['median_pcc']:.4f} / "
                     f"{v['adapter']['median_pcc']:.4f} | {'pass' if f['adapter_beats_calibration_on_validation'] else 'fail'} |")
    lines += ["", f"Spatial-map genes (fixed in config.json before testing): {', '.join(chosen)}; block of {nblock} cells."]
    (OUT_DIR / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
