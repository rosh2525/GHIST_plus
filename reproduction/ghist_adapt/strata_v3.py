"""Where does GHIST+-Adapt help and where does it hurt? (exploratory, added after the v3 protocol freeze)

For the primary config (4 cores x 1 mm2, diverse, seeds 0-2) the identical cores/models are rebuilt and the
evaluation cells are split by cell type:
  breast5k, atera : GHIST+'s own H&E-based cell-type call (ct_pred in the released predictions) -> usable at
                    deployment without measured RNA
  breast2         : the authors' measured-RNA (scClassify) labels (GHIST+ export has no ct_pred)
Per stratum: per-cell profile PCC (across genes), per-gene PCC within the stratum (median), mean RMSE and
RMSE on measured-positive entries. Per gene programme: change in per-gene PCC over the whole evaluation set.

Usage:  python strata_v3.py
Writes  outputs/<slide>/strata.csv, outputs/<slide>/programme_delta.csv, outputs/strata.md, outputs/figures/where_it_helps.png
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import adapt
from anchors import evaluation_mask, place_cores, prediction_pcs
from metrics import from_sums, gene_sums
from run_v3 import CFG, HEADLINE, MARKERS
from slides import FETCH, SLIDES, load

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
NAMES = {"breast5k": "Breast5K", "atera": "ATERA", "breast2": "BreastCancer2"}
SHORT = {"Myeloid (excluding macrophage)": "Myeloid"}
METHODS = ["GHIST+", "Adapt", "Adapt-DM"]


def labels(slide, ids):
    if slide == "breast2":
        lab = pd.read_csv(FETCH / "evaluation_data/breast2/cell_type_filtered.csv").set_index("c_id")["ct"]
        return lab.reindex(ids).fillna("no label").map(lambda s: SHORT.get(s, s)).to_numpy(), "measured scClassify label"
    with np.load(SLIDES[slide]["pred"], allow_pickle=True) as z:
        pid, ct, cls = z["cell_ids"].astype(np.int64), z["ct_pred"], z["ct_classes"].astype(str)
    order = np.argsort(pid)
    rows = order[np.searchsorted(pid[order], ids)]
    return np.array([SHORT.get(c, c) for c in cls[ct[rows]]]), "GHIST+ H&E-predicted type"


def cell_pcc(Y, Q):
    y = Y - Y.mean(1, keepdims=True)
    q = Q - Q.mean(1, keepdims=True)
    den = np.sqrt((y * y).sum(1) * (q * q).sum(1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, (y * q).sum(1) / den, np.nan)


def main():
    import sys
    want = [x for x in sys.argv[1:]] or ["breast5k", "atera", "breast2"]
    a, n, pl, _ = HEADLINE
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    if "--figure-only" in want:
        want = [w for w in want if w != "--figure-only"] or ["breast5k", "atera", "breast2"]
        return figure([pd.read_csv(OUT / s / "strata.csv") for s in want], want)
    all_rows, prog_rows, md = [], [], ["# Where GHIST+-Adapt helps and where it hurts (exploratory)", "",
                                       "Primary config (4 cores x 1 mm², diverse), mean of 3 repeats. Δ = Adapt − GHIST+, signed so that **positive = better**.", ""]
    for slide in want:
        d = load(slide)
        lab, lab_source = labels(slide, d["ids"])
        tree = cKDTree(d["xy_um"])
        pcs = prediction_pcs(d["pred"])
        rows = []
        for seed in CFG["anchors"]["repeats"]:
            centres, r, aidx, groups = place_cores(d["xy_um"], a, n, pl, seed, P=d["pred"], pcs=pcs, tree=tree)
            eidx = np.flatnonzero(evaluation_mask(d["xy_um"], centres, r, CFG["anchors"]["evaluation_buffer_um"]))
            model = adapt.fit(d["pred"][aidx], d["y"][aidx], groups, CFG["cv"]["alphas"])
            Y, P, L = d["y"][eidx], d["pred"][eidx], lab[eidx]
            preds = {m: adapt.predict(model, P, m) for m in METHODS}
            for ct in np.unique(L):
                sel = L == ct
                if sel.sum() < 1000:
                    continue
                for m, Q in preds.items():
                    gm = from_sums(gene_sums(Y[sel], Q[sel]))
                    cp = cell_pcc(Y[sel].astype(np.float64), Q[sel].astype(np.float64))
                    rows.append({"slide": slide, "seed": seed, "cell_type": ct, "n_cells": int(sel.sum()), "method": m,
                                 "cell_profile_pcc": float(np.nanmean(cp)), "gene_pcc_within": float(np.nanmedian(gm["pcc"])),
                                 "rmse": float(gm["rmse"].mean()), "rmse_pos": float(np.nanmedian(gm["rmse_pos"]))})
            del preds
        st = pd.DataFrame(rows)
        st.to_csv(OUT / slide / "strata.csv", index=False, float_format="%.5f")
        all_rows.append(st)
        # gene programmes: per-gene PCC change over the whole evaluation set (from run_v3 per-gene output)
        pg = pd.read_csv(OUT / slide / "per_gene.csv.gz")
        pg = pg[pg.config.str.startswith(f"a{a}_n{n}_{pl}_")]
        genes = json.loads((OUT / slide / "slide_info.json").read_text())["genes"]
        piv = pg.pivot_table(index="gene_idx", columns="method", values="pcc", aggfunc="mean")
        delta = (piv["Adapt"] - piv["GHIST+"]).rename("delta_pcc")
        gtab = pd.DataFrame({"gene": [genes[i] for i in piv.index], "ghist_pcc": piv["GHIST+"].values,
                             "adapt_pcc": piv["Adapt"].values, "delta_pcc": delta.values})
        gtab.to_csv(OUT / slide / "gene_delta.csv", index=False, float_format="%.5f")
        for prog, mk in MARKERS.items():
            t = gtab[gtab.gene.isin(mk)]
            if len(t):
                prog_rows.append({"slide": slide, "programme": prog, "n_genes": len(t), "mean_delta_pcc": t.delta_pcc.mean(),
                                  "up": int((t.delta_pcc > 0.005).sum()), "down": int((t.delta_pcc < -0.005).sum()),
                                  "genes": ", ".join(f"{g} {v:+.3f}" for g, v in zip(t.gene, t.delta_pcc))})
        md += [f"## {NAMES[slide]} — strata by {lab_source}", ""]
        mean = st.groupby(["cell_type", "method"]).mean(numeric_only=True).reset_index()
        md += ["| cell type | cells | profile PCC GHIST+ → Adapt | gene PCC within type | RMSE | RMSE (measured>0) | verdict |", "|---|---|---|---|---|---|---|"]
        for ct in mean.cell_type.unique():
            g = mean[(mean.cell_type == ct) & (mean.method == "GHIST+")].iloc[0]
            ad = mean[(mean.cell_type == ct) & (mean.method == "Adapt")].iloc[0]
            dm = mean[(mean.cell_type == ct) & (mean.method == "Adapt-DM")].iloc[0]
            wins = [ad.cell_profile_pcc > g.cell_profile_pcc, ad.gene_pcc_within > g.gene_pcc_within, ad.rmse < g.rmse, ad.rmse_pos < g.rmse_pos]
            verdict = "use Adapt" if sum(wins[:2]) == 2 else ("mixed" if sum(wins[:2]) == 1 else "keep GHIST+")
            md.append(f"| {ct} | {int(g.n_cells):,} | {g.cell_profile_pcc:.3f} → {ad.cell_profile_pcc:.3f} | {g.gene_pcc_within:.3f} → {ad.gene_pcc_within:.3f} | "
                      f"{g.rmse:.3f} → {ad.rmse:.3f} | {g.rmse_pos:.3f} → {ad.rmse_pos:.3f} (DM {dm.rmse_pos:.3f}) | {verdict} |")
        worst = gtab.nsmallest(5, "delta_pcc")
        best = gtab.nlargest(5, "delta_pcc")
        md += ["", "Most improved genes: " + ", ".join(f"{r.gene} ({r.ghist_pcc:.2f}→{r.adapt_pcc:.2f})" for r in best.itertuples()),
               "", "Most worsened genes: " + ", ".join(f"{r.gene} ({r.ghist_pcc:.2f}→{r.adapt_pcc:.2f})" for r in worst.itertuples()), ""]
    pr = pd.DataFrame(prog_rows)
    pr.to_csv(OUT / "programme_delta.csv", index=False, float_format="%.5f")
    md += ["## Marker-gene programmes (per-gene PCC change, whole evaluation set)", "",
           "| slide | programme | genes ↑ / ↓ | mean ΔPCC | per gene |", "|---|---|---|---|---|"]
    md += [f"| {NAMES[r.slide]} | {r.programme} | {r.up} / {r.down} | {r.mean_delta_pcc:+.3f} | {r.genes} |" for r in pr.itertuples()]
    (OUT / "strata.md").write_text("\n".join(md) + "\n")

    figure(all_rows, want)
    print("\n".join(md))


def figure(all_rows, want):
    # figure: Δ heatmap per slide (positive = better)
    st = pd.concat(all_rows)
    mean = st.groupby(["slide", "cell_type", "method"]).mean(numeric_only=True).reset_index()
    mets = [("cell_profile_pcc", "cell profile PCC", 1), ("gene_pcc_within", "gene PCC within type", 1),
            ("rmse", "RMSE", -1), ("rmse_pos", "RMSE, expressing cells", -1)]
    fig, axes = plt.subplots(1, len(want), figsize=(5 * len(want), 4.8), squeeze=False)
    for j, slide in enumerate(want):
        t = mean[mean.slide == slide]
        cts = [c for c in t.cell_type.unique() if c not in ("Unassigned", "no label")]
        M = np.zeros((len(cts), len(mets)))
        for i, ct in enumerate(cts):
            g = t[(t.cell_type == ct) & (t.method == "GHIST+")].iloc[0]
            ad = t[(t.cell_type == ct) & (t.method == "Adapt")].iloc[0]
            for k, (m, _, s) in enumerate(mets):
                M[i, k] = 100 * s * (ad[m] - g[m]) / abs(g[m])
        lim = np.nanpercentile(np.abs(M), 95)
        ax = axes[0, j]
        im = ax.imshow(np.clip(M, -lim, lim), cmap="RdBu", vmin=-lim, vmax=lim, aspect="auto")
        for (i, k), v in np.ndenumerate(M):
            ax.text(k, i, f"{v:+.0f}%", ha="center", va="center", fontsize=7)
        ax.set_xticks(range(len(mets)), [x[1] for x in mets], rotation=30, ha="right", fontsize=8)
        ax.set_yticks(range(len(cts)), cts, fontsize=8)
        ax.set_title(f"{NAMES[slide]} ({'measured labels' if slide == 'breast2' else 'GHIST+ H&E cell-type call'})", fontsize=9)
        plt.colorbar(im, ax=ax, label="% change vs GHIST+ (blue = Adapt better)")
    fig.suptitle("Where GHIST+-Adapt helps (blue) and hurts (red): primary config, mean of 3 repeats", fontsize=10)
    fig.tight_layout(); fig.savefig(OUT / "figures" / "where_it_helps.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    main()
