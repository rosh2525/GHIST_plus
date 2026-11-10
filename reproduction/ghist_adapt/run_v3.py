"""Run the frozen v3 protocol on one or more slides.

Usage:  python run_v3.py breast5k [atera breast2] [--smoke]
Writes  outputs/<slide>/{runs.csv, per_gene.csv.gz, cores.json, slide_info.json}
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata

import adapt
from anchors import evaluation_mask, occupied_area_mm2, place_cores, prediction_pcs
from metrics import (block_stats, bootstrap_median_diff, build_ssim_grid, cmd_from_corr, corr_matrix, from_sums,
                     gene_sums, giotto_scores, spearman_columns, ssim_local, ssim_spatial_map, y_block_stats)
from slides import CACHE, load

HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "config_v3.json").read_text())
METHODS = ["GHIST+", "calibration", "RRC-authors", "ridge", "Adapt", "Adapt-DM"]
# authors' Figure3.ipynb D_MARKERS (cell-type marker programmes)
MARKERS = {
    "T": ["CD3E", "CD3G", "TRAC", "CD247"],
    "B": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD19"],
    "Plasma": ["MZB1", "TNFRSF17", "PRDM1", "TENT5C"],
    "Macrophage": ["C1QA", "C1QC", "APOC1", "CD68", "CD163", "MRC1", "MMP12"],
    "Myeloid": ["FCER1A", "CLEC9A", "LILRA4", "IL3RA", "CD14", "S100A8", "MNDA", "ITGAX"],
    "Stromal": ["LUM", "PDGFRA", "PDGFRB", "ACTA2", "MYH11", "DPT", "CCDC80", "POSTN", "LRRC15", "MMP2"],
    "Epithelial": ["EPCAM", "CDH1", "KRT8", "KRT7", "KRT5", "KRT14", "GATA3", "FOXA1", "AGR3", "PIGR"],
    "Malignant": ["ERBB2", "MKI67", "TOP2A", "CCND1", "APOBEC3B", "LYPD3", "S100A14", "TACSTD2"],
}
HEADLINE = (1.0, 4, "diverse", 0)


def programme_index(genes):
    gl = list(genes)
    return {k: [gl.index(g) for g in v if g in gl] for k, v in MARKERS.items() if sum(g in gl for g in v) >= 2}


def gene_groups(d, slide):
    path = CACHE / f"{slide}_giotto.npy"
    if path.exists():
        score = np.load(path)
    else:
        score = giotto_scores(d["y"].astype(np.float32), d["xy_um"].astype(np.float32))
        np.save(path, score)
    svg = np.argsort(-score, kind="stable")
    non = np.argsort(np.abs(score), kind="stable")
    return score, {"all": np.arange(len(score)), "svg20": svg[:20], "svg50": svg[:50], "nonsvg50": non[:50]}


def run_one(d, groups, tree, pcs, area, n_cores, placement, seed, tissue, feat_stats=None):
    a = CFG["anchors"]
    centres, r, aidx, agroups = place_cores(d["xy_um"], area, n_cores, placement, seed, P=d["pred"], pcs=pcs, tree=tree)
    emask = evaluation_mask(d["xy_um"], centres, r, a["evaluation_buffer_um"])
    eidx = np.flatnonzero(emask)
    t0 = time.time()
    model = adapt.fit(d["pred"][aidx], d["y"][aidx], agroups, CFG["cv"]["alphas"])
    t_fit = time.time() - t0
    rrc = adapt.fit_rrc_authors(d["pred"][aidx], d["y"][aidx], *feat_stats)

    Y = d["y"][eidx]
    P = d["pred"][eidx]
    xy = d["xy_um"][eidx]
    rng = np.random.default_rng(1000 + seed)
    sp_idx = np.sort(rng.choice(len(eidx), min(100_000, len(eidx)), replace=False))
    Y_sp_ranks = np.column_stack([rankdata(Y[sp_idx, g], method="average") for g in range(Y.shape[1])])
    grid = build_ssim_grid(xy[:, 0], xy[:, 1], width=400)
    ssim_genes = np.concatenate([groups["svg50"], groups["nonsvg50"]])
    gt_maps = {g: ssim_spatial_map(Y[:, g], grid) for g in ssim_genes}
    valid_cmd = Y.std(0) > 0
    cmd_sets = {"all": np.flatnonzero(valid_cmd), "svg50": np.intersect1d(groups["svg50"], np.flatnonzero(valid_cmd))}
    gt_corr = {k: corr_matrix(Y, v) for k, v in cmd_sets.items()}
    blk = np.floor(xy / 400.0).astype(np.int64)
    bkey = blk[:, 0] * 1_000_000 + blk[:, 1]
    order = np.argsort(bkey, kind="stable")
    starts = np.flatnonzero(np.r_[True, np.diff(bkey[order]) != 0])
    ystats = y_block_stats(Y, order, starts)

    base = {"slide": None, "core_area_mm2": area, "n_cores": n_cores, "placement": placement, "seed": seed,
            "total_core_area_mm2": area * n_cores, "core_fraction_of_tissue": area * n_cores / tissue,
            "n_anchor_cells": int(len(aidx)), "n_eval_cells": int(len(eidx)), "alpha": model["alpha"],
            "n_genes_ridge_branch": model["n_genes_ridge"], "fit_seconds": round(t_fit, 2)}
    progs = programme_index(d["genes"])
    y_sig = {k: Y[:, v].mean(1) for k, v in progs.items()}
    rows, genes_rows, stats = [], [], {}
    for meth in METHODS:
        Q = adapt.predict_rrc_authors(rrc, P) if meth == "RRC-authors" else adapt.predict(model, P, meth)
        gm = from_sums(gene_sums(Y, Q))
        sig = {}
        for k, v in progs.items():
            q = Q[:, v].mean(1)
            sig[k] = float(np.corrcoef(y_sig[k], q)[0, 1]) if q.std() > 0 and y_sig[k].std() > 0 else np.nan
        spear = spearman_columns(Y_sp_ranks, Q[sp_idx])
        ssim = np.full(Y.shape[1], np.nan)
        for g in ssim_genes:
            t = gt_maps[g]
            ssim[g] = ssim_local(t, ssim_spatial_map(Q[:, g], grid), float(np.ptp(t)))
        cmd = {k: cmd_from_corr(gt_corr[k], corr_matrix(Q, v)) for k, v in cmd_sets.items()}
        stats[meth] = block_stats(Y, Q, order, starts) if meth in ("GHIST+", "RRC-authors", "Adapt", "Adapt-DM") else None
        row = dict(base, method=meth,
                   **{f"median_pcc_{k}": float(np.nanmedian(gm["pcc"][v])) for k, v in groups.items()},
                   mean_pcc_all=float(np.nanmean(gm["pcc"])),
                   **{f"median_spearman_{k}": float(np.nanmedian(spear[v])) for k, v in groups.items() if k != "nonsvg50"},
                   mean_rmse=float(gm["rmse"].mean()), median_rmse_pos=float(np.nanmedian(gm["rmse_pos"])),
                   median_sd_ratio=float(np.nanmedian(gm["sd_ratio"])),
                   median_ssim_svg50=float(np.nanmedian(ssim[groups["svg50"]])),
                   median_ssim_nonsvg50=float(np.nanmedian(ssim[groups["nonsvg50"]])),
                   cmd_all=cmd["all"], cmd_svg50=cmd["svg50"],
                   median_programme_pcc=float(np.nanmedian(list(sig.values()))) if sig else np.nan,
                   **{f"programme_pcc_{k}": v for k, v in sig.items()})
        if meth == "Adapt":
            lo = np.clip(Q - model["r_adapt"][None, :], 0, None); hi = Q + model["r_adapt"][None, :]
            cov = ((Y >= lo) & (Y <= hi)).mean(0)
            row.update(conformal_coverage_median=float(np.median(cov)), conformal_width_mean=float((hi - lo).mean()))
            lo = np.clip(P - model["r_ghist"][None, :], 0, None); hi = P + model["r_ghist"][None, :]
            cov = ((Y >= lo) & (Y <= hi)).mean(0)
            row.update(ghist_conformal_coverage_median=float(np.median(cov)), ghist_conformal_width_mean=float((hi - lo).mean()))
        rows.append(row)
        genes_rows.append(pd.DataFrame({"method": meth, "gene_idx": np.arange(Y.shape[1]), "pcc": gm["pcc"],
                                        "rmse": gm["rmse"], "ssim": ssim, "use_ridge": model["use_ridge"]}))
        del Q
    gsets = {"all": groups["all"], "svg20": groups["svg20"]}
    for meth in ("Adapt", "Adapt-DM"):
        ci = bootstrap_median_diff(ystats, stats["GHIST+"], stats[meth], gsets, n_boot=500, seed=seed)
        r = next(x for x in rows if x["method"] == meth)
        for k, (lo, hi) in ci.items():
            r[f"boot_diff_{k}_lo"], r[f"boot_diff_{k}_hi"] = lo, hi
    for meth in ("Adapt", "Adapt-DM"):  # head-to-head with the authors' RRC
        ci = bootstrap_median_diff(ystats, stats["RRC-authors"], stats[meth], gsets, n_boot=500, seed=seed)
        r = next(x for x in rows if x["method"] == meth)
        for k, (lo, hi) in ci.items():
            r[f"boot_vs_rrc_{k}_lo"], r[f"boot_vs_rrc_{k}_hi"] = lo, hi
    rr = next(x for x in rows if x["method"] == "RRC-authors")
    rr["rrc_alpha"], rr["rrc_blend"] = rrc["alpha"], rrc["blend"]
    ghist = next(x for x in rows if x["method"] == "GHIST+")
    for r in rows:
        r["delta_median_pcc_all"] = r["median_pcc_all"] - ghist["median_pcc_all"]
    g = pd.concat(genes_rows)
    piv = g.pivot(index="gene_idx", columns="method", values="pcc")
    for r in rows:
        dlt = piv[r["method"]] - piv["GHIST+"]
        r["genes_pcc_up"] = int((dlt > 0.005).sum()); r["genes_pcc_down"] = int((dlt < -0.005).sum())
        dr = piv[r["method"]] - piv["RRC-authors"]
        r["genes_up_vs_rrc"] = int((dr > 0.005).sum()); r["genes_down_vs_rrc"] = int((dr < -0.005).sum())
    extra = {"centres_um": centres.tolist(), "radius_um": float(np.sqrt(area * 1e6 / np.pi))}
    if (area, n_cores, placement, seed) == HEADLINE:
        extra["weights"] = {"W": model["W"], "use_ridge": model["use_ridge"], "alpha": model["alpha"]}
    return rows, g, extra


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slides", nargs="+")
    ap.add_argument("--smoke", action="store_true", help="one config only (4 x 1 mm2, diverse, seed 0)")
    ap.add_argument("--resume", action="store_true", help="skip configs already in outputs/<slide>")
    ap.add_argument("--out", default="outputs", help="output folder (e.g. outputs_rrc for the head-to-head)")
    ap.add_argument("--headline-only", action="store_true", help="only the primary, paper-matched and random-4x1 configs")
    ap.add_argument("--shard", default=None, help="k/n: run every n-th remaining config into outputs/<slide>/shard_k (merge with merge_shards.py)")
    args = ap.parse_args()
    a = CFG["anchors"]
    for slide in args.slides:
        t0 = time.time()
        d = load(slide)
        score, groups = gene_groups(d, slide)
        tissue = occupied_area_mm2(d["xy_um"])
        tree = cKDTree(d["xy_um"])
        pcs = prediction_pcs(d["pred"])
        fs = d["pred"].std(0, dtype=np.float64); fs[fs == 0] = 1.0
        feat_stats = (d["pred"].mean(0, dtype=np.float64), fs)  # authors' RRC: standardise over all slide cells (RNA-free)
        main_out = HERE / ("outputs_smoke" if args.smoke else args.out) / slide
        out = main_out / f"shard_{args.shard.split('/')[0]}" if args.shard else main_out
        out.mkdir(parents=True, exist_ok=True)
        done_main = set(json.loads((main_out / "cores.json").read_text())["cores"]) if args.shard and (main_out / "cores.json").exists() else set()
        (out / "slide_info.json").write_text(json.dumps({"n_cells": int(len(d["ids"])), "n_genes": int(len(d["genes"])),
                                                         "tissue_area_mm2_100um_grid": tissue,
                                                         "genes": d["genes"].tolist(), "giotto_score": score.tolist()}, indent=1))
        configs = [(1.0, 4, "diverse", 0)] if args.smoke else [
            (ar, n, pl, s) for ar in a["core_area_mm2"] for n in a["n_cores"] for pl in a["placement"] for s in a["repeats"]]
        # run order only (no setting changes): primary headline config, then the paper-matched one, then the rest
        prio = {(1.0, 4, "diverse"): 0, (3.5, 1, "diverse"): 1, (1.0, 4, "random"): 2}
        configs.sort(key=lambda c: prio.get(c[:3], 3))
        if args.headline_only:
            configs = [c for c in configs if c[:3] in prio]
        if args.shard:
            k, nsh = map(int, args.shard.split("/"))
            todo = [c for c in configs if f"a{c[0]}_n{c[1]}_{c[2]}_s{c[3]}" not in done_main]
            configs = todo[k::nsh]
        rows, gene_tables, cores, skipped = [], [], {}, []
        if args.resume and (out / "runs.csv").exists() and (out / "cores.json").exists():
            rows = pd.read_csv(out / "runs.csv").to_dict("records")
            gene_tables = [pd.read_csv(out / "per_gene.csv.gz")]
            cores = json.loads((out / "cores.json").read_text())["cores"]
            print(f"[{slide}] resuming: {len(cores)} configs already done", flush=True)
        for ar, n, pl, s in configs:
            if ar * n > a["max_total_area_fraction"] * tissue:
                skipped.append({"core_area_mm2": ar, "n_cores": n, "reason": f"{ar*n:.1f} mm2 > 20% of {tissue:.1f} mm2 tissue"})
                continue
            key = f"a{ar}_n{n}_{pl}_s{s}"
            if key in cores:
                continue
            try:
                rr, gg, cc = run_one(d, groups, tree, pcs, ar, n, pl, s, tissue, feat_stats)
            except RuntimeError as e:
                skipped.append({"core_area_mm2": ar, "n_cores": n, "placement": pl, "seed": s, "reason": str(e)})
                continue
            for r in rr:
                r["slide"] = slide
            rows += rr
            gene_tables.append(gg.assign(config=key))
            wts = cc.pop("weights", None)
            if wts is not None:
                np.savez_compressed(out / "headline_weights.npz", W=wts["W"], use_ridge=wts["use_ridge"],
                                    alpha=wts["alpha"], genes=d["genes"])
            cores[key] = cc
            ad = next(x for x in rr if x["method"] == "Adapt")
            gh = next(x for x in rr if x["method"] == "GHIST+")
            rc = next(x for x in rr if x["method"] == "RRC-authors")
            print(f"[{slide}] {key}: anchors={ad['n_anchor_cells']} eval={ad['n_eval_cells']} alpha={ad['alpha']:g} "
                  f"ridge-genes={ad['n_genes_ridge_branch']} | medPCC GHIST+={gh['median_pcc_all']:.4f} "
                  f"RRC={rc['median_pcc_all']:.4f} Adapt={ad['median_pcc_all']:.4f} CI[{ad['boot_diff_all_lo']:+.4f},{ad['boot_diff_all_hi']:+.4f}] "
                  f"vsRRC[{ad['boot_vs_rrc_all_lo']:+.4f},{ad['boot_vs_rrc_all_hi']:+.4f}] "
                  f"RMSE {gh['mean_rmse']:.3f}->{ad['mean_rmse']:.3f} SSIM {gh['median_ssim_svg50']:.3f}->{ad['median_ssim_svg50']:.3f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
            pd.DataFrame(rows).to_csv(out / "runs.csv", index=False, float_format="%.6g")
            pd.concat(gene_tables).to_csv(out / "per_gene.csv.gz", index=False, float_format="%.5g")
            (out / "cores.json").write_text(json.dumps({"cores": cores, "skipped": skipped}, indent=1))
        (out / "cores.json").write_text(json.dumps({"cores": cores, "skipped": skipped}, indent=1))
        print(f"[{slide}] done: {len(rows) // len(METHODS)} runs, {len(skipped)} skipped, {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
