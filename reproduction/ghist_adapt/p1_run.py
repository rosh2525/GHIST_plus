"""P0 + P1 (config_v4.json): soft hurdle expectation and hurdle-aware core correction vs the authors' RRC.

Usage:  python p1_run.py breast5k [atera] [--smoke]
Writes  outputs_v4/<slide>/{runs.csv, per_gene.csv.gz, zero_shot.json}
"""
from __future__ import annotations

import argparse
import gc
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata

import adapt
from anchors import evaluation_mask, occupied_area_mm2, place_cores, prediction_pcs
from hurdle import load_hurdle
from metrics import (block_stats, bootstrap_median_diff, build_ssim_grid, cmd_from_corr, corr_matrix, from_sums,
                     gene_sums, pcc_columns, spearman_columns, ssim_local, ssim_spatial_map, y_block_stats)
from run_v3 import CFG, gene_groups
from slides import load
from one_shot_rrc import fit_rrc as rrc_fit

HERE = Path(__file__).resolve().parent
METHODS = ["GHIST+", "GHIST+-Soft", "RRC-authors", "RRC-soft", "HAC"]
CONFIGS = [(1.0, 4, "diverse"), (3.5, 1, "diverse"), (1.0, 4, "random")]
ALPHAS = [1, 10, 100, 1000, 10000, 100000]


def slide_stats(X):
    m = X.mean(0, dtype=np.float64); s = X.std(0, dtype=np.float64); s[s == 0] = 1.0
    return m.astype(np.float32), s.astype(np.float32)


# ------------------------------------------------------------------ HAC
def _mag_fit(Y, R):
    """Per gene y ~ a + b*r on positive cells; fallback a = mean positive y, b = 0."""
    G = Y.shape[1]
    a = np.zeros(G); b = np.zeros(G)
    allpos = Y[Y > 0]
    glob = float(allpos.mean()) if allpos.size else 0.5
    for g in range(G):
        pos = Y[:, g] > 0
        n = int(pos.sum())
        if n == 0:
            a[g] = glob; continue
        y, r = Y[pos, g].astype(np.float64), R[pos, g].astype(np.float64)
        if n >= 20 and r.var() > 0:
            b[g] = np.cov(r, y, bias=True)[0, 1] / r.var(); a[g] = y.mean() - b[g] * r.mean()
        else:
            a[g] = y.mean()
    return a.astype(np.float32), b.astype(np.float32)


def _hac_predict(m, X, q, R, alpha_key="W", chunk=100_000):
    out = np.empty(q.shape, dtype=np.float32)
    for s in range(0, len(q), chunk):
        x = (X[s:s + chunk] - m["xm"]) @ m[alpha_key] + m["c"]
        qh = np.clip(q[s:s + chunk] + x, 0.001, 0.999)
        mh = np.maximum(m["a"] + m["b"] * R[s:s + chunk], 0.01)
        out[s:s + chunk] = qh * mh
    return out


def hac_fit(X, q, R, Y, groups):
    """X: standardised [S, logit q] (n x 2G); q, R=relu(mu), Y: n x G core cells."""
    T = (Y > 0).astype(np.float32) - q
    folds = np.unique(groups)
    oof = np.zeros((len(ALPHAS), *Y.shape), dtype=np.float32)
    for f in folds:
        te, tr = groups == f, groups != f
        xm = X[tr].mean(0); tm = T[tr].mean(0)
        U, S_, Vt = np.linalg.svd((X[tr] - xm).astype(np.float64), full_matrices=False)
        UtT = U.T @ (T[tr] - tm)
        a_, b_ = _mag_fit(Y[tr], R[tr])
        for k, al in enumerate(ALPHAS):
            W = (Vt.T @ ((S_ / (S_ ** 2 + al))[:, None] * UtT)).astype(np.float32)
            m = {"xm": xm, "W": W, "c": tm, "a": a_, "b": b_}
            oof[k, te] = _hac_predict(m, X[te], q[te], R[te])
    scores = [np.nanmedian(pcc_columns(Y, oof[k])) for k in range(len(ALPHAS))]
    al = ALPHAS[int(np.nanargmax(scores))]
    xm = X.mean(0); tm = T.mean(0)
    U, S_, Vt = np.linalg.svd((X - xm).astype(np.float64), full_matrices=False)
    W = (Vt.T @ ((S_ / (S_ ** 2 + al))[:, None] * (U.T @ (T - tm)))).astype(np.float32)
    a_, b_ = _mag_fit(Y, R)
    return {"xm": xm, "W": W, "c": tm, "a": a_, "b": b_, "alpha": al, "cv_scores": scores,
            "oof": oof[int(np.nanargmax(scores))]}  # out-of-fold predictions at the chosen alpha (used by P3)


# ------------------------------------------------------------------ P0 diagnostic
def rrc_loco_score(P, Y, groups, feat, r):
    """Leave-one-core-out mean PCC of the authors' RRC with its chosen alpha/blend."""
    x = ((P - feat[0]) / feat[1]).astype(np.float32)
    oof = np.zeros_like(Y, dtype=np.float32)
    for f in np.unique(groups):
        te, tr = groups == f, groups != f
        xm, ym, coef = rrc_fit(x[tr], (Y[tr] - P[tr]).astype(np.float32), r["alpha"])
        oof[te] = np.clip(P[te] + r["blend"] * ((x[te] - xm) @ coef + ym), 0, None)
    return float(np.nanmean(pcc_columns(Y, oof)))


def evaluate(Y, preds, groups, xy, seed, refs=("GHIST+", "RRC-authors")):
    rng = np.random.default_rng(1000 + seed)
    sp = np.sort(rng.choice(len(Y), min(100_000, len(Y)), replace=False))
    Yr = np.column_stack([rankdata(Y[sp, g], method="average") for g in range(Y.shape[1])])
    grid = build_ssim_grid(xy[:, 0], xy[:, 1], width=400)
    sg = groups["svg50"]
    gt_maps = {g: ssim_spatial_map(Y[:, g], grid) for g in sg}
    valid = np.flatnonzero(Y.std(0) > 0)
    gt_corr = corr_matrix(Y, valid)
    blk = np.floor(xy / 400.0).astype(np.int64)
    key = blk[:, 0] * 1_000_000 + blk[:, 1]
    order = np.argsort(key, kind="stable")
    starts = np.flatnonzero(np.r_[True, np.diff(key[order]) != 0])
    ystats = y_block_stats(Y, order, starts)
    rows, pg, stats = {}, {}, {}
    for meth, Q in preds.items():
        gm = from_sums(gene_sums(Y, Q))
        spear = spearman_columns(Yr, Q[sp])
        ssim = [ssim_local(gt_maps[g], ssim_spatial_map(Q[:, g], grid), float(np.ptp(gt_maps[g]))) for g in sg]
        rows[meth] = {f"median_pcc_{k}": float(np.nanmedian(gm["pcc"][v])) for k, v in groups.items()}
        rows[meth].update(mean_pcc_all=float(np.nanmean(gm["pcc"])), median_spearman_all=float(np.nanmedian(spear)),
                          mean_rmse=float(gm["rmse"].mean()), median_rmse_pos=float(np.nanmedian(gm["rmse_pos"])),
                          median_sd_ratio=float(np.nanmedian(gm["sd_ratio"])), median_ssim_svg50=float(np.nanmedian(ssim)),
                          cmd_all=cmd_from_corr(gt_corr, corr_matrix(Q, valid)))
        pg[meth] = gm["pcc"]
        stats[meth] = block_stats(Y, Q, order, starts)
    gsets = {"all": groups["all"], "svg20": groups["svg20"]}
    for meth in preds:
        for ref in refs:
            if meth == ref:
                continue
            ci = bootstrap_median_diff(ystats, stats[ref], stats[meth], gsets, n_boot=500, seed=seed)
            for k, (lo, hi) in ci.items():
                rows[meth][f"ci_vs_{ref}_{k}_lo"], rows[meth][f"ci_vs_{ref}_{k}_hi"] = lo, hi
        for ref in refs:
            dlt = pg[meth] - pg[ref]
            rows[meth][f"genes_up_vs_{ref}"] = int((dlt > 0.005).sum()); rows[meth][f"genes_down_vs_{ref}"] = int((dlt < -0.005).sum())
    return rows, pg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slides", nargs="+")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    a = CFG["anchors"]
    for slide in args.slides:
        t0 = time.time()
        d = load(slide)
        h = load_hurdle(slide, d["ids"])
        q, R = h["q"], np.maximum(h.pop("mu"), 0)
        S = (q * R).astype(np.float32)
        Y, P = d["y"], d["pred"]
        _, groups = gene_groups(d, slide)
        out = HERE / ("outputs_v4_smoke" if args.smoke else "outputs_v4") / slide
        out.mkdir(parents=True, exist_ok=True)
        # zero-shot on all cells
        zs = {}
        for name, Q in (("GHIST+", P), ("GHIST+-Soft", S)):
            gm = from_sums(gene_sums(Y, Q))
            zs[name] = {f"median_pcc_{k}": float(np.nanmedian(gm["pcc"][v])) for k, v in groups.items()}
            zs[name].update(mean_pcc_all=float(np.nanmean(gm["pcc"])), mean_rmse=float(gm["rmse"].mean()),
                            median_rmse_pos=float(np.nanmedian(gm["rmse_pos"])), median_sd_ratio=float(np.nanmedian(gm["sd_ratio"])))
        (out / "zero_shot.json").write_text(json.dumps(zs, indent=2))
        print(f"[{slide}] zero-shot all cells: median PCC export {zs['GHIST+']['median_pcc_all']:.4f} -> soft "
              f"{zs['GHIST+-Soft']['median_pcc_all']:.4f}; top-20 SVG {zs['GHIST+']['median_pcc_svg20']:.3f} -> "
              f"{zs['GHIST+-Soft']['median_pcc_svg20']:.3f}", flush=True)
        # features (RRC: standardised over all slide cells; HAC: [S, logit q] standardised over all slide cells)
        featP, featS = slide_stats(P), slide_stats(S)
        lq = np.log(np.clip(q, 1e-4, 1 - 1e-4) / (1 - np.clip(q, 1e-4, 1 - 1e-4))).astype(np.float32)
        featL = slide_stats(lq)
        tree = cKDTree(d["xy_um"]); pcs = prediction_pcs(P); tissue = occupied_area_mm2(d["xy_um"])
        rows, genes_tab = [], []
        configs = [(1.0, 4, "diverse", 0)] if args.smoke else [(ar, n, pl, s) for ar, n, pl in CONFIGS for s in a["repeats"]]
        for ar, n, pl, s in configs:
            centres, r, aidx, ag = place_cores(d["xy_um"], ar, n, pl, s, P=P, pcs=pcs, tree=tree)
            eidx = np.flatnonzero(evaluation_mask(d["xy_um"], centres, r, a["evaluation_buffer_um"]))
            Xa = np.hstack([(S[aidx] - featS[0]) / featS[1], (lq[aidx] - featL[0]) / featL[1]]).astype(np.float32)
            t1 = time.time()
            rrc = adapt.fit_rrc_authors(P[aidx], Y[aidx], *featP)
            rrs = adapt.fit_rrc_authors(S[aidx], Y[aidx], *featS)
            hac = hac_fit(Xa, q[aidx], R[aidx], Y[aidx], ag)
            fit_s = time.time() - t1
            # P0 diagnostic: optimism of random-cell hyperparameter selection
            x_anc = ((P[aidx] - featP[0]) / featP[1]).astype(np.float32)
            from one_shot_rrc import choose_hyperparams
            n_a = len(aidx); cap = n_a / (n_a + 1000.0)
            best, _ = choose_hyperparams(x_anc, P[aidx], Y[aidx], adapt.RRC_ALPHAS,
                                         sorted(set([b for b in adapt.RRC_BLENDS if b <= cap] + [cap])), 1, 7, 0.2)
            loco = rrc_loco_score(P[aidx], Y[aidx], ag, featP, rrc)
            Ye, Pe, Se, qe, Re = Y[eidx], P[eidx], S[eidx], q[eidx], R[eidx]
            Xe = np.hstack([(Se - featS[0]) / featS[1], (lq[eidx] - featL[0]) / featL[1]]).astype(np.float32)
            preds = {"GHIST+": Pe, "GHIST+-Soft": Se, "RRC-authors": adapt.predict_rrc_authors(rrc, Pe),
                     "RRC-soft": adapt.predict_rrc_authors(rrs, Se), "HAC": _hac_predict(hac, Xe, qe, Re)}
            del Xe
            res, pg = evaluate(Ye, preds, groups, d["xy_um"][eidx], s)
            key = f"a{ar}_n{n}_{pl}_s{s}"
            for meth, rw in res.items():
                rw.update(slide=slide, config=key, core_area_mm2=ar, n_cores=n, placement=pl, seed=s, method=meth,
                          n_anchor_cells=int(len(aidx)), n_eval_cells=int(len(eidx)), fit_seconds=round(fit_s, 1),
                          hac_alpha=hac["alpha"], rrc_alpha=rrc["alpha"], rrc_blend=rrc["blend"],
                          rrc_cv_random20=float(best["cv_mean_pcc"]), rrc_cv_loco=loco)
                rows.append(rw)
                genes_tab.append(pd.DataFrame({"config": key, "method": meth, "gene_idx": np.arange(len(pg[meth])), "pcc": pg[meth]}))
            res_rrc_eval = res["RRC-authors"]["mean_pcc_all"]
            print(f"[{slide}] {key}: anchors={len(aidx)} | medPCC GHIST+ {res['GHIST+']['median_pcc_all']:.4f} "
                  f"Soft {res['GHIST+-Soft']['median_pcc_all']:.4f} RRC {res['RRC-authors']['median_pcc_all']:.4f} "
                  f"RRC-soft {res['RRC-soft']['median_pcc_all']:.4f} HAC {res['HAC']['median_pcc_all']:.4f} "
                  f"(HAC-RRC CI [{res['HAC']['ci_vs_RRC-authors_all_lo']:+.4f},{res['HAC']['ci_vs_RRC-authors_all_hi']:+.4f}]) | "
                  f"RRC meanPCC cv-random {best['cv_mean_pcc']:.4f} / LOCO {loco:.4f} / eval {res_rrc_eval:.4f} ({time.time()-t0:.0f}s)", flush=True)
            pd.DataFrame(rows).to_csv(out / "runs.csv", index=False, float_format="%.6g")
            pd.concat(genes_tab).to_csv(out / "per_gene.csv.gz", index=False, float_format="%.5g")
            del preds, res, pg, Ye, Pe, Se, qe, Re, Xa; gc.collect()  # free this run's arrays before the next run (memory)
        print(f"[{slide}] done ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
