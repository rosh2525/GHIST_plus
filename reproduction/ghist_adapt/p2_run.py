"""P2 (config_v5.json): support-aware trust gate on HAC, with a biased-core stress test.

Usage:  python p2_run.py breast5k [atera] [--smoke]
Writes  outputs_v5/<slide>/runs.csv
"""
from __future__ import annotations

import argparse
import gc
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import adapt
from anchors import evaluation_mask, place_cores, prediction_pcs
from hurdle import load_hurdle
from metrics import block_stats, bootstrap_median_diff, from_sums, gene_sums, y_block_stats
from p1_run import hac_fit, slide_stats


def hac_parts(m, X, q, R, chunk=100_000):
    """HAC pieces: cell-specific multigene occurrence term (X - xm) W, and magnitude m_hat.
    q_hat = clip(q + c + term) reproduces p1_run._hac_predict exactly."""
    term = np.empty(q.shape, dtype=np.float32); mh = np.empty(q.shape, dtype=np.float32)
    for s in range(0, len(q), chunk):
        term[s:s + chunk] = (X[s:s + chunk] - m["xm"]) @ m["W"]
        mh[s:s + chunk] = np.maximum(m["a"] + m["b"] * R[s:s + chunk], 0.01)
    return term, mh
from run_v3 import CFG
from slides import load

HERE = Path(__file__).resolve().parent
SETTINGS = ["unbiased", "tumour_biased", "stroma_biased"]
TUMOUR_CLASSES = (4, 5)  # ct_classes order: B, Myeloid, Stromal, Macrophage, Malignant, Epithelial, Plasma, T, Unassigned


def tumour_idx(h):
    """Indices of the Malignant/Epithelial classes, by name when the cache stores class names (other checkpoints)."""
    names = [str(c) for c in h.get("ct_classes", [])]
    idx = [i for i, c in enumerate(names) if c in ("Malignant", "Epithelial")]
    return tuple(idx) if idx else TUMOUR_CLASSES


def biased_cores(xy, tree, tumour_score, area_mm2, n_cores, seed, high):
    rng = np.random.default_rng(seed)
    r = float(np.sqrt(area_mm2 * 1e6 / np.pi))
    ref = tree.query_ball_point(xy[rng.choice(len(xy), 200)], r, return_length=True)
    min_count = 0.5 * float(np.median(ref))
    cand, mem = [], []
    while len(cand) < 400:
        c = xy[rng.integers(len(xy))]
        m = tree.query_ball_point(c, r)
        if len(m) >= min_count:
            cand.append(c); mem.append(np.asarray(m))
    score = np.array([tumour_score[m].mean() for m in mem])
    order = np.argsort(-score if high else score)[:100]           # top / bottom 25%
    chosen, idx, groups = [], [], []
    for i in rng.permutation(order):
        if all(np.hypot(*(cand[i] - c)) >= 2 * r for c in chosen):
            chosen.append(cand[i]); idx.append(mem[i]); groups.append(np.full(len(mem[i]), len(chosen) - 1))
        if len(chosen) == n_cores:
            break
    if len(chosen) < n_cores:
        raise RuntimeError("could not place biased cores")
    return np.asarray(chosen), r, np.concatenate(idx), np.concatenate(groups), float(score[order].mean())


def support_lambda(Z, aidx, groups, eidx, k=10, q=0.95, power=2.0):
    Za = Z[aidx]
    loco = np.empty(len(aidx))
    for f in np.unique(groups):
        te = groups == f
        d, _ = cKDTree(Za[~te]).query(Za[te], k=k)
        loco[te] = d.mean(1)
    d_ref = float(np.quantile(loco, q))
    d, _ = cKDTree(Za).query(Z[eidx], k=k)
    d = d.mean(1)
    return np.where(d <= d_ref, 1.0, (d_ref / np.maximum(d, 1e-12)) ** power).astype(np.float32), d, d_ref


def subset_metrics(Y, preds, sel, xy, seed, ci_pairs):
    out = {}
    if sel.sum() < 500:
        return out
    Ys = Y[sel]
    blk = np.floor(xy[sel] / 400.0).astype(np.int64)
    key = blk[:, 0] * 1_000_000 + blk[:, 1]
    order = np.argsort(key, kind="stable")
    starts = np.flatnonzero(np.r_[True, np.diff(key[order]) != 0])
    ys = y_block_stats(Ys, order, starts)
    stats = {}
    for m, Q in preds.items():
        gm = from_sums(gene_sums(Ys, Q[sel]))
        out[f"{m}|median_pcc"] = float(np.nanmedian(gm["pcc"])); out[f"{m}|mean_rmse"] = float(gm["rmse"].mean())
        if any(m in p for p in ci_pairs):
            stats[m] = block_stats(Ys, Q[sel], order, starts)
    for new, base in ci_pairs:
        lo, hi = bootstrap_median_diff(ys, stats[base], stats[new], {"all": np.arange(Y.shape[1])}, n_boot=500, seed=seed)["all"]
        out[f"ci_{new}_vs_{base}_lo"], out[f"ci_{new}_vs_{base}_hi"] = lo, hi
    out["n_cells"] = int(sel.sum())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slides", nargs="+")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    a = CFG["anchors"]
    for slide in args.slides:
        t0 = time.time()
        d = load(slide); h = load_hurdle(slide, d["ids"])
        Y, P, q = d["y"], d["pred"], h["q"]
        R = np.maximum(h.pop("mu"), 0); S = (q * R).astype(np.float32)
        tidx = tumour_idx(h)
        tumour = h["ctp"][:, list(tidx)].sum(1)
        is_tumour = np.isin(h["ctp"].argmax(1), tidx)
        featP, featS = slide_stats(P), slide_stats(S)
        lq = np.log(np.clip(q, 1e-4, 1 - 1e-4) / (1 - np.clip(q, 1e-4, 1 - 1e-4))).astype(np.float32)
        featL = slide_stats(lq)
        # RNA-free state space: 10 PCs of standardised S
        rng = np.random.default_rng(0)
        sub = ((S[rng.choice(len(S), 50_000, replace=False)] - featS[0]) / featS[1]).astype(np.float64)
        _, _, Vt = np.linalg.svd(sub - sub.mean(0), full_matrices=False)
        Z = (((S - featS[0]) / featS[1]) @ Vt[:10].T.astype(np.float32)).astype(np.float32)
        tree = cKDTree(d["xy_um"]); pcs = prediction_pcs(P)
        out = HERE / ("outputs_v5_smoke" if args.smoke else "outputs_v5") / slide
        out.mkdir(parents=True, exist_ok=True)
        rows = []
        runs = [("tumour_biased", 0)] if args.smoke else [(st, s) for st in SETTINGS for s in a["repeats"]]
        for setting, seed in runs:
            if setting == "unbiased":
                centres, r, aidx, ag = place_cores(d["xy_um"], 1.0, 4, "diverse", seed, P=P, pcs=pcs, tree=tree)
                core_tumour = float(tumour[aidx].mean())
            else:
                centres, r, aidx, ag, core_tumour = biased_cores(d["xy_um"], tree, tumour, 1.0, 4, seed, setting == "tumour_biased")
            eidx = np.flatnonzero(evaluation_mask(d["xy_um"], centres, r, a["evaluation_buffer_um"]))
            Xa = np.hstack([(S[aidx] - featS[0]) / featS[1], (lq[aidx] - featL[0]) / featL[1]]).astype(np.float32)
            rrc = adapt.fit_rrc_authors(P[aidx], Y[aidx], *featP)
            hac = hac_fit(Xa, q[aidx], R[aidx], Y[aidx], ag)
            lam, dist, d_ref = support_lambda(Z, aidx, ag, eidx)
            Ye, Se = Y[eidx], S[eidx]
            Xe = np.hstack([(Se - featS[0]) / featS[1], (lq[eidx] - featL[0]) / featL[1]]).astype(np.float32)
            term, mh = hac_parts(hac, Xe, q[eidx], R[eidx]); del Xe
            qc = q[eidx] + hac["c"]                                       # globally recalibrated occurrence
            preds = {"GHIST+": P[eidx], "GHIST+-Soft": Se, "RRC-authors": adapt.predict_rrc_authors(rrc, P[eidx]),
                     "GHIST+-q-mag": (np.clip(qc, 0.001, 0.999) * mh).astype(np.float32),
                     "HAC": (np.clip(qc + term, 0.001, 0.999) * mh).astype(np.float32),
                     "HAC-trust": (np.clip(qc + lam[:, None] * term, 0.001, 0.999) * mh).astype(np.float32)}  # amendment 2
            del term, qc
            xy = d["xy_um"][eidx]
            pairs = [("HAC-trust", "HAC"), ("HAC", "GHIST+-Soft")]
            base = dict(slide=slide, setting=setting, seed=seed, n_anchor_cells=int(len(aidx)), n_eval_cells=int(len(eidx)),
                        core_tumour_prob=core_tumour, slide_tumour_prob=float(tumour.mean()), d_ref=d_ref,
                        frac_out_of_support=float((lam < 1).mean()), mean_lambda=float(lam.mean()), hac_alpha=hac["alpha"])
            for name, sel in (("all", np.ones(len(eidx), bool)), ("out_of_support", lam < 1), ("in_support", lam >= 1),
                              ("tumour_cells", is_tumour[eidx]), ("other_cells", ~is_tumour[eidx])):
                res = subset_metrics(Ye, preds, sel, xy, seed, pairs)
                rows.append(dict(base, subset=name, **res))
            rr = {x["subset"]: x for x in rows[-5:]}
            o, al = rr["out_of_support"], rr["all"]
            print(f"[{slide}] {setting} s{seed}: core tumour-prob {core_tumour:.2f} (slide {tumour.mean():.2f}); out-of-support {100*base['frac_out_of_support']:.0f}% of eval cells | "
                  f"OUT medPCC Soft {o.get('GHIST+-Soft|median_pcc', np.nan):.4f} RRC {o.get('RRC-authors|median_pcc', np.nan):.4f} HAC {o.get('HAC|median_pcc', np.nan):.4f} "
                  f"trust {o.get('HAC-trust|median_pcc', np.nan):.4f} CI[{o.get('ci_HAC-trust_vs_HAC_lo', np.nan):+.4f},{o.get('ci_HAC-trust_vs_HAC_hi', np.nan):+.4f}] | "
                  f"ALL HAC {al['HAC|median_pcc']:.4f} trust {al['HAC-trust|median_pcc']:.4f} ({time.time()-t0:.0f}s)", flush=True)
            pd.DataFrame(rows).to_csv(out / "runs.csv", index=False, float_format="%.6g")
            del preds, Ye, Se, Xa; gc.collect()  # free this run's arrays before the next run (memory)
        print(f"[{slide}] done ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
