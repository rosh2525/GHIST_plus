"""P3 (config_v7.json): conformal intervals around GHIST+, GHIST+-Soft, authors' RRC and HAC.
Usage: python p3_run.py breast5k [atera] [--smoke]   ->  outputs_v7/<slide>/runs.csv"""
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
from one_shot_rrc import fit_rrc as rrc_fit
from p1_run import _hac_predict, hac_fit, slide_stats
from p2_run import biased_cores, hac_parts, support_lambda, tumour_idx
from run_v3 import CFG
from slides import load

HERE = Path(__file__).resolve().parent
ALPHA = 0.1


def radius(res, n_min=1):
    n = len(res)
    if n < n_min:
        return None
    lvl = min(1.0, np.ceil((n + 1) * (1 - ALPHA)) / n)
    return np.quantile(res, lvl, axis=0, method="higher")


def rrc_oof(P, Y, groups, feat, r):
    x = ((P - feat[0]) / feat[1]).astype(np.float32)
    oof = np.zeros_like(Y, dtype=np.float32)
    for f in np.unique(groups):
        te, tr = groups == f, groups != f
        xm, ym, coef = rrc_fit(x[tr], (Y[tr] - P[tr]).astype(np.float32), r["alpha"])
        oof[te] = np.clip(P[te] + r["blend"] * ((x[te] - xm) @ coef + ym), 0, None)
    return oof


def cover(Y, Q, r):
    lo = np.maximum(Q - r, 0); hi = Q + r
    inside = (Y >= lo) & (Y <= hi)
    pos = Y > 0
    cov_all = inside.mean(0)
    with np.errstate(invalid="ignore"):
        cov_pos = np.where(pos.sum(0) > 0, (inside & pos).sum(0) / np.maximum(pos.sum(0), 1), np.nan)
    return cov_all, cov_pos, float((hi - lo).mean())


def calib(Y, Qp, chunk=50_000, nb=10):
    """Pooled ECE (equal-width bins) and per-gene Brier of probabilities Qp against 1[Y > 0]."""
    cnt = np.zeros(nb); sp = np.zeros(nb); so = np.zeros(nb)
    G = Y.shape[1]; bsum = np.zeros(G); n = 0
    for s0 in range(0, len(Y), chunk):
        o = (Y[s0:s0 + chunk] > 0).astype(np.float32); p = Qp[s0:s0 + chunk]
        b = np.clip((p * nb).astype(int), 0, nb - 1)
        cnt += np.bincount(b.ravel(), minlength=nb); sp += np.bincount(b.ravel(), weights=p.ravel(), minlength=nb)
        so += np.bincount(b.ravel(), weights=o.ravel(), minlength=nb)
        bsum += ((p - o) ** 2).sum(0); n += len(o)
    ece = float(np.sum(np.abs(sp - so)) / cnt.sum())
    return ece, float(np.median(bsum / n))


def mag_interval(Ya, Ma, Ye, Me, min_pos=20):
    """Conditional magnitude interval on expressing entries: radius from expressing core cells, coverage on expressing eval entries."""
    G = Ya.shape[1]; covs, widths = [], []
    for g in range(G):
        pa = Ya[:, g] > 0; pe = Ye[:, g] > 0
        if pa.sum() < min_pos or pe.sum() < min_pos:
            continue
        r = radius(np.abs(Ya[pa, g] - Ma[pa, g])[:, None])[0]
        lo = np.maximum(Me[pe, g] - r, 0); hi = Me[pe, g] + r
        covs.append(float(((Ye[pe, g] >= lo) & (Ye[pe, g] <= hi)).mean())); widths.append(float((hi - lo).mean()))
    return (float(np.median(covs)) if covs else np.nan), (float(np.mean(widths)) if widths else np.nan), len(covs)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("slides", nargs="+"); ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    a = CFG["anchors"]
    for slide in args.slides:
        t0 = time.time()
        d = load(slide); h = load_hurdle(slide, d["ids"])
        Y, P, q = d["y"], d["pred"], h["q"]; R = np.maximum(h.pop("mu"), 0); S = (q * R).astype(np.float32)
        ct = h["ctp"].argmax(1); TUMOUR_CLASSES = tumour_idx(h); tumour = h["ctp"][:, list(TUMOUR_CLASSES)].sum(1)
        featP, featS = slide_stats(P), slide_stats(S)
        lq = np.log(np.clip(q, 1e-4, 1 - 1e-4) / (1 - np.clip(q, 1e-4, 1 - 1e-4))).astype(np.float32); featL = slide_stats(lq)
        rng = np.random.default_rng(0)
        sub = ((S[rng.choice(len(S), 50_000, replace=False)] - featS[0]) / featS[1]).astype(np.float64)
        _, _, Vt = np.linalg.svd(sub - sub.mean(0), full_matrices=False)
        Z = (((S - featS[0]) / featS[1]) @ Vt[:10].T.astype(np.float32)).astype(np.float32)
        tree = cKDTree(d["xy_um"]); pcs = prediction_pcs(P)
        out = HERE / ("outputs_v7_smoke" if args.smoke else "outputs_v7") / slide; out.mkdir(parents=True, exist_ok=True)
        rows = []
        runs = [("unbiased", 0)] if args.smoke else [(st, s) for st in ("unbiased", "tumour_biased") for s in a["repeats"]]
        for setting, seed in runs:
            if setting == "unbiased":
                centres, r, aidx, ag = place_cores(d["xy_um"], 1.0, 4, "diverse", seed, P=P, pcs=pcs, tree=tree)
            else:
                centres, r, aidx, ag, _ = biased_cores(d["xy_um"], tree, tumour, 1.0, 4, seed, True)
            eidx = np.flatnonzero(evaluation_mask(d["xy_um"], centres, r, a["evaluation_buffer_um"]))
            Ya = Y[aidx]
            Xa = np.hstack([(S[aidx] - featS[0]) / featS[1], (lq[aidx] - featL[0]) / featL[1]]).astype(np.float32)
            hac = hac_fit(Xa, q[aidx], R[aidx], Ya, ag)
            rrc = adapt.fit_rrc_authors(P[aidx], Ya, *featP)
            lam, _, _ = support_lambda(Z, aidx, ag, eidx)
            Xe = np.hstack([(S[eidx] - featS[0]) / featS[1], (lq[eidx] - featL[0]) / featL[1]]).astype(np.float32)
            preds = {"GHIST+": P[eidx], "GHIST+-Soft": S[eidx], "RRC-authors": adapt.predict_rrc_authors(rrc, P[eidx]),
                     "HAC": _hac_predict(hac, Xe, q[eidx], R[eidx])}
            del Xe
            oofs = {"GHIST+": P[aidx], "GHIST+-Soft": S[aidx], "RRC-authors": rrc_oof(P[aidx], Ya, ag, featP, rrc), "HAC": hac["oof"]}
            radii = {m: radius(np.abs(Ya - o)) for m, o in oofs.items()}
            Ye = Y[eidx]; cte = ct[eidx]; is_t = np.isin(cte, TUMOUR_CLASSES); out_s = lam < 1
            # Mondrian radii for HAC (per gene x predicted cell type)
            res_h = np.abs(Ya - hac["oof"]); cta = ct[aidx]
            r_m = np.tile(radii["HAC"], (len(eidx), 1)).astype(np.float32)
            for k in np.unique(cta):
                rk = radius(res_h[cta == k], n_min=50)
                if rk is not None:
                    r_m[cte == k] = rk
            entries = list(radii.items()) + [("HAC-mondrian", r_m)]
            base = dict(slide=slide, setting=setting, seed=seed, n_anchor_cells=int(len(aidx)), n_eval_cells=int(len(eidx)))
            for m, rr_ in entries:
                Q = preds["HAC"] if m == "HAC-mondrian" else preds[m]
                rr_b = rr_ if np.ndim(rr_) == 2 else rr_[None, :]
                ca, cp, w = cover(Ye, Q, rr_b)
                row = dict(base, method=m, coverage_all=float(np.nanmedian(ca)), coverage_expressing=float(np.nanmedian(cp)), mean_width=w)
                for name, sel in (("tumour", is_t), ("other", ~is_t), ("out_of_support", out_s), ("in_support", ~out_s)):
                    if sel.sum() >= 200:
                        ca_s, cp_s, w_s = cover(Ye[sel], Q[sel], rr_b if rr_b.shape[0] == 1 else rr_b[sel])
                        row[f"coverage_all_{name}"] = float(np.nanmedian(ca_s)); row[f"coverage_expressing_{name}"] = float(np.nanmedian(cp_s))
                rows.append(row)
            # amendment 1: hurdle-aware uncertainty (occurrence calibration + conditional magnitude interval)
            Xe2 = np.hstack([(S[eidx] - featS[0]) / featS[1], (lq[eidx] - featL[0]) / featL[1]]).astype(np.float32)
            term, mh = hac_parts(hac, Xe2, q[eidx], R[eidx]); del Xe2
            qe = q[eidx]
            qhat = np.clip(qe + hac["c"] + term, 0.001, 0.999); qcal = np.clip(qe + hac["c"], 0.001, 0.999); del term
            ma = np.maximum(hac["a"] + hac["b"] * R[aidx], 0.01)
            hrow = dict(base, method="HURDLE-UNCERTAINTY")
            for nm, Qp in (("raw_q", qe), ("q_plus_c", qcal), ("hac_qhat", qhat)):
                e_, b_ = calib(Ye, Qp); hrow[f"ece_{nm}"] = e_; hrow[f"brier_{nm}"] = b_
            hrow["mag_cov_hac"], hrow["mag_width_hac"], hrow["mag_genes"] = mag_interval(Ya, ma, Ye, mh)
            hrow["mag_cov_raw"], hrow["mag_width_raw"], _ = mag_interval(Ya, R[aidx], Ye, R[eidx])
            rows.append(hrow)
            show = {x["method"]: x for x in rows[-6:-1]}
            print(f"[{slide}] {setting} s{seed}: " + " | ".join(f"{m}: cov {v['coverage_all']:.3f} (expr {v['coverage_expressing']:.3f}) w {v['mean_width']:.3f}" for m, v in show.items())
                  + f" || ECE raw {hrow['ece_raw_q']:.3f} -> HAC {hrow['ece_hac_qhat']:.3f}; Brier {hrow['brier_raw_q']:.4f} -> {hrow['brier_hac_qhat']:.4f}; "
                  f"magnitude interval cov {hrow['mag_cov_hac']:.3f} w {hrow['mag_width_hac']:.3f} (raw cov {hrow['mag_cov_raw']:.3f} w {hrow['mag_width_raw']:.3f}, {hrow['mag_genes']} genes)"
                  + f" ({time.time()-t0:.0f}s)", flush=True)
            pd.DataFrame(rows).to_csv(out / "runs.csv", index=False, float_format="%.5f")
            del preds, oofs, Ye, Ya, Xa, r_m, qhat, qcal, mh; gc.collect()  # free this run's arrays before the next run (memory)
        print(f"[{slide}] done ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
