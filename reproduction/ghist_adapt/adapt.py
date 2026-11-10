"""GHIST+-Adapt: fitted only on core (anchor) cells. No evaluation-cell RNA is used anywhere here."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from metrics import pcc_columns

# the authors' released one-shot ROI residual correction (tools/one_shot_rrc.py), imported unchanged
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from one_shot_rrc import choose_hyperparams as _rrc_choose, fit_rrc as _rrc_fit  # noqa: E402

RRC_ALPHAS = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000]          # script defaults
RRC_BLENDS = [0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
RRC_BLEND_CAP_N0, RRC_N_FOLDS, RRC_SEED, RRC_VAL_FRACTION = 1000.0, 1, 7, 0.2


def fit_rrc_authors(P_anchor, Y_anchor, feat_mean, feat_std):
    """Authors' RRC with its default config: features standardised over all slide cells, alpha/blend chosen on a
    random 20% of ROI cells, blend capped at n/(n+1000)."""
    x = ((P_anchor - feat_mean) / feat_std).astype(np.float32)
    p = P_anchor.astype(np.float32); y = Y_anchor.astype(np.float32)
    n = len(x)
    cap = float(n / (n + RRC_BLEND_CAP_N0))
    blends = sorted(set([b for b in RRC_BLENDS if b <= cap] + [cap]))
    best, _ = _rrc_choose(x, p, y, RRC_ALPHAS, blends, RRC_N_FOLDS, RRC_SEED, RRC_VAL_FRACTION)
    x_mean, y_mean, coef = _rrc_fit(x, y - p, best["alpha"])
    return dict(x_mean=x_mean, y_mean=y_mean, coef=coef, blend=float(best["blend"]), alpha=float(best["alpha"]),
                feat_mean=feat_mean.astype(np.float32), feat_std=feat_std.astype(np.float32))


def predict_rrc_authors(r, P):
    P = P.astype(np.float32)
    z = (P - r["feat_mean"]) / r["feat_std"]
    out = P + r["blend"] * ((z - r["x_mean"]) @ r["coef"] + r["y_mean"])
    return np.clip(out, 0, None).astype(np.float32)


def _svd_ridge(Z, R):
    zm, rm = Z.mean(0), R.mean(0)
    U, S, Vt = np.linalg.svd(Z - zm, full_matrices=False)
    return zm, rm, S, Vt, U.T @ (R - rm)


def _weights(svd, alpha):
    zm, rm, S, Vt, UtR = svd
    W = Vt.T @ ((S / (S ** 2 + alpha))[:, None] * UtR)
    return W, rm - zm @ W


def fit(P, Y, groups, alphas, conformal_alpha=0.1):
    P = P.astype(np.float64); Y = Y.astype(np.float64)
    G = P.shape[1]
    mu, sd = P.mean(0), P.std(0)
    sd[sd == 0] = 1.0
    Z = (P - mu) / sd
    R = Y - P
    folds = np.unique(groups)

    # pass 1: out-of-fold predictions for every alpha (pooled PCC selects alpha)
    oof = np.zeros((len(alphas), *Y.shape), dtype=np.float32)
    for f in folds:
        te, tr = groups == f, groups != f
        svd = _svd_ridge(Z[tr], R[tr])
        for a, alpha in enumerate(alphas):
            W, c = _weights(svd, alpha)
            oof[a, te] = np.clip(P[te] + c + Z[te] @ W, 0, None)
    scores = [np.nanmedian(pcc_columns(Y, oof[a])) for a in range(len(alphas))]
    best = int(np.nanargmax(scores))
    alpha = alphas[best]
    q_oof = oof[best].astype(np.float64)
    del oof
    W, c = _weights(_svd_ridge(Z, R), alpha)

    # per-gene gate: ridge if its out-of-fold PCC beats GHIST+ on the cores
    pcc_ridge = pcc_columns(Y, q_oof)
    pcc_ghist = pcc_columns(Y, P)
    use_ridge = np.nan_to_num(pcc_ridge, nan=-np.inf) > np.nan_to_num(pcc_ghist, nan=-np.inf)

    # tie-free rescale of GHIST+ (least squares through the origin)
    pp = (P * P).sum(0)
    s = np.where(pp > 0, (P * Y).sum(0) / np.where(pp > 0, pp, 1), 0.0)

    # calibration baseline
    pm, ym = P.mean(0), Y.mean(0)
    var = ((P - pm) ** 2).mean(0)
    b = np.where(var > 0, ((P - pm) * (Y - ym)).mean(0) / np.where(var > 0, var, 1), 0.0)
    a_cal = ym - b * pm

    # distribution matching (spread) parameters
    sdy = Y.std(0)
    sdq = q_oof.std(0)
    sdp = P.std(0)
    k_ridge = np.clip(np.where(sdq > 0, sdy / np.where(sdq > 0, sdq, 1), 1.0), 0, 20)
    m_ridge = q_oof.mean(0)
    k_scale = np.where(sdp > 0, sdy / np.where(sdp > 0, sdp, 1), s)

    # split-conformal radii from out-of-fold residuals (Adapt) and raw residuals (GHIST+ baseline)
    adapt_oof = np.where(use_ridge[None, :], q_oof, P * s[None, :])
    n = len(Y)
    level = min(1.0, np.ceil((n + 1) * (1 - conformal_alpha)) / n)
    r_adapt = np.quantile(np.abs(Y - adapt_oof), level, axis=0, method="higher")
    r_ghist = np.quantile(np.abs(Y - P), level, axis=0, method="higher")

    return dict(mu=mu, sd=sd, W=W.astype(np.float32), c=c, alpha=alpha, cv_scores=scores, use_ridge=use_ridge,
                s=s, a_cal=a_cal, b_cal=b, k_ridge=k_ridge, m_ridge=m_ridge, k_scale=k_scale,
                r_adapt=r_adapt, r_ghist=r_ghist, n_genes_ridge=int(use_ridge.sum()), G=G)


def predict(m, P, method):
    P = P.astype(np.float32)
    if method == "GHIST+":
        return P
    if method == "calibration":
        return np.clip(m["a_cal"] + m["b_cal"] * P, 0, None).astype(np.float32)
    ridge = np.clip(P + m["c"] + ((P - m["mu"]) / m["sd"]).astype(np.float32) @ m["W"], 0, None).astype(np.float32)
    if method == "ridge":
        return ridge
    ur = m["use_ridge"][None, :]
    if method == "Adapt":
        return np.where(ur, ridge, P * m["s"][None, :]).astype(np.float32)
    if method == "Adapt-DM":
        dm_ridge = np.clip(m["m_ridge"] + (ridge - m["m_ridge"]) * m["k_ridge"], 0, None)
        return np.where(ur, dm_ridge, P * m["k_scale"][None, :]).astype(np.float32)
    raise ValueError(method)
