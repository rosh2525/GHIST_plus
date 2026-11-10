"""Metrics for GHIST+-Adapt v3. giotto_scores, the SSIM grid/map/local functions and CMD (Eq. 17)
are copied from the authors' Figure3.ipynb so values are comparable with the paper's definitions."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter, uniform_filter
from scipy.spatial import cKDTree
from scipy.stats import rankdata


# ------------------------------------------------------------------ authors' code (Figure3.ipynb)
def giotto_scores(expression, coordinates, k=8):
    neighbours = cKDTree(coordinates).query(coordinates, k=k + 1)[1][:, 1:]
    scores = np.zeros(expression.shape[1], dtype=np.float32)
    for index in range(expression.shape[1]):
        ranks = rankdata(expression[:, index], method="average").astype(np.float32)
        near = ranks[neighbours].mean(axis=1)
        ranks -= ranks.mean()
        near -= near.mean()
        denominator = np.sqrt(np.dot(ranks, ranks) * np.dot(near, near))
        if denominator > 0:
            scores[index] = np.dot(ranks, near) / denominator
    return scores


@dataclass(frozen=True)
class SsimGrid:
    height: int
    width: int
    flat_index: np.ndarray
    count: np.ndarray
    crop: tuple


def build_ssim_grid(x, y, width=400, minimum_height=20):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    xmin, xmax = float(x.min()), float(x.max())
    ymin, ymax = float(y.min()), float(y.max())
    height = max(minimum_height, int(round(width * (ymax - ymin) / (xmax - xmin))))
    xi = np.clip((((x - xmin) / (xmax - xmin)) * (width - 1)).astype(np.int32), 0, width - 1)
    yi = np.clip((((y - ymin) / (ymax - ymin)) * (height - 1)).astype(np.int32), 0, height - 1)
    flat_index = yi * width + xi
    count = np.bincount(flat_index, minlength=height * width).reshape(height, width).astype(np.float64)
    occupied_y, occupied_x = np.where(count > 0)
    crop = (slice(int(occupied_y.min()), int(occupied_y.max()) + 1), slice(int(occupied_x.min()), int(occupied_x.max()) + 1))
    return SsimGrid(height, width, flat_index, count, crop)


def ssim_spatial_map(values, grid, sigma=1.0):
    sums = np.bincount(grid.flat_index, weights=np.asarray(values, dtype=np.float64), minlength=grid.height * grid.width).reshape(grid.height, grid.width)
    mean_map = np.divide(sums, grid.count, out=np.zeros_like(sums), where=grid.count > 0)
    return gaussian_filter(mean_map, sigma=sigma, mode="nearest", truncate=4.0)[grid.crop]


def ssim_local(target_map, prediction_map, data_range, window=7, k1=0.01, k2=0.03):
    if not np.isfinite(data_range) or data_range <= 0:
        return np.nan
    arguments = {"size": window, "mode": "reflect", "origin": 0}
    mean_target = uniform_filter(target_map, **arguments)
    mean_prediction = uniform_filter(prediction_map, **arguments)
    mean_target_sq = uniform_filter(target_map * target_map, **arguments)
    mean_prediction_sq = uniform_filter(prediction_map * prediction_map, **arguments)
    mean_cross = uniform_filter(target_map * prediction_map, **arguments)
    sample_factor = window ** 2 / (window ** 2 - 1.0)
    variance_target = sample_factor * (mean_target_sq - mean_target ** 2)
    variance_prediction = sample_factor * (mean_prediction_sq - mean_prediction ** 2)
    covariance = sample_factor * (mean_cross - mean_target * mean_prediction)
    c1 = (k1 * data_range) ** 2
    c2 = (k2 * data_range) ** 2
    numerator = (2 * mean_target * mean_prediction + c1) * (2 * covariance + c2)
    denominator = (mean_target ** 2 + mean_prediction ** 2 + c1) * (variance_target + variance_prediction + c2)
    score_map = np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator != 0)
    border = (window - 1) // 2
    interior = score_map[border:-border, border:-border]
    return float(np.mean(interior)) if np.isfinite(interior).all() else np.nan


def cmd_from_corr(target_correlation, prediction_correlation, epsilon=1e-12):
    numerator = float(np.trace(prediction_correlation @ target_correlation))
    denominator = float(np.linalg.norm(prediction_correlation, "fro") * np.linalg.norm(target_correlation, "fro") + epsilon)
    return 1.0 - numerator / denominator


# ------------------------------------------------------------------ chunked per-gene statistics
def gene_sums(Y, Q, chunk=100_000):
    G = Y.shape[1]
    s = {k: np.zeros(G) for k in ("sy", "sq", "syy", "sqq", "syq", "se", "npos", "sepos")}
    for a in range(0, len(Y), chunk):
        y = Y[a:a + chunk].astype(np.float64)
        q = Q[a:a + chunk].astype(np.float64)
        e2 = (y - q) ** 2
        pos = y > 0
        s["sy"] += y.sum(0); s["sq"] += q.sum(0)
        s["syy"] += (y * y).sum(0); s["sqq"] += (q * q).sum(0); s["syq"] += (y * q).sum(0)
        s["se"] += e2.sum(0); s["npos"] += pos.sum(0); s["sepos"] += np.where(pos, e2, 0).sum(0)
    s["n"] = len(Y)
    return s


def from_sums(s):
    n = s["n"]
    vy = s["syy"] / n - (s["sy"] / n) ** 2
    vq = s["sqq"] / n - (s["sq"] / n) ** 2
    cov = s["syq"] / n - s["sy"] * s["sq"] / n ** 2
    ok = (vy > 1e-12) & (vq > 1e-12)
    with np.errstate(invalid="ignore", divide="ignore"):
        pcc = np.where(ok, cov / np.sqrt(np.clip(vy, 1e-30, None) * np.clip(vq, 1e-30, None)), np.nan)
        sd_ratio = np.where(vy > 1e-12, np.sqrt(np.clip(vq, 0, None) / np.clip(vy, 1e-30, None)), np.nan)
        rmse_pos = np.where(s["npos"] > 0, np.sqrt(s["sepos"] / np.maximum(s["npos"], 1)), np.nan)
    return {"pcc": pcc, "rmse": np.sqrt(s["se"] / n), "rmse_pos": rmse_pos, "sd_ratio": sd_ratio}


def pcc_columns(Y, Q):
    return from_sums(gene_sums(Y, Q))["pcc"]


def spearman_columns(Y_ranks, Q):
    out = np.full(Q.shape[1], np.nan)
    for g in range(Q.shape[1]):
        rq = rankdata(Q[:, g], method="average")
        ry = Y_ranks[:, g]
        ry_c, rq_c = ry - ry.mean(), rq - rq.mean()
        den = np.sqrt((ry_c @ ry_c) * (rq_c @ rq_c))
        if den > 0:
            out[g] = (ry_c @ rq_c) / den
    return out


def corr_matrix(X, cols, chunk=100_000):
    k = len(cols)
    s = np.zeros(k); ss = np.zeros((k, k))
    for a in range(0, len(X), chunk):
        x = X[a:a + chunk][:, cols].astype(np.float64)
        s += x.sum(0); ss += x.T @ x
    n = len(X)
    cov = ss / n - np.outer(s, s) / n ** 2
    d = np.sqrt(np.clip(np.diag(cov), 1e-30, None))
    return cov / np.outer(d, d)


# ------------------------------------------------------------------ spatial block bootstrap
def block_stats(Y, Q, order, starts, gchunk=32):
    """Per-block sums (sorted by block) for PCC; (n_blocks, G) arrays, computed in gene chunks to bound memory."""
    G = Y.shape[1]
    out = {k: np.zeros((len(starts), G)) for k in ("sq", "sqq", "syq")}
    for a in range(0, G, gchunk):
        y = Y[order, a:a + gchunk].astype(np.float64)
        q = Q[order, a:a + gchunk].astype(np.float64)
        out["sq"][:, a:a + gchunk] = np.add.reduceat(q, starts)
        out["sqq"][:, a:a + gchunk] = np.add.reduceat(q * q, starts)
        out["syq"][:, a:a + gchunk] = np.add.reduceat(y * q, starts)
    return out


def y_block_stats(Y, order, starts, gchunk=32):
    G = Y.shape[1]
    out = {k: np.zeros((len(starts), G)) for k in ("sy", "syy")}
    for a in range(0, G, gchunk):
        y = Y[order, a:a + gchunk].astype(np.float64)
        out["sy"][:, a:a + gchunk] = np.add.reduceat(y, starts)
        out["syy"][:, a:a + gchunk] = np.add.reduceat(y * y, starts)
    out["n"] = np.add.reduceat(np.ones(len(order)), starts)
    return out


def bootstrap_median_diff(ystats, base, new, gene_sets, n_boot=500, seed=0):
    rng = np.random.default_rng(seed)
    nb = len(ystats["n"])
    W = rng.multinomial(nb, np.full(nb, 1.0 / nb), size=n_boot).astype(np.float64)

    def pcc(stats):
        n = (W @ ystats["n"])[:, None]
        sy, syy = W @ ystats["sy"], W @ ystats["syy"]
        sq, sqq, syq = W @ stats["sq"], W @ stats["sqq"], W @ stats["syq"]
        vy = syy / n - (sy / n) ** 2
        vq = sqq / n - (sq / n) ** 2
        cov = syq / n - sy * sq / n ** 2
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where((vy > 1e-12) & (vq > 1e-12), cov / np.sqrt(np.clip(vy * vq, 1e-30, None)), np.nan)

    pb, pn = pcc(base), pcc(new)
    out = {}
    for name, idx in gene_sets.items():
        diff = np.nanmedian(pn[:, idx], axis=1) - np.nanmedian(pb[:, idx], axis=1)
        out[name] = (float(np.percentile(diff, 2.5)), float(np.percentile(diff, 97.5)))
    return out
