"""TMA-like disk cores: random or RNA-free 'diverse' placement, CV groups and evaluation mask."""
from __future__ import annotations

import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.spatial import cKDTree


def occupied_area_mm2(xy_um, grid_um=100.0):
    cells = np.unique(np.floor(xy_um / grid_um).astype(np.int64), axis=0)
    return len(cells) * (grid_um / 1000.0) ** 2


def prediction_pcs(P, n_pc=20, n_fit=50_000, seed=0):
    rng = np.random.default_rng(seed)
    sub = P[rng.choice(len(P), min(n_fit, len(P)), replace=False)].astype(np.float64)
    mu, sd = sub.mean(0), sub.std(0)
    sd[sd == 0] = 1.0
    _, _, Vt = np.linalg.svd((sub - mu) / sd, full_matrices=False)
    V = Vt[:n_pc].T
    return lambda X: ((X.astype(np.float64) - mu) / sd) @ V


def place_cores(xy, area_mm2, n_cores, placement, seed, P=None, pcs=None, tree=None):
    rng = np.random.default_rng(seed)
    r = float(np.sqrt(area_mm2 * 1e6 / np.pi))
    tree = tree or cKDTree(xy)
    ref = tree.query_ball_point(xy[rng.choice(len(xy), 200)], r, return_length=True)
    min_count = 0.5 * float(np.median(ref))

    def ok(c, chosen):
        return all(np.hypot(*(c - s)) >= 2 * r for s in chosen)

    chosen = []
    if placement == "random":
        tries = 0
        while len(chosen) < n_cores and tries < 50_000:
            tries += 1
            c = xy[rng.integers(len(xy))]
            if ok(c, chosen) and tree.query_ball_point(c, r, return_length=True) >= min_count:
                chosen.append(c)
    elif placement == "diverse":
        cand, members = [], []
        tries = 0
        while len(cand) < 400 and tries < 50_000:
            tries += 1
            c = xy[rng.integers(len(xy))]
            m = tree.query_ball_point(c, r)
            if len(m) >= min_count:
                cand.append(c); members.append(np.asarray(m))
        cand = np.asarray(cand)
        desc = np.stack([pcs(P[m]).mean(0) for m in members])
        if n_cores == 1:
            order = np.argsort(np.linalg.norm(desc - desc.mean(0), axis=1))
            chosen = [cand[order[0]]]
        else:
            np.random.seed(seed)  # kmeans2 uses the global RNG in older SciPy
            cent, lab = kmeans2(desc, n_cores, minit="++", seed=seed)
            sizes = np.bincount(lab, minlength=n_cores)
            for k in np.argsort(-sizes):
                for i in np.argsort(np.linalg.norm(desc - cent[k], axis=1)):
                    if ok(cand[i], chosen):
                        chosen.append(cand[i]); break
    else:
        raise ValueError(placement)
    if len(chosen) < n_cores:
        raise RuntimeError(f"could only place {len(chosen)}/{n_cores} cores of {area_mm2} mm2")
    centres = np.asarray(chosen)
    idx, groups = [], []
    for k, c in enumerate(centres):
        m = np.asarray(tree.query_ball_point(c, r), dtype=np.int64)
        idx.append(m)
        if n_cores == 1:  # 4 angular sectors as CV groups
            ang = np.arctan2(xy[m, 1] - c[1], xy[m, 0] - c[0])
            groups.append(np.floor((ang + np.pi) / (np.pi / 2)).clip(0, 3).astype(int))
        else:
            groups.append(np.full(len(m), k))
    idx = np.concatenate(idx)
    groups = np.concatenate(groups)
    return centres, r, idx, groups


def evaluation_mask(xy, centres, r, buffer_um):
    d, _ = cKDTree(centres).query(xy)
    return d > r + buffer_um
