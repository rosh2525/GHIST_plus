"""GHIST+ hurdle outputs aligned to the v3 slide caches.

GHIST+ exports pred = top-K gate x ReLU(signed_mu), with K_g = round(sum_i q_ig) (model/hurdle_distribution.py,
deterministic_hurdle_prediction): every cell outside the top-K is exactly 0. The released .npz files also contain
the ungated parts: occurrence_probability q and signed_mu_scaled (same scale as pred_expr_scaled).
"""
from __future__ import annotations

import numpy as np

from slides import CACHE, SLIDES, load


def build(name: str):
    d = load(name)
    s = SLIDES[name]
    if s["pred"].suffix != ".npz":
        raise ValueError(f"{name}: export has no hurdle outputs")
    out = {}
    with np.load(s["pred"], allow_pickle=True) as z:
        pid = z["cell_ids"].astype(np.int64)
        pg = z["gene_names"].astype(str).tolist()
        order = np.argsort(pid)
        rows = order[np.searchsorted(pid[order], d["ids"])]
        cols = [pg.index(g) for g in d["genes"]]
        for key, dst in (("occurrence_probability", "q"), ("signed_mu_scaled", "mu")):
            a = z[key]
            out[dst] = a[np.ix_(rows, cols)].astype(np.float32)
            del a
        out["ctp"] = z["ct_probs"][rows].astype(np.float32)
        out["ct_classes"] = z["ct_classes"].astype("U64") if "ct_classes" in z.files else np.array([], dtype="U64")
    out["mu"] *= np.float32(s["pred_scale"])
    # sanity: the gated export must equal topK(q) x relu(mu) on the cached cohort subset's genes up to cohort effects;
    # here we only check that exported positives coincide with relu(mu) values
    pos = d["pred"] > 0
    err = float(np.abs(d["pred"][pos] - np.maximum(out["mu"][pos], 0)).max()) if pos.any() else 0.0
    np.savez(CACHE / f"{name}_hurdle.npz", ids=d["ids"], **out)
    print(f"{name}: hurdle cache written; max |pred - relu(mu)| on exported positives = {err:.2e}")


def load_hurdle(name: str, ids: np.ndarray) -> dict:
    p = CACHE / f"{name}_hurdle.npz"
    if not p.exists():
        build(name)
    with np.load(p) as z:
        if not np.array_equal(z["ids"], ids):
            raise ValueError("hurdle cache out of sync")
        return {k: z[k] for k in ("q", "mu", "ctp", "ct_classes") if k in z.files}


if __name__ == "__main__":
    import sys
    for n in sys.argv[1:]:
        build(n)
