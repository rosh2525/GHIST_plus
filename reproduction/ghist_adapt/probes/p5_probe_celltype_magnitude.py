"""DEV-ONLY probe on Breast5K: cell-type-aware magnitude and occurrence for HAC."""
import sys; sys.path.insert(0, "/Users/roshank/Documents/Codex/2026-10-02/we-x20/GHIST_plus/reproduction/ghist_adapt")
import numpy as np
from scipy.spatial import cKDTree
from slides import load; from hurdle import load_hurdle
from metrics import gene_sums, from_sums
from anchors import place_cores, evaluation_mask, prediction_pcs
from p1_run import hac_fit, slide_stats
d = load("breast5k"); h = load_hurdle("breast5k", d["ids"])
Y, P, q = d["y"], d["pred"], h["q"]; R = np.maximum(h["mu"], 0); S = (q * R).astype(np.float32); C = h["ctp"]
lq = np.log(np.clip(q, 1e-4, 1 - 1e-4) / (1 - np.clip(q, 1e-4, 1 - 1e-4))).astype(np.float32)
st = [slide_stats(S), slide_stats(lq), slide_stats(C)]
def X(idx, ct):
    b = [(S[idx] - st[0][0]) / st[0][1], (lq[idx] - st[1][0]) / st[1][1]] + ([(C[idx] - st[2][0]) / st[2][1]] if ct else [])
    return np.hstack(b).astype(np.float32)
def mag_ct(Ya, Ra, Ca, lam=1.0):
    """per gene: y ~ a + b*relu(mu) + sum_c w_c * ctp_c on expressing core cells (small ridge on w)"""
    G = Ya.shape[1]; coef = np.zeros((G, 2 + Ca.shape[1]), np.float32)
    allpos = Ya[Ya > 0]; glob = float(allpos.mean())
    for g in range(G):
        p = Ya[:, g] > 0
        if p.sum() < 20: coef[g, 0] = Ya[p, g].mean() if p.any() else glob; continue
        A = np.hstack([np.ones((p.sum(), 1)), Ra[p, g:g+1], Ca[p]]).astype(np.float64)
        reg = np.eye(A.shape[1]) * lam; reg[0, 0] = 0
        coef[g] = np.linalg.solve(A.T @ A + reg, A.T @ Ya[p, g])
    return coef
def predict(m, Xe, qe, Re, Ce, coef):
    qh = np.clip(qe + (Xe - m["xm"]) @ m["W"] + m["c"], 0.001, 0.999)
    if coef is None: mh = np.maximum(m["a"] + m["b"] * Re, 0.01)
    else: mh = np.maximum(coef[:, 0] + coef[:, 1] * Re + Ce @ coef[:, 2:].T, 0.01)
    return (qh * mh).astype(np.float32)
tree = cKDTree(d["xy_um"]); pcs = prediction_pcs(P)
for seed in (0, 1):
    c, r, aidx, ag = place_cores(d["xy_um"], 1.0, 4, "diverse", seed, P=P, pcs=pcs, tree=tree)
    eidx = np.flatnonzero(evaluation_mask(d["xy_um"], c, r, 200)); Ye = Y[eidx]
    res = {}
    for occ_ct in (False, True):
        m = hac_fit(X(aidx, occ_ct), q[aidx], R[aidx], Y[aidx], ag); Xe = X(eidx, occ_ct)
        coef = mag_ct(Y[aidx], R[aidx], C[aidx])
        for mag in ("gene", "gene+celltype"):
            Q = predict(m, Xe, q[eidx], R[eidx], C[eidx], None if mag == "gene" else coef)
            g = from_sums(gene_sums(Ye, Q)); res[f"occ{'+ct' if occ_ct else ''}/mag:{mag}"] = (np.nanmedian(g["pcc"]), g["rmse"].mean(), np.nanmedian(g["rmse_pos"]))
    gp = from_sums(gene_sums(Ye, P[eidx]))
    print(f"seed {seed}: GHIST+ export PCC {np.nanmedian(gp['pcc']):.4f} RMSEpos {np.nanmedian(gp['rmse_pos']):.3f}")
    for k, (a, b, c_) in res.items(): print(f"   {k:28s} medPCC {a:.4f}  RMSE {b:.4f}  RMSE(expressing) {c_:.3f}", flush=True)
