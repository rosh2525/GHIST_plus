"""DEV-ONLY probe on Breast5K: reference-programme prior and niche features. Not a confirmatory result."""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.spatial import cKDTree
from slides import load; from hurdle import load_hurdle; from p5_ref import load_ref
from metrics import gene_sums, from_sums, corr_matrix, cmd_from_corr
from anchors import place_cores, evaluation_mask, prediction_pcs
from p1_run import hac_fit, _hac_predict, slide_stats
from run_v3 import CFG
d = load("breast5k"); h = load_hurdle("breast5k", d["ids"]); ref = load_ref()
Y, P, q = d["y"], d["pred"], h["q"]; R = np.maximum(h["mu"], 0); S = (q * R).astype(np.float32)
genes = list(d["genes"]); rg = list(ref["genes"]); shared = [g for g in genes if g in rg]
gi = np.array([genes.index(g) for g in shared]); ri = np.array([rg.index(g) for g in shared])
print(f"shared genes {len(shared)}/{len(genes)}")
Xr = ref["X"][:, ri].astype(np.float64); Xr = (Xr - Xr.mean(0)) / np.where(Xr.std(0) > 0, Xr.std(0), 1)
_, sv, Vt = np.linalg.svd(Xr[np.random.default_rng(0).choice(len(Xr), 40000, replace=False)], full_matrices=False)
mS, sS = slide_stats(S)
Zs = ((S[:, gi] - mS[gi]) / sS[gi]).astype(np.float32)
def med(Q): return np.nanmedian(from_sums(gene_sums(Y, Q))["pcc"])
valid = np.flatnonzero(Y.std(0) > 0); gtc = corr_matrix(Y, valid)
def cmd(Q): return cmd_from_corr(gtc, corr_matrix(Q, valid))
print(f"zero-shot  export {med(P):.4f}  soft {med(S):.4f} | CMD export {cmd(P):.3f} soft {cmd(S):.3f}")
for k in (5, 10, 20, 40, 80):
    V = Vt[:k].T.astype(np.float32)
    Q = S.copy(); Q[:, gi] = np.maximum(mS[gi] + sS[gi] * ((Zs @ V) @ V.T), 0)
    print(f"  ref-projection k={k:3d}: medPCC {med(Q):.4f}  CMD {cmd(Q):.3f}", flush=True)
# cell-type expected profile from the atlas, weighted by GHIST+ cell-type probabilities
ctm = ref["celltype_major"]; Xl = ref["X"][:, ri]
prof = {c: Xl[ctm == c].mean(0) for c in np.unique(ctm)}
cls = ["B-cells", "Myeloid", "STROMAL", "Myeloid", "Cancer Epithelial", "Normal Epithelial", "Plasmablasts", "T-cells", "ALL"]
strom = np.mean([prof[c] for c in ("CAFs", "PVL", "Endothelial")], 0); allm = Xl.mean(0)
M = np.stack([strom if c == "STROMAL" else allm if c == "ALL" else prof[c] for c in cls])
Rexp = (h["ctp"] @ M).astype(np.float32)
mR, sR = slide_stats(Rexp)
print(f"  atlas cell-type profile alone (zero-shot): medPCC on shared genes {np.nanmedian(from_sums(gene_sums(Y[:, gi], Rexp))['pcc']):.4f} vs soft {np.nanmedian(from_sums(gene_sums(Y[:, gi], S[:, gi]))['pcc']):.4f}")
# core-based: HAC vs HAC + biology features (+ niche)
V20 = Vt[:20].T.astype(np.float32); prog = Zs @ V20; rprog = (((Rexp - mR) / sR) @ V20).astype(np.float32)
nn = cKDTree(d["xy_um"]).query(d["xy_um"], k=16)[1][:, 1:]
niche = np.hstack([prog[nn].mean(1), h["ctp"][nn].mean(1)]).astype(np.float32)
lq = np.log(np.clip(q, 1e-4, 1 - 1e-4) / (1 - np.clip(q, 1e-4, 1 - 1e-4))).astype(np.float32)
blocks = {"base": [S, lq], "+prog": [S, lq, prog, rprog], "+prog+niche": [S, lq, prog, rprog, niche]}
stats = {k: [slide_stats(b) for b in v] for k, v in blocks.items()}
tree = cKDTree(d["xy_um"]); pcs = prediction_pcs(P)
for seed in (0, 1):
    c, r, aidx, ag = place_cores(d["xy_um"], 1.0, 4, "diverse", seed, P=P, pcs=pcs, tree=tree)
    eidx = np.flatnonzero(evaluation_mask(d["xy_um"], c, r, 200))
    out = []
    for k, bl in blocks.items():
        X = lambda idx: np.hstack([(b[idx] - s[0]) / s[1] for b, s in zip(bl, stats[k])]).astype(np.float32)
        m = hac_fit(X(aidx), q[aidx], R[aidx], Y[aidx], ag)
        Q = _hac_predict(m, X(eidx), q[eidx], R[eidx])
        g = from_sums(gene_sums(Y[eidx], Q))["pcc"]
        out.append(f"{k} {np.nanmedian(g):.4f}")
    print(f"  seed {seed} HAC: " + " | ".join(out), flush=True)
