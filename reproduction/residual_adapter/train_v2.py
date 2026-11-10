"""v2 training: calibration, ridge and an MLP residual adapter (Apple Metal / MPS) over a
learning curve of local training-set sizes. Uses train + validation cells only.

Usage:  python train_v2.py          (needs torch; see README)
Writes  outputs_v2/seed{s}/n{N}/{calibration.csv, ridge_weights.csv.gz, scaler.csv, fit.json, mlp_curve.csv}
        local_artifacts/adapter_v2_models/seed{s}_n{N}_mlp.pt
"""
from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd
import torch

from common import load_cache, load_split, pcc_per_gene, rmse_per_gene, rows_for
from train_adapter import apply_calibration, fit_calibration, ridge_path
from v2_common import MODEL_DIR_V2, OUT_V2, build_pool, features, load_config_v2, load_extras, make_mlp

BLOCK_UM = 400.0


def torch_median_pcc(Y: torch.Tensor, Q: torch.Tensor) -> float:
    Yc, Qc = Y - Y.mean(0), Q - Q.mean(0)
    den = torch.sqrt((Yc * Yc).sum(0) * (Qc * Qc).sum(0))
    pcc = (Yc * Qc).sum(0) / den
    pcc = pcc[den > 0]
    return float(pcc.median().cpu())


def train_mlp(cfg, Xtr, Rtr, Xva, Pva, Yva, seed, device):
    m = cfg["methods"]["mlp"]
    torch.manual_seed(1000 + seed)
    net = make_mlp(Xtr.shape[1], Rtr.shape[1], m["hidden"], m["dropout"]).to(device)
    opt = torch.optim.AdamW(net.parameters(), lr=m["lr"], weight_decay=m["weight_decay"])
    Xtr_t, Rtr_t = torch.from_numpy(Xtr).to(device), torch.from_numpy(Rtr).to(device)
    Xva_t, Pva_t, Yva_t = (torch.from_numpy(a).to(device) for a in (Xva, Pva, Yva))
    gen = torch.Generator(device="cpu").manual_seed(1000 + seed)
    n, bs = len(Xtr_t), min(m["batch_size"], len(Xtr_t))
    best, best_state, bad, curve, step = -np.inf, None, 0, [], 0
    perm, pos = torch.randperm(n, generator=gen), 0
    while step < m["max_steps"]:
        if pos + bs > n:
            perm, pos = torch.randperm(n, generator=gen), 0
        idx = perm[pos:pos + bs].to(device)
        pos += bs
        net.train()
        loss = torch.mean((net(Xtr_t[idx]) - Rtr_t[idx]) ** 2)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        step += 1
        if step % m["eval_every_steps"] == 0:
            net.eval()
            with torch.no_grad():
                Q = torch.clamp(Pva_t + net(Xva_t), min=0)
                score = torch_median_pcc(Yva_t, Q)
            curve.append({"step": step, "train_loss": float(loss.detach().cpu()), "val_median_pcc": score})
            if score > best:
                best, bad = score, 0
                best_state = {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}
            else:
                bad += 1
                if bad >= m["early_stopping_patience_evals"]:
                    break
    net.load_state_dict(best_state)
    net.eval()
    return net, pd.DataFrame(curve)


def summarise(Y, Q):
    pcc = pcc_per_gene(Y, Q)
    return {"median_pcc": float(np.nanmedian(pcc)), "mean_pcc": float(np.nanmean(pcc)), "mean_rmse": float(rmse_per_gene(Y, Q).mean())}


def main() -> None:
    cfg = load_config_v2()
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    print(f"device: {device}  torch {torch.__version__}")
    d = load_cache()
    ex = load_extras(d["ids"], d["genes"])
    ids, genes = d["ids"], d["genes"]
    MODEL_DIR_V2.mkdir(parents=True, exist_ok=True)
    for seed in (0, 1, 2):
        sp = load_split(seed)
        va = rows_for(ids, np.sort(sp.loc[sp.split == "val", "cell_id"].to_numpy()))
        pool = build_pool(d, seed, BLOCK_UM)
        Pva, Xva = features(d, ex, va)
        Yva = d["gt"][va].astype(np.float32)
        for N in cfg["train_sizes"]:
            t0 = time.time()
            tr = rows_for(ids, np.sort(pool.loc[pool["rank"] < N, "cell_id"].to_numpy()))
            Ptr, Xtr = features(d, ex, tr)
            Ytr = d["gt"][tr].astype(np.float32)
            mu, sd = Xtr.mean(0), Xtr.std(0)
            sd[sd == 0] = 1.0
            Ztr, Zva = ((Xtr - mu) / sd).astype(np.float32), ((Xva - mu) / sd).astype(np.float32)
            out = OUT_V2 / f"seed{seed}" / f"n{N}"
            out.mkdir(parents=True, exist_ok=True)

            a, b = fit_calibration(Ptr.astype(np.float64), Ytr.astype(np.float64))
            pd.DataFrame({"gene": genes, "a": a, "b": b}).to_csv(out / "calibration.csv", index=False)

            best = None
            ridge_curve = []
            for alpha, W, c in ridge_path(Ztr.astype(np.float64), (Ytr - Ptr).astype(np.float64), cfg["methods"]["ridge"]["alphas"]):
                s = summarise(Yva, np.clip(Pva + c + Zva @ W, 0, None))
                ridge_curve.append({"alpha": alpha, **s})
                if best is None or s["median_pcc"] > best[0]["median_pcc"]:
                    best = (s, alpha, W, c)
            _, alpha, W, c = best
            wdf = pd.DataFrame(np.vstack([c[None, :], W]), columns=genes)
            wdf.to_csv(out / "ridge_weights.csv.gz", index=False, float_format="%.6g")
            pd.DataFrame({"mu": mu, "sd": sd}).to_csv(out / "scaler.csv", index=False)
            pd.DataFrame(ridge_curve).to_csv(out / "ridge_val_curve.csv", index=False)

            net, curve = train_mlp(cfg, Ztr, (Ytr - Ptr).astype(np.float32), Zva, Pva, Yva, seed, device)
            curve.to_csv(out / "mlp_curve.csv", index=False)
            torch.save(net.state_dict(), MODEL_DIR_V2 / f"seed{seed}_n{N}_mlp.pt")
            with torch.no_grad():
                Qmlp = np.clip(Pva + net(torch.from_numpy(Zva).to(device)).cpu().numpy(), 0, None)

            val = {"GHIST+": summarise(Yva, Pva), "calibration": summarise(Yva, apply_calibration(Pva, a, b)),
                   "ridge": summarise(Yva, np.clip(Pva + c + Zva @ W, 0, None)), "mlp": summarise(Yva, Qmlp)}
            best_step = int(curve.loc[curve.val_median_pcc.idxmax(), "step"])
            fit = {"seed": seed, "n_train": int(len(tr)), "ridge_alpha": alpha, "mlp_best_step": best_step,
                   "mlp_steps_run": int(curve.step.max()), "device": str(device), "validation": val,
                   "seconds": round(time.time() - t0, 1)}
            (out / "fit.json").write_text(json.dumps(fit, indent=2))
            print(f"seed {seed} N={N:>6}: val medPCC GHIST+={val['GHIST+']['median_pcc']:.4f} calib={val['calibration']['median_pcc']:.4f} "
                  f"ridge={val['ridge']['median_pcc']:.4f} (a={alpha:g}) mlp={val['mlp']['median_pcc']:.4f} (step {best_step}) "
                  f"[{fit['seconds']}s]", flush=True)


if __name__ == "__main__":
    main()
