"""Summarise P2 (config_v5.json) and P3 (config_v7.json). Usage: python p23_report.py -> outputs_v5/summary.md, outputs_v7/summary.md"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SLIDES = [("breast5k", "Breast5K (dev)"), ("atera", "ATERA (confirmation)")]


def p2():
    L = ["# P2: support-aware trust gate (config_v5.json, incl. amendments 1–2)", "",
         "Median per-gene PCC (all genes) on evaluation cells, mean of 3 repeats. OUT = cells outside the cores' support (λ < 1).", ""]
    verdict = {}
    for s, lab in SLIDES:
        p = HERE / "outputs_v5" / s / "runs.csv"
        if not p.exists():
            continue
        r = pd.read_csv(p)
        L += [f"## {lab}", "", "| setting | out-of-support cells | OUT: Soft / RRC / HAC / HAC-trust | CI(trust−HAC) > 0 on OUT | ALL: HAC / HAC-trust | HAC worse than Soft on OUT? |", "|---|---|---|---|---|---|"]
        for st in ("unbiased", "tumour_biased", "stroma_biased"):
            o = r[(r.setting == st) & (r.subset == "out_of_support")]
            al = r[(r.setting == st) & (r.subset == "all")]
            if o.empty:
                continue
            m = lambda df, k: df[f"{k}|median_pcc"].mean()
            ci_pos = int((o["ci_HAC-trust_vs_HAC_lo"] > 0).sum())
            problem = m(o, "HAC") < m(o, "GHIST+-Soft")
            L.append(f"| {st} | {100*al.frac_out_of_support.mean():.0f}% | {m(o,'GHIST+-Soft'):.4f} / {m(o,'RRC-authors'):.4f} / {m(o,'HAC'):.4f} / {m(o,'HAC-trust'):.4f} | "
                     f"{ci_pos}/{len(o)} | {m(al,'HAC'):.4f} / {m(al,'HAC-trust'):.4f} | {'yes' if problem else 'no'} |")
            verdict[(s, st)] = dict(problem=problem, ci_pos=ci_pos, all_drop=m(al, "HAC") - m(al, "HAC-trust"))
        L.append("")
    L += ["## Pre-declared decision", ""]
    ok = True
    for s, lab in SLIDES:
        b = [verdict.get((s, st)) for st in ("tumour_biased", "stroma_biased")]
        u = verdict.get((s, "unbiased"))
        if not u or None in b:
            continue
        problem = any(x["problem"] for x in b)
        helps = all(x["ci_pos"] >= 2 for x in b if x["problem"])
        neutral = u["all_drop"] <= 0.002
        ok &= problem and helps and neutral
        L.append(f"- {lab}: problem exists in a biased setting = {problem}; gate fixes it where it exists (CI>0 in ≥2/3) = {helps if problem else 'n/a'}; "
                 f"unbiased all-cell cost {u['all_drop']:+.4f} (allowed ≤ 0.002) → {'PASS' if problem and helps and neutral else 'FAIL'}")
    L += ["", f"**Overall: P2 {'is useful' if ok else 'does not meet the bar'} on these slides.**", ""]
    (HERE / "outputs_v5" / "summary.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def p3():
    L = ["# P3: uncertainty (config_v7.json, incl. amendment 1)", ""]
    res = {}
    for s, lab in SLIDES:
        p = HERE / "outputs_v7" / s / "runs.csv"
        if not p.exists():
            continue
        r = pd.read_csv(p)
        L += [f"## {lab}", "", "**Marginal per-gene intervals (original rule):** median coverage (all entries / expressing entries), mean width", "",
              "| setting | GHIST+ | GHIST+-Soft | RRC-authors | HAC | HAC-mondrian |", "|---|---|---|---|---|---|"]
        for st in ("unbiased", "tumour_biased"):
            t = r[(r.setting == st) & (r.method != "HURDLE-UNCERTAINTY")].groupby("method")[["coverage_all", "coverage_expressing", "mean_width"]].mean()
            if t.empty:
                continue
            L.append(f"| {st} | " + " | ".join(f"{t.loc[m,'coverage_all']:.3f} / {t.loc[m,'coverage_expressing']:.3f}, w {t.loc[m,'mean_width']:.3f}"
                                              for m in ["GHIST+", "GHIST+-Soft", "RRC-authors", "HAC", "HAC-mondrian"]) + " |")
        L += ["", "**Hurdle-aware uncertainty (amendment 1):**", "",
              "| setting | ECE raw q → q+c → HAC q̂ | median Brier raw → HAC | magnitude interval coverage (expressing) / width: HAC vs raw relu(μ) |", "|---|---|---|---|"]
        for st in ("unbiased", "tumour_biased"):
            h = r[(r.setting == st) & (r.method == "HURDLE-UNCERTAINTY")].mean(numeric_only=True)
            if h.empty or np.isnan(h.get("ece_raw_q", np.nan)):
                continue
            L.append(f"| {st} | {h.ece_raw_q:.3f} → {h.ece_q_plus_c:.3f} → **{h.ece_hac_qhat:.3f}** | {h.brier_raw_q:.4f} → **{h.brier_hac_qhat:.4f}** | "
                     f"**{h.mag_cov_hac:.3f} / {h.mag_width_hac:.3f}** vs {h.mag_cov_raw:.3f} / {h.mag_width_raw:.3f} |")
            res[(s, st)] = h
        L.append("")
    L += ["## Pre-declared decisions", ""]
    for s, lab in SLIDES:
        h = res.get((s, "unbiased"))
        if h is None:
            continue
        ok = h.ece_hac_qhat < h.ece_raw_q and h.brier_hac_qhat < h.brier_raw_q and 0.85 <= h.mag_cov_hac <= 0.95
        L.append(f"- {lab}, hurdle-aware rule: ECE {h.ece_raw_q:.3f}→{h.ece_hac_qhat:.3f}, Brier {h.brier_raw_q:.4f}→{h.brier_hac_qhat:.4f}, magnitude coverage {h.mag_cov_hac:.3f} → **{'PASS' if ok else 'FAIL'}**")
        r = pd.read_csv(HERE / "outputs_v7" / s / "runs.csv")
        t = r[(r.setting == "unbiased") & (r.method != "HURDLE-UNCERTAINTY")].groupby("method")[["coverage_all", "coverage_expressing", "mean_width"]].mean()
        best = min(("HAC", "HAC-mondrian"), key=lambda m: t.loc[m, "mean_width"])
        okm = 0.87 <= t.loc[best, "coverage_all"] <= 0.93 and t.loc[best, "mean_width"] < min(t.loc["GHIST+", "mean_width"], t.loc["RRC-authors", "mean_width"])
        L.append(f"- {lab}, original marginal rule: {best} coverage {t.loc[best,'coverage_all']:.3f}, width {t.loc[best,'mean_width']:.3f} → {'PASS' if okm else 'FAIL'} "
                 f"(but expressing-entry coverage {t.loc[best,'coverage_expressing']:.3f}: marginal intervals are not informative for expressing cells)")
    (HERE / "outputs_v7" / "summary.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    p2(); print(); p3()
