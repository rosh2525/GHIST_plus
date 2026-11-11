"""Summarise v8 (fixes 2/3) and v9 (BioHAN) runs against their frozen decision rules.
Usage: python v89_summary.py [--root DIR]  -> prints tables, writes outputs_v9/summary.md
"""
import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
SLIDES = ["breast5k", "atera", "breast5k_pan"]
COLS = ["median_pcc_all", "median_pcc_svg20", "median_ssim_svg50", "cmd_all", "mean_rmse", "median_rmse_pos"]


def table(version, methods, ref):
    out, verdict = [], []
    for s in SLIDES:
        p = HERE / f"outputs_{version}" / s / "runs.csv"
        if not p.exists():
            continue
        r = pd.read_csv(p)
        m = r.groupby("method")[COLS].mean().loc[[x for x in methods if x in set(r.method)]]
        m.insert(0, "slide", s); m.insert(1, "runs", r.config.nunique())
        out.append(m)
        for meth in methods:
            if meth == ref or meth not in set(r.method):
                continue
            rr = r[r.method == meth]
            lo, hi = f"ci_vs_{ref}_all_lo", f"ci_vs_{ref}_all_hi"
            if lo in rr:
                verdict.append({"slide": s, "method": meth, "vs": ref, "runs": len(rr),
                                "CI>0": int((rr[lo] > 0).sum()), "CI<0": int((rr[hi] < 0).sum()),
                                "dPCC": rr.median_pcc_all.mean() - r[r.method == ref].median_pcc_all.mean(),
                                "dSSIM%": 100 * (rr.median_ssim_svg50.mean() / r[r.method == ref].median_ssim_svg50.mean() - 1),
                                "dRMSE%": 100 * (rr.mean_rmse.mean() / r[r.method == ref].mean_rmse.mean() - 1),
                                "dCMD%": 100 * (rr.cmd_all.mean() / r[r.method == ref].cmd_all.mean() - 1)})
    return (pd.concat(out) if out else pd.DataFrame()), pd.DataFrame(verdict)


def main():
    ap = argparse.ArgumentParser(); ap.parse_args()
    md = []
    # v9_ablate: a component "contributes" if removing it lowers PCC vs full BioHAN (CI<0) in the majority of runs
    # on both breast5k and atera (config_v9.json, ablations rule). breast5k_pan is reported but not part of the rule.
    for version, methods, ref in [("v8", ["GHIST+-Soft", "Niche-Soft", "RRC-soft", "Niche-RRC", "Prog-RRC", "Niche-Prog-RRC"], "RRC-soft"),
                                  ("v9", ["GHIST+-Soft", "RRC-soft", "HAC", "BioHAN", "BioHAN-noniche"], "RRC-soft"),
                                  ("v9_ablate", ["RRC-soft", "BioHAN", "BioHAN-nobottleneck", "BioHAN-noweight", "BioHAN-noniche"], "BioHAN")]:
        t, v = table(version, methods, ref)
        for name, df in ((f"{version} means", t), (f"{version} vs {ref}", v)):
            s = df.to_string(float_format=lambda x: f"{x:.4f}") if len(df) else "(none yet)"
            print(f"\n## {name}\n{s}")
            md.append(f"## {name}\n\n```\n{s}\n```\n")
    (HERE / "outputs_v9").mkdir(exist_ok=True)
    (HERE / "outputs_v9" / "summary.md").write_text("\n".join(md))


if __name__ == "__main__":
    main()
