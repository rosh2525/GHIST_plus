#!/usr/bin/env bash
# P4a: does the zero-shot soft expectation and HAC (and P2's gate) hold for a DIFFERENT released checkpoint?
# Same Breast5K slide, predicted by the pan-cancer GHIST+ model (Breast5K is its held-out test slide 11).
# Protocols unchanged: config_v4.json (P0/P1) and config_v5.json (P2). CPU-bound; run on a CPU node or alongside GPU jobs.
#SBATCH --job-name=ghist_p4a
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=06:00:00
set -euo pipefail
cd "$(dirname "$0")/.."
: "${GHIST_BUNDLE_ROOT:?}"; export GHIST_V3_FETCH=${GHIST_V3_FETCH:-$GHIST_BUNDLE_ROOT}
export GHIST_ADAPT_CACHE=${GHIST_ADAPT_CACHE:-/scratch/$USER/ghist_adapt_cache}
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
python slides.py breast5k_pan
python hurdle.py breast5k_pan
python p1_run.py breast5k_pan > p4a_p1.log 2>&1
python p2_run.py breast5k_pan > p4a_p2.log 2>&1
python - <<'PY'
import json, pandas as pd
z = json.load(open("outputs_v4/breast5k_pan/zero_shot.json"))
print("zero-shot all cells:", {k: round(v["median_pcc_all"], 4) for k, v in z.items()})
r = pd.read_csv("outputs_v4/breast5k_pan/runs.csv")
print(r.groupby(["config", "method"]).median_pcc_all.mean().unstack().round(4))
PY
