#!/usr/bin/env bash
# Run the v3 grid on a DGX/HPC node. v3 is CPU-bound (closed-form ridge + metrics), so it needs CPU cores and RAM,
# not a GPU: one slide per process, in parallel. Adjust the scheduler header to your cluster, e.g. for SLURM:
#SBATCH --job-name=ghist_adapt_v3
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=04:00:00
set -euo pipefail
cd "$(dirname "$0")/.."
export GHIST_BUNDLE_ROOT=${GHIST_BUNDLE_ROOT:?set to the GHIST-Plus-bundle folder}
export GHIST_V3_FETCH=${GHIST_V3_FETCH:-$GHIST_BUNDLE_ROOT}       # breast2 files, if downloaded elsewhere
export GHIST_ADAPT_CACHE=${GHIST_ADAPT_CACHE:-$PWD/../../local_artifacts/ghist_adapt_cache}
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
python -m pip install -q -r requirements.txt
for s in breast5k atera breast2; do
  python run_v3.py "$s" > "run_${s}.log" 2>&1 &
done
wait
python report_v3.py
