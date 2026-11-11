#!/bin/bash -l
#SBATCH --job-name=spill_test
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/data/p26bb0003/roshanBTP/logs/spill_%x_%j.log
# usage: sbatch run_spill.sh breast5k   (then atera, breast5k_pan)
set -e
ROOT=/scratch/data/p26bb0003/roshanBTP
export GHIST_ADAPT_CACHE=$ROOT/cache OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd $ROOT/code/ghist_adapt
date; hostname
$ROOT/env/bin/python -u spillover/spill_test.py "${1:-breast5k}" --out $ROOT/outputs/spillover
date
