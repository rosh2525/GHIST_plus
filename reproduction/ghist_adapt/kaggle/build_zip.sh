#!/usr/bin/env bash
# Builds the code bundle to upload as a private Kaggle Dataset named "ghist-adapt-code".
set -euo pipefail
cd "$(dirname "$0")/../../.."          # GHIST_plus root
OUT=${1:-local_artifacts/ghist_adapt_kaggle.zip}
rm -f "$OUT"
zip -qr "$OUT" model tools/one_shot_rrc.py \
  reproduction/ghist_adapt/*.py reproduction/ghist_adapt/*.json reproduction/ghist_adapt/requirements.txt \
  reproduction/ghist_adapt/dgx reproduction/ghist_adapt/kaggle -x '*/__pycache__/*'
ls -lh "$OUT"
