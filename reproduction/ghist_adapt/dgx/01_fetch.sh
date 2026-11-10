#!/usr/bin/env bash
# Download what the DGX runs need from the public GHIST+ bundle (pinned revision) into $GHIST_BUNDLE_ROOT.
#   P4a (second-checkpoint confirmation): pan-cancer Breast5K export (8.2 GB) + Breast5K measured RNA/coords (5.8 GB)
#   P1/P2 reruns: breast_multi Breast5K + ATERA exports and measured data
#   Optional: released checkpoints + configs (only needed if you later run GHIST+ inference yourself)
set -euo pipefail
: "${GHIST_BUNDLE_ROOT:?export GHIST_BUNDLE_ROOT=/scratch/$USER/ghist_bundle first}"
REV=37a937677c5a8ebe5c7b8fbd486fbc5bb50ce1de
FILES=(
  GHIST_plus/predictions/pancancer/breast_5k/test_slide11_pred_expr_scaled.npz
  GHIST_plus/predictions/breast_multi/breast5k/test_slide6_pred_expr_scaled.npz
  GHIST_plus/predictions/breast_multi/atera/test_slide9101_pred_expr_scaled.npz
  GHIST_plus/predictions/breast_multi/test_slide3_pred_expr_scaled.csv
  evaluation_data/breast5k/cell_gene_matrix_filtered.csv
  evaluation_data/breast5k/cell_coords_histology_filtered.csv
  evaluation_data/atera/cell_gene_matrix_filtered.csv
  evaluation_data/atera/cell_coords.csv
  evaluation_data/breast2/cell_gene_matrix_filtered.csv
  evaluation_data/breast2/cell_coords.csv
  evaluation_data/breast2/cell_type_filtered.csv
)
if [[ "${WITH_CHECKPOINTS:-0}" == "1" ]]; then
  for m in breast_multi pancancer; do
    FILES+=(GHIST_plus/models/$m/genes.txt GHIST_plus/models/$m/standardisation_hist_fold_1.npy)
  done
  FILES+=(GHIST_plus/models/breast_multi/ghist_plus_breast_multi_checkpoint.pth GHIST_plus/models/breast_multi/config_multi_breast.json
          GHIST_plus/models/breast_multi/config_inference_breast5k_domain0.json GHIST_plus/models/breast_multi/config_inference_atera_domain0.json
          GHIST_plus/models/pancancer/ghist_plus_pancancer_checkpoint.pth GHIST_plus/models/pancancer/config.json)
fi
huggingface-cli download SydneyBioX/GHIST-Plus-bundle "${FILES[@]}" --repo-type dataset --revision "$REV" --local-dir "$GHIST_BUNDLE_ROOT"
export GHIST_V3_FETCH="$GHIST_BUNDLE_ROOT"   # BreastCancer2 files live in the same tree on the DGX
du -sh "$GHIST_BUNDLE_ROOT"
echo "Done. Also export GHIST_V3_FETCH=$GHIST_BUNDLE_ROOT and GHIST_ADAPT_CACHE=/scratch/\$USER/ghist_adapt_cache before running."
