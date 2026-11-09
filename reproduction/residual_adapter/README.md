# Residual adapter on frozen GHIST+ predictions (Breast5K, within-slide)

**Question.** Can a small amount of measured RNA from a few tissue regions train a lightweight model that improves frozen GHIST+ predictions in other regions of the same slide?

**Answer on this run: no.** The multigene residual ridge adapter failed the pre-specified validation gate on all three spatial splits, so it was not developed further. On the locked test regions it *lowered* median per-gene PCC from 0.071 to 0.066 and reduced PCC for about 106 of 258 genes. Per-gene calibration cut RMSE by about 3×, but mostly by shrinking predictions toward zero (predicted/measured SD ratio 0.065). That is the near-zero artefact the plan said to check for, not a real improvement. This is a negative result for one donor and one slide.

This is supervised adaptation with extra local RNA, not GHIST+ retraining. GHIST+ itself was not re-run: the authors' released predictions (`test_slide6_pred_expr_scaled.npz`) are used unchanged.

## Data

| Item | Source (SydneyBioX/GHIST-Plus-bundle, rev `37a9376`) |
|---|---|
| Measured counts | `evaluation_data/breast5k/cell_gene_matrix_filtered.csv` (5.76 GB) |
| GHIST+ predictions | `GHIST_plus/predictions/breast_multi/breast5k/test_slide6_pred_expr_scaled.npz` (2.7 GB; the bundle ships `.npz`, not `.csv`) |
| Coordinates | `evaluation_data/breast5k/cell_coords_histology_filtered.csv` (22 MB) |

These were already downloaded to `local_artifacts/bundle/` by the earlier reproduction run and were reused here. The dataset is the 10x Xenium Prime FFPE human breast cancer sample (one donor).

**Loader check** (`outputs/scale_check.json`). Following `Figure3.ipynb`, the target is log1p(raw counts) and the prediction is `pred_expr_scaled × 0.5`. This gives 564,723 cells × 258 genes, the authors' expected cohort. Full-slide per-gene PCC matches the released `results/reproduction/Figure3/3i.csv` to within 1e-8 (full-slide median PCC 0.070). Measured expression is very sparse: mean log1p is 0.019, against a mean prediction of 0.150.

## Protocol (frozen in `config.json` before any test metric)

- **Splits** (`splits/`): 400 µm blocks, keeping blocks with at least 600 cells. Each seed shuffles the blocks and assigns them greedily to train, validation and test. Blocks in different splits are never adjacent, so every block is separated by at least one 400 µm block. Cells are then subsampled to exactly 5,000 train, 5,000 validation and 20,000 test. The secondary `rest` set is every cell in blocks at least one block away from train and validation (about 480k cells). Seeds 0, 1 and 2 are used. Block choice uses coordinates only. See `outputs/figures/split_layout.png`.
- **Per-gene calibration**: OLS `y_g = a_g + b_g·p_g` on training cells, clipped at 0.
- **Multigene residual adapter**: ridge regression from all 258 standardised predictions to the residual `y − p`. The output is `clip(p + correction, 0)`. Alpha comes from {0.1 … 1e5}, chosen by validation median per-gene PCC.
- **Gate**: the adapter must beat calibration on validation median PCC. If it does not, stop expanding it.
- **Metrics**: per-gene PCC, RMSE (log1p), RMSE on measured-positive cells, predicted/measured SD ratio, and counts of genes improved or degraded (|ΔPCC| > 0.005).

## Results

Mean ± SD over three spatial splits (`outputs/summary.md`, `outputs/summary.csv`):

| set | method | median PCC | mean RMSE | median RMSE, measured>0 | SD ratio | genes PCC ↑/↓ vs GHIST+ |
|---|---|---|---|---|---|---|
| test (20k) | GHIST+ frozen | **0.0712 ± 0.0010** | 0.324 | 0.667 | 2.75 | – |
| test (20k) | per-gene calibration | 0.0708 ± 0.0006 | **0.107** | 0.725 | 0.065 | 17 / 26 |
| test (20k) | residual adapter | 0.0660 ± 0.0019 | 0.202 | **0.649** | 1.81 | 12 / 106 |
| rest (~480k) | GHIST+ frozen | **0.0719** | 0.335 | 0.664 | 3.07 | – |
| rest (~480k) | per-gene calibration | 0.0717 | **0.100** | 0.725 | 0.072 | 10 / 28 |
| rest (~480k) | residual adapter | 0.0645 | 0.207 | **0.644** | 2.01 | 11 / 106 |

Validation gate: **failed for all seeds**. On validation, median PCC for GHIST+ / calibration / adapter was 0.065 / 0.062 / 0.057 (seed 0), with similar values for the other seeds. Alpha was always selected at the top of the grid (1e5), which means validation preferred the least-changed adapter.

### Interpretation

1. **The adapter did not improve gene ranking.** Mixing in information from other genes reduced PCC for far more genes than it helped, both on test and on the roughly 480k remaining cells. Validation predicted this correctly.
2. **The calibration RMSE gain is mostly an artefact.** GHIST+ overestimates magnitude: its predictions have about 3× the measured SD. Calibration fixes this by flattening predictions toward the mean (SD ratio about 0.07). Error on measured-positive cells gets worse (0.667 → 0.725), and 45–53 genes get negative slopes (MKI67 in `outputs/figures/spatial_maps_seed0.png`). This is why the plan required the positive-cell and variance checks.
3. **The adapter has the lowest error on positive cells**, but with worse PCC and higher overall RMSE than calibration. It is a trade-off, not a win.
4. **Probable cause.** The 258-gene Breast5K panel is extremely sparse at the cell level, so most genes' median PCC is near noise (about 0.07). Within one slide, 5,000 local cells carry too little gene–gene signal for a 258×258 linear map to generalise across regions. Strong markers (ERBB2 0.54, EPCAM 0.46 in the mapped block) are already ranked well by GHIST+.

### Limits

These are within-slide results for one donor, not evidence about other patients. The gate failure was decided on validation alone; test numbers were computed once with the frozen config. Train-only fitting was used (validation cells were not refitted in).

## Reproduce

Local (measured on a 16 GB, 8-core Mac with the bundle already present: about 45 s for prepare, 4 s for train, 16 s for evaluate; writes a 1.2 GB cache to `local_artifacts/adapter_cache/`):

```sh
cd reproduction/residual_adapter
PY=../../../work/ghist_plus_reproduction/.venv/bin/python   # numpy, pandas, scipy, matplotlib
$PY prepare_data.py      # match cells/genes, check against 3i.csv, write splits
$PY train_adapter.py     # calibration + ridge; selects alpha on validation
$PY evaluate.py          # locked test + rest metrics, figures, summary.md
```

Set `GHIST_BUNDLE_ROOT` if the bundle is stored elsewhere. Google Colab: `colab_residual_adapter.ipynb` downloads the three files (about 8.5 GB) and runs the same scripts. The comparison against `3i.csv` is skipped there because that file is a local reproduction output.

## Files

`config.json` holds the frozen settings. `splits/` holds the split manifests (cell IDs per split) and `split_summary.json`. `outputs/seed*/` holds the fitted calibration, adapter weights (`adapter_weights.csv.gz`, intercept plus 258×258), scaler, validation curve, fit summary and per-gene metrics. `outputs/figures/` holds the split layout, PCC distributions and spatial maps.
