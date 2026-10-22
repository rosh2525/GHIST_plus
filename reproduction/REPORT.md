# GHIST+ manuscript-only reproduction

Scope updated **22 October 2026**. Source: *Generalisable tissue-wide molecular reconstruction from histology*, bioRxiv v1, posted 12 June 2026, DOI 10.64898/2026.06.09.731252.

**Exact numerical reproduction of the manuscript has not been achieved.** The completed local run reproduces the authors' later released notebooks. Those notebooks differ from the manuscript in reported values, panel numbering, normalization and configuration. A successful notebook run does not establish a successful paper reproduction.

## Metrics defined in this paper

| Metric | Paper page | Meaning |
| --- | --- | --- |
| PCC | 22 | one gene across matched cells |
| Spearman correlation | 22,24 | gene across cells, or pathway across spatial regions |
| RMSE | 22 | one gene across matched cells; manuscript size-factor log normalization |
| SSIM | 23 | binned/smoothed spatial expression map for one gene |
| CMD | 23 | predicted versus observed gene-gene correlation matrices |
| TISSUE score | 23 | local-variability-normalized cell-gene error; signed definition and displayed aggregation must be resolved |
| AUROC | 23 | one cell type versus the remaining cell types using pyUCell signatures |
| Recall | 23,24 | class-average recall; defined in Methods, no exact main-text result transcribed |
| Spatial IoU | 24 | per-cell-type spatial density maps; 150 micrometre bins and sigma 1.2 bins |

SSIM and AUROC should not be called percent accuracy. TISSUE is an error/variability score, not confidence-interval coverage. The paper also presents downstream cell proportions, ECM scores, regression slopes and P-values for GTEx; these are not single-cell RNA prediction accuracy.

## Values actually reported by the paper

The following are manuscript targets, **not new measurements from our machine**. Values come from the Results text on PDF pages 5–6. The full transcription, including pathway results, is in `paper_reported_values.csv`.

| experiment | method | metric | cohort | value | pdf_page |
| --- | --- | --- | --- | --- | --- |
| Within slide | GHIST+ | PCC | Top 20 SVGs | 0.7440 | 5 |
| Within slide | GHIST | PCC | Top 20 SVGs | 0.7140 | 5 |
| Within slide | SpatialEx | PCC | Top 20 SVGs | 0.5210 | 5 |
| Across slide | GHIST+ | PCC | Independent slide headline | 0.6530 | 5 |
| Imputation | GHIST+ | PCC | Held-out genes result 1 | 0.7430 | 5 |
| Imputation | GHIST+ | PCC | Held-out genes result 2 | 0.7470 | 5 |
| Imputation | GHIST+ | RMSE | Held-out genes result 1 | 2.9040 | 5 |
| Imputation | GHIST+ | RMSE | Held-out genes result 2 | 3.0100 | 5 |
| Imputation | GHIST+ | TISSUE score | Held-out genes | 0.3430 | 6 |
| PancreasADC2 to PancreasADC3 | Single slide | PCC | Top 20 SVGs | 0.0030 | 6 |
| PancreasADC2 to PancreasADC3 | Pan-cancer assisted | PCC | Top 20 SVGs | 0.1020 | 6 |
| PancreasADC2 to PancreasADC3 | Single slide | PCC | Top 50 SVGs | 0.0030 | 6 |
| PancreasADC2 to PancreasADC3 | Pan-cancer assisted | PCC | Top 50 SVGs | 0.1720 | 6 |
| PancreasADC to PancreasADC3 | Single slide | PCC | Top 20 SVGs | -0.0010 | 6 |
| PancreasADC to PancreasADC3 | Pan-cancer assisted | PCC | Top 20 SVGs | 0.2360 | 6 |

The manuscript lists two PCC and two RMSE imputation values together. Their complete protocol assignments are not specified in that sentence; they are retained as separate reported results rather than assigning an invented average or cohort. SSIM, CMD, AUROC, recall and IoU are included in the metric specification; no precise scalar was invented for a boxplot without a reported numerical value.

## Cell-type values printed on the paper figure

These GHIST+ values are rounded labels from Figure3e/f (PDF page30), not results obtained by training here. The corresponding file retains their figure provenance.

| Cell type | AUROC | Spatial IoU |
| --- | --- | --- |
| B | 0.7000 | 0.2500 |
| Epithelial | 0.9300 | 0.7100 |
| Macrophage | 0.7300 | 0.5400 |
| Malignant | 0.9500 | 0.6500 |
| Myeloid excluding macrophage | 0.7700 | 0.3400 |
| Plasma | 0.8700 | 0.4400 |
| Stromal | 0.8500 | 0.5200 |
| T | 0.8900 | 0.8000 |

## Current released-code results differ

| Experiment | Method | Paper-reported PCC | Released-code PCC | Status |
| --- | --- | --- | --- | --- |
| Within slide, top 20 SVGs | GHIST+ | 0.7440 | 0.7623 | Does not reproduce manuscript value |
| Within slide, top 20 SVGs | GHIST | 0.7140 | 0.7746 | Does not reproduce manuscript value |
| Within slide, top 20 SVGs | SpatialEx | 0.5210 | 0.6377 | Does not reproduce manuscript value |

The current released within-slide cohort has 12,719 matched cells. The paper describes five 80/20 spatial folds, while the checked-in config supplies only one 0–20% validation division. These differences have not been reconciled, so the released values cannot be labelled as matched manuscript results. We do not tune or replace observed results just to obtain a target number.

## Figure correspondence

| Paper panels | Current counterpart | Status |
| --- | --- | --- |
| 2b-f | Figure2.ipynb | Paper-relevant plots exist; numerical mismatch prevents an exact reproduction claim. |
| 3a,b,c | Figure3 tables a,b,c | Current released cohort; three paper methods retained; five-fold paper result not reconstructed. |
| 3d,e,f | Figure3 tables f,d,e respectively | CMD/AUROC/IoU remapped to manuscript numbering; results remain released-code results. |
| 4a,b,d | Figure4.ipynb | Available partial outputs; paper normalization, eligibility and complete baselines unresolved. |
| 4c | Later Figure4 table f | Current frozen-checkpoint intervention has not been established as the manuscript no-VQ experiment. |
| 5a | Current Figure5.ipynb is different | Paper-specific PancreasADC and PancreasADC2 transfer has not been reproduced. |
| 5b-f | Figure5_new.Rmd | GTEx sample objects, compartments and molecular predictions are missing; analyses not run. |

The main paper compares GHIST+, GHIST and SpatialEx for expression prediction. Phoenix, ATERA/Breast5K external benchmarks, the extra common-finite-gene analysis and the SpatialEx+ arm are excluded from this manuscript-only report. No new architecture or improvement experiment is introduced. Earlier broader output is archived separately.

## Specific issues preventing a full match

1. **Normalization.** Manuscript RMSE uses library-size-normalized log expression (Eq.15). The released imputation notebook uses log1p of raw counts. Its previous RMSE is therefore excluded from the manuscript result table.
2. **Folds and settings.** Paper: five folds, learning rate 0.001, weight decay 0.0001. Released config: one division, learning rate 0.00115, weight decay 0.00001. Templates restore the two stated optimizer settings; they do not invent the missing original fold assignments.
3. **Panel selection.** The imputation release evaluates 45 common slide–gene cases; the paper describes masking the top ten genes on each of five slides. Original exclusions and complete measured panels are needed before claiming the same cohort.
4. **Variants.** The later VQ/composition intervention and manuscript ablation are not established as the same experiment. Baseline outputs for SpaGE, SpaIM, stImpute and Tangram described in the text are also absent locally.
5. **Figure5.** The manuscript includes pancreas transfer and 77-sample normal-breast/age analyses. The current Python notebook's three transfer directions do not reproduce that full figure. The R script refers to unavailable author-side objects.
6. **Model execution.** Aligned H&E/masks/matched inputs, approved UNI2-H weights and a suitable model environment are still needed. No fresh model training or image inference has run.

The paper itself has caption inconsistencies: Figure3b plots HVGs but its caption says non-SVGs; Figure4c plots non-SVGs but its caption says HVGs. These are recorded rather than silently equating HVGs with non-SVGs. The Methods define a signed TISSUE error while the released notebook aggregates absolute errors; that convention also needs reconciliation.

## Workflow

Read [WORKFLOW.md](WORKFLOW.md) for the exact data → preprocessing → model → experiments → metrics sequence and commands. The original PDF is preserved unmodified at `local_artifacts/paper/GHIST_plus_2026_06_09_731252v1.pdf` with SHA-256 `e67d6d6db08144737c75a553d96ec875269e98224e0c157eb1c1681839abc938`.

Execution receipts and source hashes are in `paper_scope_receipt.json`. All exported numeric rows use metrics mentioned in this manuscript; every source is labelled as manuscript text or later released-code evaluation.
