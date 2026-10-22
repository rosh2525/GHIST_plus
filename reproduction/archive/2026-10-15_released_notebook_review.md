# Archived broad released-notebook review — superseded

This October15 review is outside the current paper-only scope. Use [the current report](../REPORT.md). Its numeric results are later released-code results, not verified manuscript matches.

# GHIST+ reproduction report

Run completed **15 October 2026 (Asia/Kolkata)** in this workspace.

**Status: all four released figure notebooks executed successfully. This is a reproduction of evaluation from released predictions and summary tables. Fresh H&E-to-RNA inference and training have not been performed.**

## What was run

| Notebook | Scope reproduced | Evidence level |
| --- | --- | --- |
| Figure2 | Panels b–f: ROI/pathway maps and correlations | Replotting authors' released CSV summaries; original input checksums passed |
| Figure3 | Panels a–i, including ATERA and Breast5K | Recomputed metrics/maps from released cell-level predictions and measured RNA |
| Figure4 | Panels a, b, c, e and f provided by notebook | Recomputed imputation/structure metrics and VQ/composition ablation from released predictions |
| Figure5 | Three transfer directions | Replotting authors' paired per-gene PCC table; coverage checked; raw transfer predictions not rerun |

The original code cells were executed unchanged. One final cell per notebook preserves temporary plots and tables. No synthetic cells, RNA profiles or images were generated. The published randomized VQ control is an existing released ablation on real cells, not newly generated synthetic data.

## Source and environment

- Official source: https://github.com/SydneyBioX/GHIST_plus
- Source commit: `c5c895174c19903cc27d2661934b9b726bdb9306`.
- Official data: https://huggingface.co/datasets/SydneyBioX/GHIST-Plus-bundle
- Dataset revision: `37a937677c5a8ebe5c7b8fbd486fbc5bb50ce1de`.
- 45 downloaded files; 10.923 GB. Every file was checked against the source LFS SHA-256 or Git blob hash.
- Apple Silicon macOS, 16 GiB RAM; analysis environment uses Python 3.12.14. This differs from the authors' Python 3.10/CUDA training environment.
- Main numeric package versions follow the authors' tested versions. Matplotlib inline backend pinned to 0.1.7 for compatibility with Matplotlib 3.8.2.
- Full package versions: [pip_freeze.txt](../../results/reproduction/pip_freeze.txt). Training seeds/configs remain in the unchanged upstream files.

## Expression prediction results

PCC measures correlation, **not percentage accuracy**. SVGs are genes ranked as spatially variable in measured target RNA. These selected genes do not represent the full gene panel.

| Scope | Gene set | Method | Median PCC | Defined genes | Requested genes |
| --- | --- | --- | --- | --- | --- |
| Within slide | Top 50 SVGs | GHIST+ | 0.6410 | 50 | 50 |
| Within slide | Top 50 SVGs | GHIST | 0.6346 | 50 | 50 |
| Within slide | Top 50 SVGs | SpatialEx | 0.4957 | 50 | 50 |
| Across slide | Top 50 SVGs | GHIST+ | 0.5508 | 50 | 50 |
| Across slide | Top 50 SVGs | GHIST | 0.2759 | 42 | 50 |
| Across slide | Top 50 SVGs | SpatialEx | 0.1101 | 50 | 50 |
| Across slide | Top 50 SVGs | Phoenix | 0.4804 | 50 | 50 |

The within-slide cohort contains 12,719 matched cells and 280 shared genes. The across-slide cohort contains 76,779 cells and 276 shared genes. The numbers above reproduce the notebook's per-method median, excluding undefined correlations. In the across-slide top-50 comparison, GHIST has only 42 defined correlations; GHIST+, SpatialEx and Phoenix each have 50.

For a separate comparison on **exactly the same genes with finite PCC for every method**, the across-slide top-50 results are:

| Gene set | Method | Median PCC | Common defined genes | Requested genes |
| --- | --- | --- | --- | --- |
| Top 50 SVGs | GHIST+ | 0.5659 | 42 | 50 |
| Top 50 SVGs | GHIST | 0.2759 | 42 | 50 |
| Top 50 SVGs | SpatialEx | 0.1058 | 42 | 50 |
| Top 50 SVGs | Phoenix | 0.4859 | 42 | 50 |

This additional coverage diagnostic does not replace or modify the official notebook results. The all-shared-gene across-slide median for GHIST+ is **0.2031** (276/276 defined), well below its selected top-50 median of 0.5508. In the within-slide top-20 SVG group, original GHIST slightly exceeds GHIST+ (0.7746 vs 0.7623). Phoenix also has lower/better correlation-matrix distance in the released across-slide comparison. These results do not support a claim that GHIST+ wins every task.

![Across-slide released-prediction evaluation](../../results/reproduction/across_slide_summary.png)

## External datasets

| Dataset | Gene set | Median PCC | Defined genes | Requested genes | Matched cells |
| --- | --- | --- | --- | --- | --- |
| ATERA | Top 20 SVGs | 0.4707 | 20 | 20 | 168335 |
| ATERA | Top 50 SVGs | 0.3949 | 50 | 50 | 168335 |
| ATERA | All shared genes | 0.1085 | 403 | 413 | 168335 |
| Breast5K | Top 20 SVGs | 0.3826 | 20 | 20 | 564723 |
| Breast5K | Top 50 SVGs | 0.2434 | 50 | 50 | 564723 |
| Breast5K | All shared genes | 0.0704 | 258 | 258 | 564723 |

ATERA has 168,335 unique matched cells and 413 shared genes; 403 genes have defined PCC. Breast5K has 564,723 matched cells and 258 shared genes, all with defined PCC. These are the released breast multi-slide model's external evaluations, not the Figure5 PanCancer comparisons. Counts refer to cells, not independent patients.

## Gene imputation and ablation

Median metrics over **45 shared held-out slide–gene cases across five real slides**:

| Method | PCC | RMSE | SSIM |
| --- | --- | --- | --- |
| GHIST+ | 0.7371 | 0.7848 | 0.7085 |
| PASTA | 0.5387 | 1.4730 | 0.3178 |
| SpatialEx | 0.0901 | 5.0790 | 0.2938 |
| SpatialEx+ | 0.4662 | 0.8568 | 0.6274 |

Higher PCC/SSIM is better; lower RMSE is better. These use the authors' saved prediction transforms and evaluation protocol. The 45 cases include different genes on different slides; they are not 45 independent patients.

The VQ/composition ablation used 13,621 cells and 280 genes. The full-model top-20 median PCC was 0.7586, versus 0.5730 when the VQ/composition pathway was off and 0.5810 for the published randomized control. Top-50 medians were 0.6425, 0.5183 and 0.5122. This evaluates saved frozen-checkpoint interventions; no control was retrained here.

The Figure3 within-slide cohort has 12,719 cells, while this Figure4 ablation uses 13,621. Their scores must not be treated as identical evaluation cohorts.

## PanCancer transfer: coverage matters

Figure5 is generated from a released per-gene summary table. Its paired comparisons have only **5/13** genes for pancreas top-20/top-50, **5/9** for breast, and **2/4** for lung. The union target panels are larger. Plotted pair counts should not be described as 20 or 50 jointly evaluable genes in every direction. Raw per-cell PanCancer predictions for these case studies were not independently re-evaluated by this notebook run.

| Direction | Gene_set | target_n | paired_n | single_median_paired | pancancer_median_paired |
| --- | --- | --- | --- | --- | --- |
| Pancreas 1 → Pancreas 3 | Top 20 | 20 | 5 | 0.2716 | 0.4321 |
| Pancreas 1 → Pancreas 3 | Top 50 | 50 | 13 | 0.1841 | 0.3750 |
| Breast 1 → Breast 5k | Top 20 | 20 | 5 | 0.3341 | 0.4919 |
| Breast 1 → Breast 5k | Top 50 | 50 | 9 | 0.3156 | 0.4729 |
| Lung 1 → Lung 5k | Top 20 | 20 | 2 | 0.5092 | 0.8026 |
| Lung 1 → Lung 5k | Top 50 | 50 | 4 | 0.6146 | 0.7365 |

## Verification

- Every original source code cell matches the executed notebook exactly.
- All four notebooks completed with no execution errors; 16 displayed plot outputs were saved.
- The notebooks' cohort, gene-coverage and input assertions passed.
- Independent SciPy PCC checks recomputed three SVG-ranked genes across four methods from the real count/prediction files on 76,779 matched cells. Maximum difference from saved notebook PCC: 5.34e-13.
- External cohort checks passed for ATERA and Breast5K; the five-slide imputation check passed.
- [verification.json](../../results/reproduction/verification.json) contains the check receipts; [run_manifest.json](../../results/reproduction/run_manifest.json) contains timings and source hashes. This verifies execution/numeric consistency, not independent model training or exact correspondence to every value in a manuscript version.

## Remaining requirements for a full model reproduction

The official figure bundle does not contain the aligned training H&E images, nuclei masks and matched-nuclei files required by the training configuration. Preflight found 31 missing path fields for `config_multi_breast.json` (some fields are related inputs for the same slide).

UNI2-H is a gated encoder. The GHIST+ released checkpoints exclude its tensors; its official weights require an approved Hugging Face account and institutional-email access: https://huggingface.co/MahmoodLab/UNI2-h . No encoder weights were substituted and no access application was submitted. The current analysis environment has no PyTorch/CUDA installation; this Mac has no NVIDIA CUDA GPU. The authors' tested training setup uses Linux/CUDA.

Consequently, **new image predictions, released-checkpoint inference, fresh training, baseline retraining, and patient-level replication remain unfinished**. The original training/inference code and configurations are cloned locally, and the readiness check names each missing input. A full rerun must retain the published slide splits and settings; a reduced Mac run would be a different experiment.

## Open the results

- [Figure2 executed notebook](../../results/reproduction/Figure2/Figure2.executed.ipynb)
- [Figure3 executed notebook](../../results/reproduction/Figure3/Figure3.executed.ipynb)
- [Figure4 executed notebook](../../results/reproduction/Figure4/Figure4.executed.ipynb)
- [Figure5 executed notebook](../../results/reproduction/Figure5/Figure5.executed.ipynb)
- [Figure3 across-slide plot](../../results/reproduction/Figure3/3c.png)
- [Figure4 ablation plot](../../results/reproduction/Figure4/figure4f_vq_ablation.png)
- [Figure5 transfer plot](../../results/reproduction/Figure5/Figure5.png)

The original Figure4 ablation strip titles overlap slightly in the authors' plotting layout; its numerical data are intact. The separate summary figure above has a readable layout without altering the original notebook.
