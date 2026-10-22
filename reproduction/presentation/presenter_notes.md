# GHIST+ workflow and reproduction evidence

22 October 2026. Six slides.

## Slide 1

Purpose and status
This presentation explains the published GHIST+ model, the real datasets used for the local evaluation, commands and saved evidence, and the steps still required for an exact reproduction. GHIST+ is the direct successor to GHIST; this deck is about GHIST+ rather than a new architecture. The target manuscript is Andrew Zhang et al., Generalisable tissue-wide molecular reconstruction from histology, bioRxiv v1, posted 12 June 2026. It is a preprint.

What has actually been done
On 15 October the authors' released Figure2–Figure5 notebooks executed successfully in an isolated analysis environment. Some figures replot authors' summary CSVs, while Figure3 and parts of Figure4 recalculate metrics from measured RNA and released predictions. On 22 October the source/cohort verification and the presentation evidence commands were run again. No model was newly trained and no new H&E image was passed through the complete GHIST+ model locally. Do not say that we independently reproduced every paper experiment or attained the paper's accuracy.

Beginner explanation
H&E is a coloured microscope image of tissue. Xenium measures selected RNA molecules at cell locations. A gene-expression matrix contains one row per cell and one column per gene. GHIST+ learns an association between image appearance and measured RNA during training, then predicts a spatial expression matrix from an image and nucleus masks. Its inferred RNA is not itself a molecular measurement.

Sources
https://doi.org/10.64898/2026.06.09.731252
https://github.com/SydneyBioX/GHIST_plus
https://huggingface.co/datasets/SydneyBioX/GHIST-Plus-bundle
Local evidence: reproduction/receipts/run_manifest.json; reproduction/receipts/verification.json; reproduction/REPORT.md.
Code commit c5c895174c19903cc27d2661934b9b726bdb9306. Data revision 37a937677c5a8ebe5c7b8fbd486fbc5bb50ce1de.

## Slide 2

What data did we use?
Only real public data, obtained from the official GHIST-Plus-bundle. The source is processed human 10x Xenium RNA and the model outputs released by the authors. The selected download has 45 source-verified files totalling 10.922663646 GB, about 10.923 GB in decimal units. This is an evaluation subset, not the entire roughly 38 GB bundle. It does not contain every aligned image, mask and original training split needed for fresh model training. No synthetic RNA, cells or fabricated prediction arrays were used.

Why three BreastCancer2 counts?
86,344 is the number of cells in the processed BreastCancer2 RNA input, across 280 genes. Within-slide comparison uses a smaller set of 12,719 cells shared by the available within-slide outputs and targets, across 280 shared genes. Across-slide comparison has 76,779 matched cells and 276 common genes. These are different experiments and different matching filters; do not add their cell counts together or treat them as separate patients. Gene correlation is computed across matched cells for each gene.

Cohort qualification
The authors' current across-slide notebook forms a common intersection including a later Phoenix output. Our manuscript-scoped export retains GHIST+, GHIST and SpatialEx, the three methods in this manuscript comparison, but preserves the current released notebook's already frozen cohort rather than silently redefining it. Thus the released cohort is not certified as the manuscript's original cohort. The large number of cells does not establish multi-patient generalization. Donor identity is not inferred from slide names.

Imputation
Imputation means predicting genes intentionally withheld from the observed training panel. These are real measured genes retained separately as truth. The available comparison contains 45 slide–gene cases, rather than the intended 50 cases implied by ten genes on five slides. The released aligned sets contain: BreastCancer1 (slide0), 4,229 cells and 10 genes; IDC (slide1), 75,232 cells and 9 genes; ILC (slide2), 19,883 cells and 10 genes; IDC3 (slide4), 17,151 cells and 6 genes; IDC2 (slide5), 46,493 cells and 10 genes. These are evaluation cohorts and not the complete slide cell counts. A slide–gene case means one gene evaluated on one slide. Missing overlap must not be invented or replaced with synthetic measurements.

Paper data versus local data
The manuscript's Table1 describes 15 processed Xenium slides spanning breast, colorectal, lung and pancreas. Its cross-breast model trains on BreastCancer1, IDC, IDC2, ILC and IDC3 and tests BreastCancer2. The pancreas transfer study holds out PancreasADC3. The application to 77 normal female GTEx breast samples is a separate downstream analysis and is not locally reproduced. GTEx does not supply measured single-cell RNA truth for every image nucleus.

Sources
Manuscript pp.15–16, Table1. https://huggingface.co/datasets/SydneyBioX/GHIST-Plus-bundle
Local count source: results/paper_only/presentation_dataset_summary.json generated by reproduction/presentation_evidence.py data. Data/download manifest and verification receipts are retained in the repository.

## Slide 3

Start with the tissue
The training data pair an H&E image with measured Xenium RNA. Registration puts the image and RNA cell coordinates into the same coordinate system. Nucleus instance masks identify individual nuclei: each nucleus has its own label. The workflow follows HoVerNet segmentation and nucleus/RNA correspondence referenced from GHIST, with scClassify/reference-derived cell-type labels. Training pairs require reliable spatial correspondence; a nucleus mask is not an RNA measurement.

Prepare model inputs
The paper uses Macenko stain normalization, per-channel RGB standardization and 224 by 224 image tiles. Training tiles do not overlap; validation/inference uses 30-pixel overlap. Training augmentation includes flips, right-angle rotations and HED colour changes. The GTEx segmentation pipeline has different window sizes; those should not be confused with model tiles.

UNI2h
UNI2h is a pretrained pathology image encoder. Frozen means its weights are kept fixed. It converts image patches into numerical features; cell masks pool features for individual nuclei and tile features represent surrounding tissue. Access to UNI2-H is gated and must follow the official model's approval and license terms. Random encoder weights would not reproduce this model.

ECRM
The cell-context branch passes information between nearby cells, using 12 within-patch and 8 cross-patch neighbours in the paper settings. A cell's neighbourhood can help interpret visually similar nuclei. It does not directly observe a missing gene.

VQ
Vector quantization represents recurring tissue contexts with 64 learned prototypes. It provides tissue-level structure and supports the selection of representative regions in the tissue-core application. A prototype is a learned numerical pattern, not a synthetic biological cell.

Reconstruction and MAGCI
The model reconstructs cell gene expression, with auxiliary information about cell types and tissue composition. MAGCI handles completion across heterogeneous measured gene panels. A gene absent from a measurement panel must be masked rather than treated as an observed zero. This makes it possible to learn from panels with different measured genes, while predictions still require independent validation. The diagram simplifies the full multi-task architecture; see the paper's detailed architecture and losses for implementation.

Training versus deployment
Training consumes measured RNA as a target, together with labels and patch composition. In ordinary zero-shot image inference, the trained model uses H&E and nucleus information, without measured RNA from the target tissue. The separate core-to-whole-slide broadcasting experiment additionally uses measured RNA from selected tissue cores as anchors. Its zero-shot and anchored results must not be combined. Output is a predicted cell-by-gene matrix with tissue locations, not directly measured transcripts.

Manuscript settings
30 epochs; batch64; learning rate0.001; weight decay0.0001; first/second moments0.9/0.999; frozen UNI2h; five within-slide spatial folds with80/20 splits. The authors report an RTX A6000 with48GB VRAM. Their approximate training times are2h for one breast slide,72h for multi-slide breast and120h for pan-cancer. These are author-reported times, not local training measurements.

Sources
Manuscript Methods, especially p.21 and associated architecture descriptions. https://doi.org/10.64898/2026.06.09.731252
https://github.com/SydneyBioX/GHIST_plus
https://huggingface.co/MahmoodLab/UNI2-h

## Slide 4

What the screenshot proves
This is a real browser screenshot of a local page displaying the exact captured stdout of a command executed in the repository. It is not a screenshot of the native Terminal app. The browser rendering was used because the capture tool cannot access Terminal. The transcript, command, return code, capture timestamp and script SHA256 are stored in work/ghist_workflow_deck/evidence/transcripts.json. The screenshot does not show a fabricated training session.

Commands
From GHIST_plus, use the existing isolated analysis environment at ../work/ghist_plus_reproduction/.venv/bin/python. The presentation evidence command is:
python reproduction/presentation_evidence.py checks
It reads the pinned download/run records and reruns reproduction/verify_outputs.py before reporting the result. The displayed Figure2–Figure5 success messages refer to the recorded notebook execution on15October, not a new training run on22October.

Step1: verify source data
The 45 selected downloaded files have a retained manifest and source provenance. They total10.923GB. We use real targets and saved predictions.

Step2: align identifiers
A prediction is compared with the same measured cell ID and gene name. Aligning by row order alone could compare different cells and produce a misleading score. The evaluated cohort is saved rather than guessed from the full input dimensions.

Step3: execute notebooks
The official Figure2, Figure3, Figure4 and Figure5 notebooks completed successfully. Figure3 recomputes metrics from released predictions and measured RNA. Figure4 includes metric computation from released outputs. Figure2 and Figure5 primarily replot released summary tables. Successfully executing notebooks is narrower than reproducing the full figure experiments from original model training. The current Figure5 notebook is also not identical to the full manuscript Figure5 pancreas/GTEx analysis.

Step4: verify outputs
The verification checks execution counts/errors, cell/gene cohort integrity, imputation case counts and an independent SciPy PCC spot check. The independent check covers three SVG-ranked genes, rather than every possible gene/model metric. The paper-only report filters comparison columns to GHIST+, GHIST and SpatialEx. Undefined gene correlations remain undefined.

Useful follow-up commands
python reproduction/presentation_evidence.py data
python reproduction/presentation_evidence.py compare
python reproduction/presentation_evidence.py status
python reproduction/paper_only.py
All use the installed analysis environment. See reproduction/WORKFLOW.md for exact environment paths and setup.

Result locations
results/reproduction/Figure3/ and Figure4/ hold released notebook outputs. results/paper_only/ holds manuscript-scoped tables. reproduction/receipts/run_manifest.json and verification.json preserve the execution and checking evidence.

Sources
Local recorded receipts and scripts. Official notebooks: https://github.com/SydneyBioX/GHIST_plus

## Slide 5

How to read the table
Each row compares the manuscript's reported within-slide median PCC over its top20 spatially variable genes with the value calculated from the current released notebook outputs. Paper values are reference targets. Local values are measurements from released predictions against measured RNA; they are not substituted with the paper numbers.

Values
GHIST+: paper0.744, released0.7623199106216736.
GHIST: paper0.714, released0.7745567466037253.
SpatialEx: paper0.521, released0.637727915893977.
Rounded display values use three decimals for the paper and four for the local release. In this released within-slide comparison GHIST has the highest displayed value; do not hide or reinterpret that result as a GHIST+ improvement. The differing aggregation/cohort/configuration must be resolved before attributing any model ranking to the manuscript experiment.

What PCC means
For one gene, compare predicted expression with measured expression across the same cells. PCC measures whether they vary together. It ranges from−1 to1; higher is better in this use. A value of0.7623 does not mean76.23percent of cells or genes are correctly predicted. The median over20genes is a summary of20per-gene correlations. SVGs are genes whose expression changes spatially; top20is a selected subset and not all genes.

Why this is not an exact paper match yet
The manuscript uses five distinct spatial80/20folds. The current released configuration lists one validation division; its loader can reuse that division when iterated, so running fold IDs1–5 alone would not fix the discrepancy. Paper optimizer settings are learning rate0.001 and weight decay0.0001; current release settings are0.00115 and0.00001. Original cohort/fold definitions and some figure aggregation/labels need reconciliation. A later official code commit is not proof it generated an earlier manuscript figure.

Other paper metrics
The manuscript describes PCC and Spearman correlation (higher is better), RMSE (lower), SSIM for spatial pattern similarity (higher), CMD for co-expression structure difference (lower), TISSUE error score using its exact definition/aggregation, cell-type AUROC and Recall (higher), and Spatial IoU for spatial region overlap (higher). These are different endpoints, not one overall accuracy. This compact slide shows the directly auditable within-slide PCC comparison; the detailed registry and exported per-gene/per-class tables remain in the repo.

Critical RMSE caveat
Manuscript Eq15 uses log(1+10000×expression/library_size). The released imputation code uses a different raw-count/log scaling and method-specific transforms. Its RMSE must not be presented as matching the paper's normalized RMSE. The full measured panel and library size are needed to implement Eq15. Do not add a released RMSE next to a paper value on a different scale.

Screenshot provenance
Actual stdout of python reproduction/presentation_evidence.py compare, displayed in a browser and captured on22October. It reads results/paper_only/paper_vs_released.csv.

Sources
Manuscript results pp.5–6 and metric Methods pp.22–24. reproduction/paper_protocol/paper_values.csv includes page/figure references. results/paper_only/paper_vs_released.csv preserves exact values. https://doi.org/10.64898/2026.06.09.731252

## Slide 6

What remains
The current machine has an isolated analysis environment, not the full supported GPU training environment. No CUDA training device or approved UNI2-H weights were supplied to the preflight. It finds31missing aligned-data path fields. This is a count of configuration fields, not necessarily31distinct files. The release lists one validation division while the manuscript requires five distinct within-slide folds. This readiness check deliberately reports incomplete status.

1. Recover the exact prepared data and splits
For each slide obtain the registered H&E image; nucleus instance mask; nucleus coordinates and nucleus-to-Xenium matches; complete measured cell-by-gene panel and original library sizes; cell-type labels; panel-availability masks; and original fold manifests. Preserve missing genes separately from observed zeros. Confirm original spatial boundaries before running five folds. The existing evaluation bundle subset does not provide every required item. Do not synthesize any missing expression.

2. Prepare training environment
Follow the official repository's Linux/CUDA setup in a separate environment. Obtain UNI2-H through its official approved access flow and place genuine weights in the expected cache. The paper's authors used an RTX A6000 with48GB VRAM. The local notebook environment does not imply that model training is ready. No paid cloud job has been started.

3. Train the published experiments separately
Examples after real data paths and original folds are supplied:
export GHIST_DATA_ROOT="/path/to/prepared/real/data"
python train.py --config_file reproduction/paper_protocol/config_templates/within_slide.json --fold_id 1 --gpu_id 0
python train.py --config_file reproduction/paper_protocol/config_templates/breast_transfer.json --fold_id 1 --gpu_id 0
python train_gene_mask_imputer.py --config_file reproduction/paper_protocol/config_templates/imputation.json --fold_id 1 --gpu_id 0
python train.py --config_file reproduction/paper_protocol/config_templates/pancancer.json --fold_id 1 --gpu_id 0
These are conditional example commands, not commands that already completed. Templates restore the paper's explicitly stated optimizer settings but are not archived original paper configurations. Fold1is an example; all five distinct original spatial folds are required for the within-slide experiment. Never repeat one division and call it five-fold validation. The file config_pancancer.json is the correct public source for the pancreas template; the public filename config_all_cancers.json currently lists only breast sources.

4. Infer and evaluate
Use tools/inference.py with the actual training run directory, saved checkpoint, target slide and run-specific imputation cache. Do not guess paths before a run exists. Evaluate held-out cells/slides using the paper's fixed gene ranking, normalization, class definitions and aggregation. Save every eligible cell ID, gene name, seed, fold, checkpoint hash and per-gene result. Compare results with the manuscript without replacing local numbers.

Experiment checklist
Within-slide: BreastCancer2, five spatial folds.
Cross-breast: train five named breast slides, test BreastCancer2.
Imputation: mask real measured top10Giotto genes on each of five slides; preserve ground truth for evaluation.
VQ: reproduce the paper's full/no-VQ/random-prototype ablation, not a later unrelated checkpoint perturbation.
Pancreas: hold out PancreasADC3 and compare pancreas-only versus pan-cancer assistance.
Core-to-whole-slide: separate zero-shot predictions from measured-core anchors.
GTEx: repeat the77-sample normal breast application and downstream analyses; it is not a cell-RNA accuracy benchmark.

Definition of completion
Fresh model training and inference must succeed using a documented paper-matched protocol, followed by the manuscript's metric calculations. Differences should be explained with reproducible evidence. Merely executing figure notebooks is a useful evaluation milestone, but not the full reproduction.

Source
reproduction/WORKFLOW.md; reproduction/paper_protocol/protocol.json; reproduction/paper_protocol/config_templates/; reproduction/preflight.py; current status output captured22October.
https://github.com/SydneyBioX/GHIST_plus
https://huggingface.co/MahmoodLab/UNI2-h
