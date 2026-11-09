# GHIST+ implementation and experiment plan
Date: 30 October 2026
Status: proposed research plan; no new model training, inference or improvement metrics were produced while preparing this plan.

## Decision
Investigate learning gene expression under uncertain H&E nucleus-to-Xenium cell correspondence. Treat this as a gated feasibility project. Start with a data/geometry audit and a small controlled pilot; do not launch a full pan-cancer run first.

Research question: Can explicitly accounting for uncertain physical cell correspondence improve single-cell expression prediction, relative to improved registration and filtering/weighting, while retaining difficult cells?

This addresses an H&E nucleus being assigned to the wrong molecular cell. It does not by itself correct transcript spillover inside an otherwise correctly matched Xenium cell or recover molecular states that are invisible in histology.

Potential novelty is the training treatment of unresolved physical correspondence, not a new generic uncertainty head. GHIST+ already has a hurdle likelihood in the released implementation, and reference-based residual prediction. PHARAOH improves registration; VOICE filters nucleus matches. Neither the novelty nor an accuracy gain of our proposed method has been established.

## Current starting point
- Repository: GHIST_plus, branch reproduction/ghist-plus-local.
- Official implementation inspected: c5c895174c19903cc27d2661934b9b726bdb9306.
- Current local HEAD at inspection: 29c9c7a, with reproduction/presentation additions; model files match the official implementation.
- Existing work executed the authors' figure notebooks using real released predictions/measurements. Fresh model training and image inference have not been completed.
- Exact manuscript metrics have not been reproduced. See ../REPORT.md for differences in folds, normalization and later release outputs.
- Prior bundle pin: 37a937677c5a8ebe5c7b8fbd486fbc5bb50ce1de. The current web listing references 9d8ede79035c172e1ca14cb55692150158295fb8. Compare revisions and pin each experiment; do not silently mix old outputs and new files.
- The local BreastCancer2 measured matrix was a macOS cloud placeholder on the last read attempt. Validate complete file availability/checksums, preferably downloading chosen files directly to HPC scratch.

## Dataset order and acquisition
1. Official GHIST+ bundle: https://huggingface.co/datasets/SydneyBioX/GHIST-Plus-bundle
   Role: checkpoint/metric sanity checks, not a new external validation cohort. It includes evaluation files and predictions; it does not supply every raw aligned training input. Current release includes breast5k and ATERA evaluation, so those are not novel untouched benchmarks.

2. 10x breast: https://www.10xgenomics.com/datasets/ffpe-human-breast-with-pre-designed-panel-1-standard
   Role: image/coordinate smoke test using one selected sample. Official paired resources:
   H&E: https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FFPE_Human_Breast_IDC/Xenium_V1_FFPE_Human_Breast_IDC_he_image.ome.tif
   Xenium bundle: https://s3-us-west-2.amazonaws.com/10x.files/samples/xenium/1.0.2/Xenium_V1_FFPE_Human_Breast_IDC/Xenium_V1_FFPE_Human_Breast_IDC_outs.zip
   These are the vendor IDC sample, not an asserted mapping to paper alias BreastCancer2. Establish source-to-paper alias mapping separately if reproducing that exact experiment.
   Original paper preview source: https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast/
   Paper six breast sections must not be counted as six independent donors without source verification. Some 10x releases reuse samples/reprocessed experiments.

3. GSE311609: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311609
   Preferred independent research source, from a single-center study containing 27 donors/41 sections overall. Public H&E was added in June 2026. Do not assume all 41 sections have usable matched H&E.
   Main compact pilot: six NSCLC Prime5K donors L1,L3,L5,L6,L7,L8, samples GSM9509134–GSM9509139 respectively.
   Two H&E NDPI files contain these six donor regions:
   - GSE311609_hne_NSCLC_5k_L1_L1_NSCLC_5k_L3_L3_NSCLC_5k_L6_L6.ndpi
   - GSE311609_hne_NSCLC_5k_L5_L5_NSCLC_5k_L7_L7_NSCLC_5k_L8_L8.ndpi
   Sample outputs include RNA counts, morphology/DAPI images and cell/nucleus geometry. The series lists selected IHC images too. Ready-to-use transformation files have not been verified. Region extraction and registration are required until demonstrated otherwise.
   Second tissue: 17 breast donors/19 sections, B1 has three sections. Proposed initial breast audit cases: B10/GSM9509164, B12/GSM9509162, B14/GSM9509159, B8/GSM9509167. Raw files for the first three were checked in indexed sample listings; B8 per-file completeness remains to be checked.
   The raw archive is listed at approximately 179 GB. Acquire selected sample files and corresponding scanner images first; do not download the whole series for the first audit.

4. Optional high-panel breast stress test: https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-breast-cancer
   One donor; 5K plus custom panel; post-Xenium H&E and spatial measurements. Useful for realistic gene withholding if the project pivots. Already represented in GHIST+ release evaluation: describe as a shared public benchmark, not an independent novelty claim.

5. Optional external lung case: https://www.10xgenomics.com/datasets/xenium-human-lung-cancer-post-xenium-technote
   Two experiments from one donor, not two independent patients. Use only after core results and source overlap audit.

6. Optional later benchmarking resources:
   HEST: https://huggingface.co/datasets/MahmoodLab/hest
   CELLO: https://huggingface.co/datasets/gaozijun/cello_data
   HEST includes multiple assay types; select paired human Xenium entries. HEST access requires terms acceptance. CELLO's 52 pairs are not automatically 52 patients. Audit duplicate source studies/donors against 10x and GSE before calling either external.

Required per section:
H&E at known physical scale; DAPI/nuclear morphology; nucleus and cell boundaries/IDs; real cell-by-gene counts; transcripts where available; panel list; coordinate transform and direction; donor, block, section and scanner-image/scene mapping; assay version; raw and processed hashes. DAPI is not H&E. CZI/NDPI conversion must preserve physical scale, orientation and scene origin.

## Phase 0: freeze provenance and protocol (planning estimate: 1–2 working days)
- Preserve official files and existing manuscript reproduction.
- Create a separate research branch when implementation begins.
- Pin code, weights, bundle revision, dataset files, Python/PyTorch/CUDA and container/environment.
- Obtain approved UNI2-H access: https://huggingface.co/MahmoodLab/UNI2-h
- Reconstruct a COPY of released GHIST+ checkpoints using official instructions; preserve original checkpoint/hash.
- Create raw-source and processed-sample manifests. Verify real data are downloaded rather than cloud placeholders.
- Define output genes/order and scale. Intersect known genes for a checkpoint smoke test; a 5K dataset does not expand a checkpoint's trained vocabulary.
- Fix an H&E-only prediction task for this project. Target RNA/DAPI can support training or independent evaluation, but target RNA cannot enter prediction.
- Keep manuscript replication and new controlled research comparisons as separate protocols.

Deliverable: versions.json, dataset manifest, split proposal, data-ready report.

## Phase 1: verify pairing before any model training (3–5 working days)
- Start with one vendor breast case plus development-only GSE lung L1/L3 regions.
- Establish same-physical-section H&E/Xenium pairing. Adjacent sections are unsuitable for ground-truth one-nucleus identity.
- Verify pixel/µm conversion, origin, axes, pyramid level, affine direction, tissue region and nuclear ID mapping with overlays.
- Inspect approximately 300–500 nuclei per development section across 20–30 regions, stratified by tissue density and compartment. These are workload targets, not a formal power calculation.
- Review images/geometry without seeing model expression errors. Classes: clearly correct, clearly wrong, ambiguous, split/merge, no counterpart. Ask a second reviewer to independently assess a subset and adjudicate disagreements.
- Record region-stratified mismatch/ambiguity rates, duplicate claims, unmatched fraction, displacement and candidate recall. Ambiguity is not automatically a wrong match.
- Compare existing registration with a stronger registration procedure (PHARAOH or an appropriate validated pipeline).
- Assess whether a small local candidate set can recover adjudicated correct partners. Choose physical-radius thresholds from training/development image geometry rather than arbitrary pixel counts across scanners.
- Do not diagnose mismatch solely from an RNA prediction residual. Contamination and weak biological observability are separate explanations.

Gate A: continue only if material, potentially recoverable correspondence uncertainty remains after sensible preprocessing. If pairing is already reliable, or the dominant error is unresolved segmentation/transcript contamination, stop this branch.

## Phase 2: GPU preflight and ordinary baseline (1–3 working days plus queue)
- Run a few real image tiles through an official checkpoint; verify finite outputs, cell IDs, gene order, scales and coordinate overlays.
- Compare against released results only where image, mask, cell cohort and preprocessing actually match. Do not chase an exact manuscript number by tuning.
- Run one ordinary training step, save/resume, and a short validation pass before modifying the loss.
- Log peak GPU memory, host RAM, throughput and disk reads; estimate later job sizes from these measurements.
- Run identical inference under different batch partitions once. ECRM builds cross-tile neighbors within batches; document behavior and fix batch/graph policy.
- Start with few loader workers because current dataset code loads full images/masks into host RAM.
- Use test-independent preprocessing for the new benchmark. The released config's stain-reference and external-evaluation paths need auditing; fit stain references and make model-selection decisions using development data only.

Gate B: no full training until the baseline loads, trains/resumes and scores correctly with a fixed protocol and complete input provenance.

## Phase 3: implement the narrow contribution
Keep UNI2-H, ECRM, VQ, gene vocabulary and inference interface unchanged initially. Use training-only correspondence uncertainty.

Data:
- he_nucleus_id identifies the image prediction and graph node.
- xenium_cell_id identifies the molecular measurement.
- A sparse candidate table records plausible pairs and geometric evidence.
- Small candidate groups enforce one-to-one use of molecular cells. Initially handle groups of 2–5 nuclei; explicitly report excluded larger/unresolvable groups.
- Calibrate assignment probabilities from development-only reviewed geometry, then freeze them. Never use gene-expression similarity to choose a physical match.
- Fix unmatched status/cardinality through image QC initially. A learned null assignment must not be able to discard all hard examples.
- Keep graph nodes at H&E coordinates.

Loss:
Expose the existing hurdle likelihood without reducing across cells. For a group, sum likelihood over permitted one-to-one assignments using a log-sum-exp weighted by fixed geometric assignment probabilities. Normalize consistently by observed gene entries/group size. For a singleton certain match, recover the original loss.
Do not substitute an averaged neighboring RNA vector as ground truth. Compare simple weighted loss as a baseline.

Auxiliary supervision:
- In hurdle mode, morphology cell-type, pixel-map and composition supervision still remain.
- Use a common trusted subset for these RNA-derived labels across controlled comparison arms, or implement properly uncertainty-consistent labels.
- Compute cell-type expression priors from training-only RNA/reference data, with no test donor RNA.
- Preserve use_gt_ct=false during inference and avoid feeding evaluation RNA-derived types into ECRM.
- Keep MAGCI changes out of the primary experiment.

Proposed modules (not yet implemented):
- preprocessing/match_audit.py: overlays and reviewed matches.
- dataio/correspondence.py: candidate groups, probabilities and constraints.
- dataio/dataset_input.py and dataset_input_union.py: separate IDs, candidate targets, masks.
- model/hurdle_distribution.py: unreduced likelihood and assignment marginalization.
- train.py: candidate supervision and consistent auxiliary masking.
- dataio/references.py: training-only prior audit.
- evaluation/correspondence_benchmark.py: fixed cohort and stratified metrics.
- configs/research/: fixed split and variant configurations.
- jobs/: scheduler templates after HPC scheduler/partition details are known.

Required correctness checks using small real-data fixtures:
singleton equivalence; no duplicate molecular assignment; permutation invariance; masked NaNs excluded; fixed unmatched cases create no invented target; no cross-split candidate groups; overlapping tiles do not repeatedly count a cell's supervised loss; checkpoint resume.

## Phase 4: inexpensive controlled pilot (about one week)
Use development donors and modest spatial regions. A short 1–3 epoch/limited-update run is a feasibility screen, not a paper metric or adequate converged comparison.

Ablations:
A. Ordinary released GHIST+ workflow.
B. Better registration plus hard matching.
C. B plus geometry-based filtering.
D. B plus confidence-weighted loss.
E. Proposed assignment marginalization.
F. Hard-matching control with the SAME trusted auxiliary-label masks as E, isolating the new loss.

Keep genes, data exposure, preprocessing, initialization, encoder, batch/graph policy, optimizer and updates comparable. Filtering necessarily reduces supervised cells: report this and test coverage-matched controls.

Use fixed frozen features for a fast pilot only if cache parity is established. Cache pre-VQ nucleus features and raw tile features; do not freeze the output of trainable VQ/ECRM. Augmentation precedes encoding, so either precompute identical fixed views for all arms or label the no-augmentation pilot explicitly. Restore comparable augmentation in final online training.

Gate C:
- Benefit must repeat on independent development donors and exceed between-run noise.
- E must add value beyond B/C/D/F at comparable coverage.
- Do not advance solely because global PCC rises after dropping difficult cells.
- After estimating baseline variability, freeze a practical improvement target before final tests. An example planning target is ~5% lower normalized RMSE with non-decreasing PCC, not a promised or universal publication threshold.
- If better preprocessing/filtering performs equally well, prefer it operationally and stop the novel-loss claim.

## Phase 5: donor-held-out experiments
Proposed compact lung split, conditional on successful file/registration QC:
training L1/L3/L5; validation L6; untouched test L7/L8.
Keep all targeted/custom/5K sections from each donor in that donor's group. These two test donors provide a preliminary result, not broad population generalization.

For breast confirmation, after all paired regions are verified, reserve an approximately 10/3/4 donor train/validation/test partition from 17 donors, balanced using available clinical/acquisition metadata before examining model outcomes. All B1 sections stay together.
Scanner files can contain multiple donor regions. Retain scanner-file identifiers, crop correctly and assess whether shared scanner batches confound results; use an additional acquisition-group sensitivity split where feasible.

Train baseline and selected improvement on each tissue under identical budgets. Use three seeds for finalists if resource estimates permit; avoid a full model grid. If data/compute allow, add grouped donor cross-validation with the architecture and selection rules already frozen.
Only now consider expanding training to the paper's other tissues. Cross-tissue transfer is a separate question from the primary pairing intervention.

Do not use final test RNA for checkpoint selection, reference construction, gene selection, preprocessing tuning or pairing chosen to maximize expression agreement. Test H&E/DAPI can establish the fixed physical correspondence used by the evaluator; deployment prediction itself remains H&E plus nucleus information.

## Metrics and claims
Primary: normalized expression RMSE on a fixed, independently checked cell/gene cohort, with per-gene PCC alongside it.
Secondary: per-donor metrics; crowded/boundary strata; within-cell-type performance; low-abundance/rare-type groups where adequately supported; spatial map agreement; supervised/evaluated coverage; runtime and peak memory.
Report all supported genes plus a predeclared training-selected informative set. Document constant/near-zero-variance genes and undefined correlations. Do not select the best genes after seeing final test scores.
Use the same normalization for every arm; evaluate both absolute error and association. Define aggregation before fitting.
Bootstrap/report intervals by donor, not millions of cells as independent patients. With only a few donors, show each donor's result and acknowledge broad uncertainty.
Evaluate against real measured RNA but recognize that Xenium itself has noise/contamination. Optional IHC validates only its measured protein markers/cell identity, not thousands of RNA genes.
Independent adjudicated correspondences are the primary matching test. Larger automatically paired cohorts are secondary and must not silently replace trusted truth.
Measure effect of each intervention through ablation. Do not call gains from extra data or cleaned auxiliary labels a gain from the new likelihood.

## HPC resources and job sequence
The author implementation was tested on CUDA Linux/Python3.10/PyTorch2.6/CUDA12.4 with an RTX A6000-class GPU. The paper reports ~502 GB host RAM; that is an experimental configuration, not a measured minimum.
Request an available 48 GB-class GPU for the first faithful baseline if possible. Host-memory requirements must be profiled on selected images; the unmodified full-slide loader may require much more RAM than cropped/cached experiments.
Do not promise that a 24 GB GPU or 128 GB RAM is sufficient. Start small, record peak usage, then request measured peak plus headroom.
Lowering batch size and accumulating gradients does not restore identical cross-tile graph neighborhoods. Fix one comparable policy for all research arms.
Storage: inventory selected raw downloads first. Budget separately for expanded H&E, image pyramids, masks, frozen features, checkpoints and logs; compressed download size understates working space.
Download on an appropriate transfer/login service, compute through the HPC scheduler, and stage hot data on permitted local scratch. Keep credentials out of scripts, logs and Git.
Suggested job dependency order: manifest/preprocess -> match audit -> feature extraction/parity -> baseline smoke -> controlled pilot -> full finalists -> frozen evaluation -> report.
Do not assume SLURM/partition names before HPC details are available. No jobs have been submitted.

## Fallback if matching is not the bottleneck
Study support under real panel changes, as a separate hypothesis.
Use measured high-panel cells; hide authentic commercial panels or whole gene groups; compare random masking, structured masking and a specific support-conditioned correction. Compare with MAGCI and relevant VOICE behavior.
Distinguish H&E-only prediction from panel completion with observed RNA.
Mask hidden genes before input normalization/reference/label generation; hidden counts cannot contribute to observed-input library totals. Fixed target normalization for scoring is separate.
Gene withholding tests missingness, not probe chemistry. Use independently acquired GSE lung panels for an additional real-assay evaluation while keeping all versions of a donor in one split.
No fabricated expression data are required.
Generic masking/fusion alone is already crowded and is not an established novel contribution.

## Timeline and deliverables
Planning estimates exclude data transfer, queue delays and expert availability:
Week1: manifests, paired-image QC, baseline preflight.
Week2: blind matching audit, diagnostic baseline, go/no-go.
Weeks3–4: candidate loss + auxiliary masks, correctness checks, controlled pilots.
Weeks5–7: finalist training, independent donor and second-tissue evaluation.
Week8: failure analysis, uncertainty, literature recheck and professor/paper report.
Stop or pivot after Week2 if the mechanism is absent; do not spend the full HPC allocation to rescue a weak premise.

Required final artifacts: exact dataset/split manifests; reviewed-match subset; source/weight hashes; environments; reproducible preprocessing/training/evaluation entrypoints; per-donor results; ablation/coverage plots; limitations and failure cases; explicit distinction from manuscript reproduction.

## Sources
- GHIST+ paper: https://doi.org/10.64898/2026.06.09.731252
- Official code: https://github.com/SydneyBioX/GHIST_plus
- Official bundle: https://huggingface.co/datasets/SydneyBioX/GHIST-Plus-bundle
- UNI2-H: https://huggingface.co/MahmoodLab/UNI2-h
- GSE311609 cohort paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC13259927/
- GSE311609: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311609
- Author H&E deposition confirmation: https://github.com/bdsc-tds/xenium_analysis_pipeline/issues/111#issuecomment-4620206770
- Author data-processing notebooks: https://github.com/bdsc-tds/xenium_analysis_pipeline/tree/main/resources/notebooks
- PHARAOH: https://pmc.ncbi.nlm.nih.gov/articles/PMC13308372/
- VOICE, including nucleus matching Appendix D: https://arxiv.org/html/2608.08366v1
- CELLO: https://arxiv.org/abs/2609.36429
- 10x image alignment instructions: https://www.10xgenomics.com/support/software/xenium-explorer/latest/tutorials/xe-image-alignment
