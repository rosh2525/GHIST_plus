# GHIST+ implementation plan: auditing and improving RNA supervision

Date: 31 October 2026  
Status: implementation specification; experiments have not been run.  
Main execution target: one NVIDIA GPU with 48 GB VRAM. DGX access is optional.  
This revision supersedes the correspondence-first plan dated 30 October, preserved in archive/.

## 1. Project decision and faculty problem statement

**Research question:** Does uncertain cell-level RNA supervision limit GHIST+ at tumor–immune boundaries, and can ambiguity-aware supervision improve predictions beyond registration, auxiliary-label correction, and ordinary target purification?

We will first measure the problem. If the audit supports it, the primary candidate is **selective interval supervision using GHIST+'s native hurdle distribution**. On disputed targets, the model is allowed a fixed range between the observed measurement and a purification estimate instead of being forced to reproduce either value exactly. This is a proposed weak-supervision method, not an established decontamination algorithm or a verified novel contribution.

Faculty statement:
> We will reproduce a pinned GHIST+ release and test whether uncertain RNA assignment and purification affect its predictions on public human cancer tissue. Using donor-separated data and independently reviewed cell correspondence, we will compare ordinary training, better registration, label correction, target purification, and a training objective that represents disputed targets as ranges. The deployed predictor will use H&E and nucleus information. We will pursue the new method only if it adds value beyond simple cleaning and masking without removing rare cells or real expression variation.

Required project: one tissue, a defensible baseline, supervision audit, controlled comparisons, and a reproducible report. A second tissue, whole-cohort training, encoder adaptation, and DDP are extensions, not prerequisites.

No synthetic expression datasets will be used. Mathematical unit tests may use analytic scalar cases; model fitting and biological evaluation use real measurements.

## 2. Exact task and starting state

### Prediction task
Input at deployment: H&E, H&E-derived nucleus masks/coordinates, a fixed gene vocabulary, and reference statistics learned or selected using training data. Output: a cell-by-gene prediction at H&E nucleus locations.

DAPI, measured RNA, transcript coordinates, and IHC may support training-data preparation and independent evaluation. They must not supply target-slide expression or cell-type labels to the prediction forward pass. An evaluator can use fixed test H&E–DAPI correspondence; the predictor cannot use test RNA to choose matches or features.

Panel completion with observed target RNA is a separate task, excluded from the primary claim. A 5K assay does not give an existing checkpoint a 5K output vocabulary.

### Verified local baseline
- Official implementation inspected: c5c895174c19903cc27d2661934b9b726bdb9306.
- Local repository HEAD at inspection: 29c9c7a; branch reproduction/ghist-plus-local. Inspected model/training files match the official pin.
- Only analysis of real released predictions/measurements has been completed. No fresh image inference or model training has been completed.
- Exact manuscript reproduction is unresolved: see ../REPORT.md. For example, the within-slide top-20-SVG PCC is 0.744 in the manuscript and 0.7623 in our released-prediction analysis; this is not an improvement result.
- Prior bundle revision: 37a937677c5a8ebe5c7b8fbd486fbc5bb50ce1de; a later listing references 9d8ede79035c172e1ca14cb55692150158295fb8. Pin the chosen revision and audit changes.
- Some existing local data are cloud placeholders. Check byte availability and hashes; do not assume a file listing proves a readable dataset.
- Obtain authorized UNI2-h access and reconstruct copies of stripped checkpoints. Keep original hashes and third-party weights outside Git.

### Baseline behavior to preserve
Frozen UNI2-h image features, nucleus/tile pooling, VQ tissue context, ECRM graph context, reference-guided residual prediction, and the native hurdle head remain unchanged initially.

The released hurdle expression objective models log1p values with a separate occurrence probability and truncated-Gaussian positive magnitude. Alternative correlation, immune/invasive expression, and variance penalties are not all active in hurdle mode. Morphology cell-type/map/composition auxiliaries remain. Do not add a new correlation term and call it baseline preservation.

The multi-breast config uses predicted cell-type probabilities for graph inputs and same-predicted-type edge restriction; use_gt_ct is false. Preserve this in the controlled experiment. Any graph change is a separate ablation.

## 3. Data acquisition and initial split

Use dataset_sources.json as a discovery registry, not as proof of downloaded or fully usable inputs. All donor assignments remain provisional until image/source QC, then freeze before model outcomes.

| Priority | Source | Role |
|---|---|---|
| 1 | [Official bundle](https://huggingface.co/datasets/SydneyBioX/GHIST-Plus-bundle) | Released checkpoint and metric parity; selected files first |
| 2 | [10x breast standard](https://www.10xgenomics.com/datasets/ffpe-human-breast-with-pre-designed-panel-1-standard) | One raw-image/coordinate smoke test; not an assumed BreastCancer2 alias |
| 3 | [GSE311609](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311609) | Main research source; six Prime5K lung donors for compact development |
| 4 | GSE311609 breast donors | Optional confirmation or alternative main cohort if lung QC fails |
| 5 | [10x lung technical note](https://www.10xgenomics.com/datasets/xenium-human-lung-cancer-post-xenium-technote) | Optional external case; two experiments are not two donors |
| 6 | [HEST](https://huggingface.co/datasets/MahmoodLab/hest), [CELLO data](https://huggingface.co/datasets/gaozijun/cello_data) | Optional comparator data after source/donor overlap audit |

GSE311609 reports 27 donors and 41 sections overall. Do not assume all sections have usable matched H&E. The series includes multi-sample NDPI/CZI scans, so extraction and registration are real work. Approximately 179 GB is the listed raw archive, not a full working-storage budget.

Proposed uniform-panel lung split:

| Donor | GEO sample | Role |
|---|---|---|
| L1 | [GSM9509134](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9509134) | Training/development |
| L3 | [GSM9509135](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9509135) | Training/development |
| L5 | [GSM9509136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9509136) | Training |
| L6 | [GSM9509137](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9509137) | Validation |
| L7 | [GSM9509138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9509138) | Reserved test |
| L8 | [GSM9509139](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9509139) | Reserved test |

All other panels, serial sections and replicates from each donor stay in that donor's split. Two test donors support preliminary case-level evidence, not strong population-level significance. More training seeds do not increase the biological sample size.

Two scanner images associated with these samples are:
- GSE311609_hne_NSCLC_5k_L1_L1_NSCLC_5k_L3_L3_NSCLC_5k_L6_L6.ndpi
- GSE311609_hne_NSCLC_5k_L5_L5_NSCLC_5k_L7_L7_NSCLC_5k_L8_L8.ndpi

Verify exact crop/scene-to-GSM mapping, same physical section, transformation availability, and IHC marker/image linkage. Record scanner-image IDs because shared scans can confound acquisition effects. Add a scanner-group sensitivity analysis if sample counts permit.

For breast expansion, the study has 17 donors/19 sections. A proposed 10/3/4 donor split must be assigned from verified metadata before outcomes, with all B1 replicates grouped. Do not automatically use all 41 sections or call serial sections paired single-cell ground truth.

Prime breast is already represented in the GHIST+ released evaluation. It can be a shared benchmark, not an untouched new test dataset. No ancestry claim is supported merely by the study's hospital location.

### Required files and provenance per section
H&E; DAPI/nuclear image; nucleus and cell geometry; real counts; gene panel; transcript coordinates/quality when available; cell IDs; donor/block/section IDs; assay/software version; physical pixel size; transform and its direction; scanner/scene/crop metadata; download URL, file bytes, SHA-256, and access terms.

Use selected GEO sample files and scanner images before downloading the whole archive. Record exact resolved file URLs in a download manifest; do not invent missing direct links. Convert images without losing pixel size, orientation, origins or pyramid information.

Raw counts, nuclear counts, and SPLIT estimates remain separate arrays. Nuclear counts measure a different compartment. SPLIT estimates are not clean ground truth. Restrict nuclear transcripts using verified nuclear assignments/masks, and preserve the transcript quality threshold across reconstructions.

## 4. Weeks 1–2: mandatory audit and baseline preflight

### Work package P0 — freeze versions and data policy
1. Create a research branch when code implementation begins; keep the manuscript reproduction intact.
2. Pin code, checkpoints, encoder, purification and registration versions, and environment.
3. Build download/sample manifests and the donor split registry.
4. Inspect actual input sizes and availability. Hydrate cloud placeholders or download directly to HPC scratch.
5. Define H&E-only inference and the exact prediction scale.
6. Freeze a pilot gene list from training data. For checkpoint parity use its original ordered vocabulary; for new lung training, use the common measured panel, subject to memory profiling. Any reduced-gene pilot is labelled as such and is identical across arms.
7. Record raw data terms and approved encoder access. Store no credentials or redistributed encoder weights in Git.

### Work package P1 — image correspondence audit
Start with development L1/L3 regions and one vendor smoke-test case.
- Verify same-physical-section pairing, units, axes, flips, pyramid level, crop origin and transform direction.
- Compare provider registration with PHARAOH or another verified registration procedure.
- Review approximately 300–500 nuclei across 20–30 regions per development section; include tumor core, immune aggregates, boundaries and sparse stroma. These are workload targets, not a statistical power claim.
- Score correct, wrong, ambiguous, split/merge and unmatched correspondence without viewing GHIST+ residuals.
- Have a second reviewer assess a subset and adjudicate. Save contours, coordinates, evidence and reviewer agreement.
- Record displacement, duplicate assignments, missing cells, geometric ambiguity and candidate recall by region.

Do not assume fast registration eliminates all uncertainty. Do not relabel nuclei by whichever RNA profile fits the model best.

### Work package P2 — contamination-associated supervision audit
On geometrically reliable matches:
- Construct original observed, nuclear-compartment and pinned SPLIT-derived targets separately.
- Use training-donor or independently sourced reference profiles; audit donor overlap.
- Record correction magnitude, reference sensitivity and transcript locations for a prespecified marker subset.
- Seek orthogonal support from available registered IHC or independently reviewed lineage. Verify exact marker/sample availability first.
- Compare boundary and interior cells with matching or stratification for lineage, RNA abundance, nuclear size and local density. These are observational diagnostics, not proof of a single causal error source.
- Compare baseline predictions with disputed signal and examine both false-positive and positive-marker sensitivity. A model that predicts zero for every marker must not pass.
- Treat cytoplasmic localization as biological information, not automatically contamination. Nuclear restriction can lose real signal.
- Residual error in apparently clean cells does not establish that the biology is inherently invisible in H&E.

### Work package P3 — baseline execution
- Reconstruct and load one author checkpoint with approved encoder weights.
- Run real tiles and then an identified evaluation cohort. Check IDs, scales, genes, finite outputs and coordinate overlays.
- Compare outputs to the released predictions only when checkpoint, cohort, images, masks and preprocessing match. Diagnose discrepancies; do not tune to force a manuscript PCC.
- Train a small real-data baseline, save/resume, and run validation.
- Profile warm-up plus 50 batches at the actual gene dimension: GPU peak memory, host memory including workers, loader wait, step time, validation time and checkpoint size.
- Test a no-target inference route: no target RNA/cell-type labels in model inputs or cell-selection side effects.
- The impute cache contains per-source expression CSVs and masks, not just model weights. The inspected inference forward pass supplies batch_expr=None and batch_ct=None, but loader filtering, reference construction and metrics still need a dataflow audit. Cache presence alone is not proof of leakage.

### Gate G1 — written decision by end of week 2
Proceed with interval-supervision experiments only if:
- enough usable same-section pairs and independent lineage evidence exist;
- disputed targets are reproducibly identifiable beyond geometry alone;
- baseline predictions show a relevant failure on development data;
- the intervention can retain adequate trustworthy supervision for rare and common lineages;
- the baseline executes correctly and actual resource costs fit the allocation.

Freeze audit thresholds and acceptable marker-sensitivity/coverage margins from development evidence before final tests. There is no universal PCC cutoff proving publishability.

If geometry dominates, use the correspondence branch in Section 11. If neither hypothesis is supported, deliver the audit and use the panel-transfer fallback; do not force a new loss.

## 5. Conditional primary method: selective interval supervision

### Intended interpretation
Test whether a model should be less certain about disputed supervision. The interval is a **heuristic acceptance range**, not a confidence interval, a physical lower/upper bound, or a guarantee that true intracellular expression is contained within it. Because it includes the original observation, it does not automatically remove contamination.

The key comparison is against ordinary GHIST+ trained on exactly the same purification estimates and auxiliary labels. No new neighbor mixing or count-conservation claim is made.

### Offline metadata
For each matched training cell i and measured gene g:
- x_obs: original nonnegative counts.
- x_alt: alternative nonnegative purification estimate, with method/version/reference provenance.
- m_ig: whether that gene was actually measured.
- a_ig: fixed binary flag for disputed supervision.
- w_i: fixed image-derived correspondence confidence.
- trusted auxiliary label/mask and evidence provenance.

Define y_obs=log(1+x_obs), y_alt=log(1+x_alt).
For flagged entries, l=min(y_obs,y_alt), u=max(y_obs,y_alt).

Build flags from a frozen development-only rule combining target disagreement with independent image/transcript/lineage evidence and reference sensitivity. Start with auditable lineage-informative genes; do not expand to every gene based on poor model prediction. Preserve quantitative thresholds and reasons in metadata. Nuclear counts support diagnostics, not universal interval endpoints.

The flag and interval generators never consume model residuals or test-donor RNA to tune rules. If a correction increases expression, retain endpoint ordering and provenance; do not assume purification only removes counts.

### Distribution and objective
Use the existing hurdle parameters:
- q_ig=sigmoid(occurrence_logit_ig)
- mu_ig: signed magnitude location in unscaled log1p units
- sigma_g: existing positive-magnitude scale

Z has mass 1-q at zero; conditional on Z>0 it follows a Normal(mu,sigma^2) truncated to positive values. For u>l>=0:

P(Z in [l,u]) =
  (1-q) * 1(l=0)
  + q * {Phi((u-mu)/sigma)-Phi((l-mu)/sigma)} / Phi(mu/sigma).

Compute probability in stable log space using tail-aware log-CDF differences. Handle zero mass explicitly. Equal endpoints use the native point likelihood. Exclude unobserved/padded entries before arithmetic; reject malformed or nonfinite endpoints on eligible supervised entries.

For the exact controlled comparison, D uses the native point NLL on y_alt. E uses that same y_alt point NLL on unflagged entries and interval probability on flagged entries:

L_E =
 sum_ig w_i m_ig [(1-a_ig) * NLL_native(y_alt_ig)
                  + a_ig * {-log P_theta(Z_ig in [l_ig,u_ig])}]
 / sum_ig w_i m_ig.

For M, normalize the expression loss by its surviving effective weights sum(w*m*(1-a)); for W, normalize by the declared downweighted effective weights. Log these denominators and effective supervision coverage. If an expression denominator is zero, return a differentiable zero expression loss and retain any eligible auxiliary losses.

All other losses, labels, reference priors, model initialization and graph policies are identical between D and E. To recover D, set every a_ig=0. To recover the original supervised objective, additionally restore observed targets, original labels/weights/masks and all original settings.

This is a weak-constraint training objective, not proof of an interval-censoring data-generation process. Do not compare its raw training-loss magnitude with point-NLL runs to pick checkpoints; use the same validation outcomes.

The current loader applies expr_scale * log1p; one inspected config uses expr_scale=2.0. The native likelihood wrapper divides signed magnitude and targets by this scale. Convert interval endpoints exactly once and keep sigma in the likelihood's units. Never feed cohort-gated exported predictions into this objective.

Monitor whether wider sigma, excessive zeros or loss of within-type variation explains apparent gains. Trusted point-supervised entries anchor the model. Keep sigma constraints identical across arms; any new regularizer needs an ablation.

### Auxiliary labels, unmatched cells and priors
- Use a common trusted auxiliary-label subset for C/D/E/M/W. Keep disputed labels masked rather than accepting every SPLIT reassignment as truth.
- Mask cell-type, pixel-map, expression-to-type and composition targets consistently wherever active. Recompute tile composition over eligible labels and record its denominator.
- Retain unmatched H&E nuclei as graph context, with no direct RNA/label supervision. Record them explicitly; do not fill missing targets with supervised zeros.
- Context-only nodes can receive indirect gradients through graph messages; this is expected.
- Fix priors across the primary D/E/M comparison. Use training-only corrected or original profiles according to one declared common policy. Separate prior-cleaning changes from loss changes.
- Use a uniform measured-panel primary experiment. Keep gene-mask imputation auxiliary reconstruction off for all controlled arms initially, with that shared change documented. Later restore it for a compatibility check, routing all target-consuming reconstruction branches through the same supervision policy. The exact author configuration is retained separately as A0.
- Do not change ECRM edge features or feed target labels to inference.

## 6. Controlled experiments

The full variant definitions are in EXPERIMENT_MATRIX.csv. Every controlled arm uses the same input cell/image cohort, genes, geometric graph, augmentation policy, optimizer budget and fixed evaluation cohort. Effective RNA-supervision exposure intentionally differs for masking/weighting arms and must be reported.

| ID | Supervision | Purpose |
|---|---|---|
| A0 | Pinned author workflow | Release reproducibility, separate protocol |
| B | Improved registration; observed RNA; original auxiliary labels | Controlled research baseline |
| C | B with trusted/corrected auxiliary-label policy | Isolate auxiliary-label effect |
| D | C with purification estimates; ordinary hurdle | Mandatory target-cleaning baseline |
| E | D with selective interval objective | Proposed method |
| M | D but mask flagged gene targets | Mandatory simple alternative |
| W | D but downweight flagged entries | Optional weighting baseline |
| N | Nuclear-compartment targets | Separate sensitivity analysis, not whole-cell truth |

B–E/M/W use the same controlled uniform-panel configuration, including the shared imputation-auxiliary policy. If graph nodes or training cells differ between A0 and B, describe that whole preprocessing change instead of attributing it solely to registration.

Primary mechanistic comparisons: E versus D and E versus M; C versus B isolates auxiliary changes. Use a simple training-reference mean/cell-type predictor as a diagnostic for excessive reference shrinkage; distinguish oracle-type and image-predicted-type versions.

Short pilots screen behavior; they do not establish converged superiority. Adequately train B/D/E/M before final claims. W is required if it is competitive during screening. N does not compete for the same biological estimand.

Gate G2 after pilots:
- no further model branch if E is indistinguishable from D or M within development variability;
- no advancement if specificity improves only by losing positive-marker sensitivity, rare cells or prediction variance;
- no advancement on SPLIT-agreement alone;
- select one configuration, a fixed training budget and final evaluation rules before unsealing tests.

## 7. Code implementation work packages

Paths below are repository-relative planned changes; proposed files do not yet exist.

| File / component | Work | Acceptance condition |
|---|---|---|
| tools/audit_supervision.py — new | Registration overlays, reviewed pairing, transcript/marker diagnostics | Reproducible review table; no model-driven matching |
| tools/prepare_supervision.py — new | Versioned observed/alternative targets, flags, bounds, geometry weights and masks | Same source IDs/scales; reference and split provenance |
| dataio/supervision.py — new | Typed metadata reader, ID joins, bounds/mask validation | Fails on duplicate IDs, gene-order mismatch or invalid bounds |
| dataio/dataset_input.py::get_input_data / DataProcessing | Separate all H&E context nodes from supervised nodes; preserve transforms | Unmatched nuclei retained without invented targets |
| dataio/dataset_input_union.py::DataProcessingUnion | Return cell-gene metadata alongside measured-gene masks | Collation respects cell/gene order and padding |
| model/hurdle_distribution.py | Unreduced native NLL plus stable hurdle interval log-probability | Point equivalence and valid gradients |
| train.py | Central variant routing, reductions, auxiliary masks, resume/config logging | B/D/E differ only as declared |
| dataio/references.py | Trace build_train_region_expression_fallbacks and build_train_region_avgexp_df_by_slide | Test donors cannot enter priors |
| tools/inference.py | Audit target-free path; only fix hidden evaluation dependencies if found | Same H&E predictions when evaluation RNA is unavailable |
| tools/evaluate_supervision.py — new | Fixed-cell, per-donor/lineage metrics and coverage | Identical evaluation cohorts; undefined correlations retained |
| configs/research/ — new | Explicit B–N configs derived from pinned baseline | Manifest and hypothesis IDs in every run |
| jobs/ — new later | Preprocess, pilot, train and evaluate job wrappers | Scheduler-independent commands until allocation details are known |

Store cell metadata in Parquet; store dense or sparse count/target arrays in a chunked format after profiling, with a separate ordered genes file. Avoid a dense cell-by-gene-by-method JSON or per-cell file explosion.

Required sample schema:
study_id, donor_id, block_id, section_id, assay_run_id, scanner_image_id, scene_id,
crop_origin, physical_pixel_size, transform_hash, he_nucleus_id, xenium_cell_id,
match_status, geometry_weight, auxiliary_label, auxiliary_valid, reference_version.

Separate measured_gene_mask, direct_supervision_mask, ambiguity_mask and padding masks. A missing measured gene is not a biological zero. Context nodes and evaluation coverage cannot be selected by the method's prediction quality.

### Correctness checks
Use a small consent-compatible real-data fixture, plus analytic distribution tests:
- Original native point loss/gradients unchanged after refactoring.
- Zero flags reproduce D exactly; complete baseline restoration reproduces B.
- Interval probability lies in [0,1], grows with interval inclusion, handles zero mass and extreme tails.
- Zero-width intervals dispatch to the defined point rule; tiny intervals are numerically stable.
- Masked NaN/Inf values do not contaminate arithmetic.
- expr_scale conversion occurs once.
- Padding/context-only cells have zero direct target loss; indirect graph gradients are allowed.
- Augmentations transform image, mask, IDs and coordinate metadata consistently.
- No duplicate supervised cell from overlapping tiles in a cohort reduction.
- Target-free inference and fixed graph policy give stable cell IDs/outputs.
- Save/resume restores optimizer, random state, epoch, sampler and method settings.
- No held-out donor enters reference, flag calibration, normalization, gene ranking or checkpoint selection.

Do not require every noisy training loss to decrease at every step. Check finite gradients, reproducibility under controlled settings and appropriate aggregate learning.

## 8. Evaluation and anti-leakage rules

### Outcomes
Before final evaluation, freeze one primary independently supported biological outcome and one reconstruction guardrail based on actual assay coverage. If there is no adequate independent outcome, the claim is limited to robustness across alternative targets.

For available markers, compare false-positive signal in independently supported negative lineages while maintaining sensitivity in positive lineages. Freeze thresholds from development data. Report marker-level results, not only a composite score that rewards predicting zeros.

Also report:
- per-gene PCC and Spearman, plus log-scale RMSE/MAE;
- raw measured, nuclear and purified-reference agreement in separate tables;
- within-lineage correlation/error where nonzero variance supports it;
- rare/low-RNA lineages, boundary versus interior and density strata;
- per-gene occurrence rates, predicted/observed variance and reference-mean baseline;
- cell coverage and supervised-entry coverage, including exclusions;
- runtime, host memory, GPU memory and I/O.

Preserve the native log1p scale as one evaluation. For library-normalized metrics, first recover count-scale predictions according to their recorded export scale: divide scaled-log1p values by expr_scale where applicable, then apply expm1; do not repeat this transform for exports already in counts. Never sum log1p values as though they were counts. Normalize observed and predicted count profiles using their own declared totals on the same vocabulary; never scale predictions by a test cell's true RNA total. Any train-fitted scaling is fixed. Describe how zero-total profiles are handled.

Select informative genes using training data only. Report all supported genes and constant-gene exclusions. Do not rank final test genes to improve the headline.

Use donor summaries first, seeds second. With L7/L8 only, show both donor effects and uncertainty limits; do not manufacture statistical power from cells or seeds. More donors or a confirmation cohort are needed for broader claims. Shared scanner images require an acquisition sensitivity analysis where feasible.

### Independent truth and target leakage
Fix test correspondence before comparing model outputs. Do not let each model pick its own matches. Keep a reviewed cohort and a wider automatic cohort separate.

IHC supports selected protein/lineage outcomes, not transcriptome-wide RNA truth. Reference-purified RNA may be evaluated as a secondary estimate but cannot be the sole arbiter for a model trained on it. Compare reference choices on development data and report sensitivity.

Keep test data out of model selection, stain-reference fitting, reference means, label-model fitting, flag thresholds and gene selection. Purely image-based test QC is allowed under predeclared rules. For the no-target forward test, also ensure removing RNA does not change the H&E node cohort or graph.

## 9. Single 48 GB GPU execution plan

Authors reported a 48 GB RTX A6000, 16 CPU cores/32 threads and 502 GB host RAM. Their reported training times (within 2 h single-slide, 72 h multi-breast, 120 h pan-cancer) describe their experiments, not guaranteed runtimes for GSE311609 or minimum requirements.

Initial request: one 48 GB NVIDIA GPU; preferably 128–256 GB host RAM for a cropped pilot; approximately 0.5 TB fast scratch if available. These are planning requests, not proven sufficiency. Full-slide loading can require more; request measured peak plus headroom.

Execution order:
1. Transfer and checksum selected files.
2. CPU image/metadata preparation; registration and audit.
3. GPU segmentation only where required.
4. Baseline checkpoint inference and 50-batch profile.
5. Limited-update pilot variants, with an initial spending cap such as 12 GPU-hours.
6. Converged finalists after G2.
7. Frozen evaluation and exports.

The pilot cap is a budget, not a claim that full training fits in that time. Allocate equal updates/regions for screened arms and carry promising arms to convergence.

Memory/I/O:
- Begin with few loader workers; account for images/masks replicated or shared across workers.
- Stream large transcript tables and process scanner regions rather than holding all slides in RAM.
- Stage hot data on local scratch. Inventory compressed, expanded, temporary, mask, feature and checkpoint bytes separately.
- Profile batch64; reduce only when necessary. Gradient accumulation does not restore cross-tile graph neighbors lost with smaller batches.
- Cache only frozen pre-VQ features for a declared fixed-view pilot. A single cache removes on-the-fly image augmentation. Do not cache trainable VQ/ECRM outputs.
- AMP/BF16 is a measured option; verify numerical stability/parity. Do not promise a speedup.

If a full run takes t GPU-hours, four core finalists B/D/E/M with three seeds cost 12t, plus pilots and evaluation. At an illustrative t=72 h that is 864 GPU-hours, or 36 continuous GPU-days. A smaller budget can first converge all four once and add two seeds to D/E (8t total, 576 hours or 24 continuous days at t=72 h); claims must acknowledge fewer baseline repeats.

At G2, compare the measured total with actual allocated GPU-hours and availability. Weeks 5–7 provide only 504 hours even with uninterrupted single-GPU access; the illustrative 864-hour and 576-hour budgets do not fit that window. If needed, extend training through weeks 9–10 and move frozen evaluation/reporting to weeks 11–12, reduce all arms' cohort/gene scope equally, or use independent DGX jobs. Reserve measured validation/checkpoint overhead and queue contingency. If the allocation still cannot fit, reduce the promised project scope before committing.

## 10. Optional DGX execution plan

First use task parallelism: separate single-GPU jobs for seeds and variants. Give each job unique experiment, cache and output paths. Bound concurrency by shared host RAM, CPU and scratch throughput. A DGX label alone does not specify the assigned GPU model or memory.

DDP is not required for the first project. If one run becomes the bottleneck, add it on a separate branch:
- process group/local device binding and distributed sampling;
- slide-homogeneous/group-aware batches and explicit global versus per-rank graph policy;
- differentiable feature exchange if preserving a graph across ranks, or documented changed neighborhoods;
- correct observed-entry/global weighted loss reductions;
- synchronization of codebook statistics only if the selected implementation uses mutable/EMA updates;
- synchronized sampler state, random seeds, unused branches and resume behavior;
- rank-safe caches, checkpoint writes and unique-cell evaluation aggregation;
- one-GPU versus multi-GPU parity and measured scaling.

Splitting batch64 across GPUs does not preserve the original cross-tile graph automatically. Independent jobs also do not scale linearly if storage or CPU is saturated. No guaranteed DGX speedup is assumed.

## 11. Fallback decisions

### F1: correspondence dominates
Compare strong registration + hard matching against geometry-confidence weighting. If local ambiguity remains, test small-group one-to-one assignment marginalization with fixed geometry-only priors:

L_group = -log sum_A q_geo(A) exp(-sum_(i,j in A),g m_jg * NLL_native(theta_i,y_jg)).

Combine group log-likelihoods before normalization by the common observed-entry count. Use feasible fixed-cardinality assignments; handle unmatched/split/merge cases explicitly and keep auxiliary labels trusted. Small groups only at first. Weighted candidate NLL is a distinct baseline, not the same objective; it can still produce averaging behavior.

No expression-residual-driven candidate selection. Test matches remain independent. If hard matching or weights solve the problem, stop the new matching-loss claim.

### F2: no convincing supervision bottleneck
Use the existing data for a separately defined panel-transfer study. Compare random masking, actual panel masks and structured gene-set withholding. Distinguish genes missing from one panel from genes never supervised anywhere.

Real cross-panel sections from the same donor are not the same cells. Validate cell-level withholding on measured high-panel cells; use separate acquisitions for donor/lineage/distribution-level transfer questions. Preserve donor grouping, mask before input/reference normalization and do not infer probe chemistry from gene masking.

MAGCI/VOICE already cover related mechanisms. Structured masking alone is not a verified novel contribution.

### F3: corrected targets help but E does not
Use the simpler successful preprocessing method. Produce a reproducible supervision-quality benchmark and failure analysis. Publication value depends on evidence and scope; do not rename a null loss result as a new architecture.

## 12. Schedule and deliverables

| Week | Deliverable / decision |
|---|---|
| 1 | Version/data manifests, access, real sample conversion, baseline preflight |
| 2 | Blind correspondence and supervision audit; G1 decision; frozen development protocol |
| 3 | Metadata loader, unreduced native loss, interval objective and correctness tests |
| 4 | Controlled pilots, simple baselines and G2; measured final budget |
| 5–7, or 5–10 | Converged finalists and seeds; G2 selects the window from measured runtime and allocated hours |
| 8, or 11 | Frozen donor-level evaluation, coverage and marker/state failure analysis |
| 9–10, or 12 | Reference/threshold/acquisition sensitivity, nearest-prior-art recheck, writing |
| 11–12 on the faster path only | Contingency or additional donors/second tissue if the main study is complete |

Eight weeks is possible for a narrow completed study; twelve provides contingency and stronger validation. Queue delays, data registration and expert review can move dates.

Required outputs:
- exact source and donor/section/crop manifests;
- pinned environment and source/weight hashes;
- independently reviewed correspondence/marker audit;
- executable preprocessing/train/evaluate entrypoints once implemented;
- complete run configs/checkpoints/logs and per-donor results;
- B/D/E/M comparisons, coverage and failure plots;
- explicit limits on purified-reference truth and publication novelty.

No result is called improved biological accuracy based solely on a lower training loss, agreement with SPLIT, higher PCC after exclusions, or additional compute.

## 13. Source and evidence index

This plan's method, gates, resource requests and schedule are proposals. Scientific facts are grounded in the sources below; dataset readiness beyond inspected listings remains to be established.

- [GHIST+ manuscript](https://doi.org/10.64898/2026.06.09.731252)
- [Pinned training implementation](https://github.com/SydneyBioX/GHIST_plus/blob/c5c895174c19903cc27d2661934b9b726bdb9306/train.py)
- [Pinned hurdle implementation](https://github.com/SydneyBioX/GHIST_plus/blob/c5c895174c19903cc27d2661934b9b726bdb9306/model/hurdle_distribution.py)
- [Official bundle](https://huggingface.co/datasets/SydneyBioX/GHIST-Plus-bundle)
- [UNI2-h model and access terms](https://huggingface.co/MahmoodLab/UNI2-h)
- [SPLIT/cohort paper](https://www.nature.com/articles/s41592-026-03089-8), published online 30 April 2026
- [SPLIT implementation](https://github.com/bdsc-tds/SPLIT): pin chosen version; current repository describes newer annotation-agnostic and residual-cleaning options, not just legacy RCTD
- [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311609)
- [Author image-deposition discussion](https://github.com/bdsc-tds/xenium_analysis_pipeline/issues/111#issuecomment-4620206770)
- [PHARAOH paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13308372/) and [implementation](https://github.com/yaosicong1999/PHARAOH)
- [VOICE](https://arxiv.org/html/2608.08366v1), [CELLO](https://arxiv.org/abs/2609.36429), [SPATIA](https://arxiv.org/abs/2507.04704)
- [Proseg](https://github.com/dcjones/proseg), [ResolVI](https://pmc.ncbi.nlm.nih.gov/articles/PMC13645626/)

Before a novelty claim, record an evidence table comparing interval/weak supervision, reference correction, robust losses and nearest H&E methods. A targeted search not finding this exact combination does not prove novelty. No later-paper absence or publication acceptance is guaranteed.
