# Slide content: GHIST → GHIST+ → our improvements (status 2026-11-11, 17:40 IST)

Approaches, in the order we built them:
- **Approach 1:** residual adapters (v1–v3)
- **Approach 2:** soft output + core correction + uncertainty
- **Approach 3:** BioHAN neural network + biology-driven extensions (on the IITJ HPC)

Every comparison uses identical cells, a protocol frozen before the run, and block-bootstrap 95% CIs. **[PENDING]** marks results still running.

---

## Slide 1: About GHIST
**GHIST: single-cell gene expression from routine H&E**
- **Paper:** Fu et al., *Nature Methods*, 2025.
- **Code:** https://github.com/SydneyBioX/GHIST
- **Input:** an ordinary H&E slide. **Output:** expression of every gene in a targeted panel (about 280 genes on BreastCancer2), in every cell.
- **Training data:** 10x Xenium spatial transcriptomics measured on the same section, so each nucleus is paired with its measured RNA.
- **Why it matters:** spatial transcriptomics is expensive and rare; H&E is in every pathology archive. It enables virtual molecular maps of tumours at scale.
- **Biological idea:** a cell's expression depends on what it looks like, what type it is, and who its neighbours are.

## Slide 2: GHIST architecture (checked against the released code)
1. **Preprocessing** (scripts released):
   - nuclei segmented on H&E with HoverNet
   - H&E aligned to Xenium
   - nuclei matched to Xenium cells
   - cell types labelled from measured RNA (scClassify)
2. **Image encoder:** a UNet3+-style CNN (full-scale skip connections, 64–1024 channels), trained from scratch on 256×256 H&E patches.
3. **Per-cell features:** encoder features averaged inside each nucleus mask, plus patch-level context features.
4. **Expression head:** each cell's expression is a learned mix of **cell-type reference expression profiles**.
   - The mix is refined by cross-attention with the patch's **estimated neighbourhood cell-type composition**.
   - There are separate branches for **immune** and **invasive tumour** cells.
5. **Auxiliary heads, trained jointly:**
   - cell-type classification (from image features, and from the predicted genes)
   - neighbourhood composition
6. **Output:** a single-cell expression map of the whole slide.
- **Caveat:** cell-type labels are themselves predictions from RNA, not ground truth.

## Slide 3: Successor GHIST+: what it fixed, and metrics
**GHIST+** (Zhang, Yu et al., bioRxiv, June 2026)
- **Code:** https://github.com/SydneyBioX/GHIST_plus
- **Our fork with all our code:** https://github.com/rosh2525/GHIST_plus (branch `ghist-plus-soft-adapt`)

| GHIST limitation | GHIST+ change |
|---|---|
| CNN trained from scratch on limited data | **UNI2-h** pathology foundation model (ViT-H/14, frozen), 224×224 patches |
| Neighbours used only through a composition estimate | **ECRM neighbour graph**: each cell reads about 12 neighbours' features, also across patch borders |
| No notion of recurring tissue states | **VQ tissue prototypes** (64 codes), also used to choose TMA cores |
| Locked to one gene panel | **MAGCI imputer** fills genes missing from a slide's panel |
| Breast-only training | **Pan-cancer training** (breast, lung, colorectal, pancreas) |
| Zero-inflated RNA treated as plain regression | **Hurdle output**: P(expressed) × amount |
| No way to adapt to a new slide | **Ridge correction (RRC)** from a few measured cores (`tools/one_shot_rrc.py`) |

**Metrics (median per-gene PCC)**

| Setting | GHIST+ | GHIST | SpatialEx |
|---|---|---|---|
| Within-slide, top-20 SVGs: paper | 0.744 | 0.714 | 0.521 |
| Within-slide, top-20 SVGs: our rescoring of the authors' released predictions (12,719 cells) | 0.762 | **0.775** | 0.638 |
| Within-slide, all 280 genes (our rescoring) | 0.312 | 0.308 | 0.153 |
| Unseen slide, top-20 SVGs (paper) | 0.653 | – | – |
| Masked-gene imputation (paper) | 0.743–0.747 | – | – |
| External slides, all genes (our measurement) | Breast5K 0.070, ATERA 0.109 | – | – |

- The headline uses the 20 easiest genes. Whole-panel accuracy on new slides is low.
- On the released files, **GHIST scores above GHIST+**, so the paper's ranking does not reproduce exactly.
- **Our GHIST+ retrain with an open encoder** (Phikon-v2 instead of the gated UNI2-h; authors' code and config, 30 epochs, A100):
  - training-log validation, top-20 SVGs: 0.687 (epoch 1) → 0.722 (epoch 3) → **0.732 (epoch 10)**
  - training-log validation continues: 0.739 (epoch 13), 0.748 (epoch 14), 0.737 (epoch 19)
  - final epoch-30 score, rescored with our checked script on the same 12,719 cells: **[PENDING, about 2 h]**

## Slide 4: Limitations of GHIST+ and our goal
**Technical**
- **Hard on/off export:** only the top-K cells per gene keep a value; the rest become 0 (86% zeros on Breast5K). That creates ties and destroys ranking.
- **Miscalibrated:**
  - predicted spread is 3× the measured spread
  - "30% chance expressed" isn't 30% (calibration error 0.12)
- **No per-cell uncertainty.**
- **Core correction (RRC):** a linear shift of the gated output, tuned on random cells (mildly optimistic).
- **Low whole-panel accuracy** on new slides (0.07–0.11).
- **Costly to retrain** (72–120 GPU-hours) and needs gated UNI2-h.

**Biological**
- Expression is **graded** (activation, partial EMT, infiltration gradients), but the export forces yes/no.
- **Xenium misses transcripts**, so a measured zero is often not true absence.
- Each tumour runs **its own programmes** (receptor status, immune infiltration); RRC applies one global shift.
- **Rare immune cells** and **rare signalling genes** (chemokines in under 1% of cells) are poorly predicted.
- The **measured ground truth is contaminated:** tumour-marker RNA appears 6.1× more often in immune cells touching tumour than in those 30–50 µm away (transcript spillover; our measurement).

**Goal:** improve on GHIST+ without retraining it.
1. **No RNA (zero-shot):** beat the GHIST+ export.
2. **A few measured cores:** beat the **authors' own correction (RRC)**, not just raw GHIST+.

**Rules:**
- cores placed without looking at RNA
- 200 µm buffer
- protocol frozen before each run
- Breast5K for development; ATERA and the pan-cancer model for confirmation
- results reported either way

## Slide 5: Approach 1: residual adapters (v1 → v3)
**Technical**
- **v1:** ridge regression from all 258 GHIST+ predictions to the residual, trained on 5k cells.
- **v2:** the same ridge with 5k–150k cells, plus an MLP on the Mac GPU.
- **v3 (GHIST+-Adapt):** ridge from realistic cores (4 × 1 mm²), per-gene fallback, spread matching, conformal intervals.

**Biological idea:** the slide's own measured cells reveal which genes co-vary in *this* tumour.

**Better than GHIST+**
- **Technical:**
  - v2 with 150k cells: 0.071 → 0.126
  - v3: better on 9–11 of 12 paper metrics on 3 slides (ATERA 0.106 → 0.159; BreastCancer2 0.199 → 0.271)
  - sharper spatial maps (SSIM 0.376 → 0.677 on Breast5K)
- **Biological:** cell-type programme scores improve on all 3 slides (Breast5K 0.247 → 0.278, ATERA 0.398 → 0.439, BreastCancer2 0.376 → 0.433), and gene–gene structure (CMD) improves.

**Worse**
- v1 (5k cells) is worse than GHIST+ (0.071 → 0.066).
- The v2 MLP failed (0.064): squared-error loss on mostly-zero RNA plus clipping at 0 created ties.
- **v3 lost to the authors' RRC in 18/18 head-to-heads** (ATERA 0.157 vs 0.166; BreastCancer2 0.265 vs 0.274). The authors ship the same idea with better-tuned shrinkage.
- Magnitudes collapse (SD ratio 3.08 → 0.17), and error on expressing cells rises on Breast5K and ATERA.

**When to use:** don't. It's superseded by Approach 2. It's kept as the negative result that made the authors' RRC the baseline to beat.

## Slide 6: Approach 2: soft output + core correction + uncertainty
**Technical**
1. **Soft output (no RNA):** q · ReLU(μ), GHIST+'s own probability × amount, instead of the top-K export.
2. **RRC-soft:** the authors' RRC applied to the soft output.
3. **HAC:** correct "is it on?" and "how much?" separately from cores.
4. **Uncertainty:** calibrated P(expressed) plus 90% intervals for the amount in expressing cells.

**Biological reason:** expression is graded and detection imperfect, so the probability is the right output. Cores reveal this tumour's programmes.

**Better than GHIST+**
- **Soft output, no RNA, 27/27 runs:**
  - +25% Breast5K (0.070 → 0.088), +35% ATERA (0.109 → 0.146), +31% pan-cancer model (0.054 → 0.071)
  - 216/258 and 358/413 genes improve
- **Biological:**
  - **all 8 cell-type marker programmes improve on both slides** (malignant +0.081)
  - largest gains in luminal/ER genes (CCND1 0.44 → 0.53, GATA3 0.47 → 0.56), epithelial identity (CDH1 0.41 → 0.59), T/lymphoid genes (RORC 0.16 → 0.33) and vessel/smooth muscle (MYH11 0.32 → 0.43)
- **RRC-soft** beats the authors' RRC in **27/27 runs** (+12–15%).
- **HAC:**
  - beats the authors' RRC in 20/27 runs
  - largest biological gains for **scattered immune cells** (ATERA plasma +0.022, myeloid +0.018, B +0.017, macrophage +0.017, T +0.015)
  - best gene–gene structure of all methods
- **Uncertainty:** calibration error 0.120 → 0.005 (Breast5K) and 0.058 → 0.003 (ATERA); 90% coverage with intervals 20–48% narrower.

**Worse**
- The soft output alone **blurs spatial maps** (SSIM 0.376 → 0.279) and gene–gene structure.
- **HAC fails on the pan-cancer model.**
- All core corrections, the authors' included, err more on the amount in expressing cells (0.71 vs 0.66).
- **No gain** for rare chemokines (CCL20, CCL8) or myeloid cells (+0.004–0.006).
- **A few genes get worse** (CD36, TSPAN8, ITGB1).
- **Failed add-ons:**
  - trust gate for unseen tissue (−0.011)
  - Phikon image features (worse in 7/9 runs)
  - atlas and niche priors (no gain)

**When to use**
- **No RNA:** the soft output, for gene ranking and marker scores.
- **Sharp maps or yes/no calls without RNA:** the original GHIST+ export.
- **With cores:** RRC-soft for accuracy; HAC when co-expression matters.
- **When you need per-cell confidence:** the uncertainty layer.

## Slide 7: Approach 3: BioHAN + biology-driven extensions (HPC)
**Technical: BioHAN** (neural network trained on core cells; 20–100 s per fit on an A100)
- **Two-part output:** "is it on?" × "how much?". It's strictly positive, which fixes the v2 MLP failure.
- **64-factor gene-programme bottleneck.**
- **Niche inputs:** same-type neighbours within 30 µm, and the cell-type mix within 50 µm.
- **Rare-cell-type up-weighting.**

**Biological reasons:**
- detection and expression level are separate processes
- genes act in programmes
- the microenvironment shapes cell state
- rare immune cells matter

**Results (9 runs per combination)**

| | Raw GHIST+ | RRC-soft (best previous) | **BioHAN** | BioHAN beats RRC-soft (CI > 0) |
|---|---|---|---|---|
| Breast5K: gene PCC | 0.071 | 0.112 | **0.118** | **9/9** |
| ATERA: gene PCC | 0.109 | 0.186 | **0.192** | **5/9** |
| Pan-cancer model: gene PCC | 0.054 | 0.082 | **0.092** | **9/9** |
| Spatial maps (SSIM), Breast5K / ATERA / pan | 0.37 / 0.66 / 0.40 | 0.66 / 0.79 / 0.59 | **0.78 / 0.81 / 0.79** | – |

**Better than GHIST+**
- **Technical:**
  - +66% / +76% / +70% gene PCC over raw GHIST+
  - top-20 SVG PCC +0.10 to +0.23
  - **sharpest spatial maps of any method**
  - lowest overall error
- **Against the authors' correction:** it **passes the pre-declared rule on all 3 combinations**.
- **Biological:** tissue-structure genes (top-50 SVG maps) are reconstructed much better. Gene–gene structure improves over raw GHIST+ on ATERA (0.404 → 0.314) and pan-cancer (0.522 → 0.499).

**Worse**
- **Gene–gene co-expression** is 6–18% worse than RRC-soft, and worse than raw GHIST+ on Breast5K (0.445 → 0.501).
- **Amount in expressing cells:**
  - error is higher than raw GHIST+ (0.664 → 0.713 on Breast5K)
  - predictions are compressed (SD ratio 3.1 → 0.17)
- The ATERA margin is thin (5/9).
- **Ablation: none of the biology components passes its pre-declared test** (removing a component must lower PCC with CI below 0 in most runs on both Breast5K and ATERA).
  - **Niche inputs:** formal test fails (better in 1/9 ATERA runs, worse in 2/9 Breast5K runs, 0/9 either way on pan-cancer).
  - **Gene-programme bottleneck:** small mean gain on ATERA (+0.003) but significant in only 1/6 runs; neutral on Breast5K and pan-cancer.
  - **Rare-cell weighting:** no measurable effect on any of the 3 combinations.
  - **Pan-cancer ablation (6/6 runs):** every component within ±0.001 PCC of the full model, CI never excludes 0 against it.
  - The gain comes mainly from the two-part neural output.

**Other biology extensions tested (linear)**
- **Transcript spillover:**
  - real contamination (6.1×)
  - modelling it didn't help (+0.0013 vs control +0.0019)
  - kept as an evaluation caveat
- **Same-type niche smoothing:** a null result. The cores chose no smoothing in 27/27 runs.
- **Gene-programme projection (Prog-RRC):**
  - PCC +0.010 (Breast5K), +0.010 (ATERA), +0.019 (pan-cancer, 0.082 → 0.102), CI above 0 in **27/27 runs**
  - but co-expression is 16–22% worse, so it **fails its pre-declared rule** (which required co-expression to improve)
  - spatial maps also blur on Breast5K (−6% SSIM) and pan-cancer (−10%)

**When to use**
- **With cores, best gene accuracy and spatial maps:** BioHAN.
- **A simple linear alternative for gene accuracy:** Prog-RRC (secondary finding).
- **Gene networks, co-expression or absolute expression levels:** HAC, the authors' RRC, or raw GHIST+.
- **No cores:** BioHAN can't be used; use the Approach 2 soft output.

**Pending overall**
- GHIST+ open-encoder final score vs 0.744 (training at epoch 20/30 at 17:35 IST; about 2 h left)
- our soft-output fix tested on that retrained GHIST+ (after training)
