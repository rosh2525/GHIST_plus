# Improving GHIST+ predictions: what is new, what is not

**Setting.** This is the paper's core-to-WSI scenario: H&E everywhere, RNA measured only in a few TMA-like cores of the slide, and GHIST+ predictions everywhere. All work here uses the authors' released GHIST+ predictions; GHIST+ itself is not retrained.

## Summary

1. **New, and the main result: GHIST+ discards information in its own export, and recovering it needs no RNA.** GHIST+'s exported prediction applies a hard top-K gate (`model/hurdle_distribution.py::deterministic_hurdle_prediction`); on Breast5K, 86% of cell-gene values are exactly 0. Using the model's own soft expectation, q·ReLU(μ), from the released hurdle outputs improves median per-gene PCC on **three model/slide combinations**, better in **27/27 runs**:
   - **Breast5K, multi-breast model (development):** 0.070 → 0.088 (+25%)
   - **ATERA, multi-breast model (untouched):** 0.109 → 0.146 (+35%)
   - **Breast5K, pan-cancer model (untouched second checkpoint, 700 genes; run on Kaggle):** 0.054 → 0.071 (+31%)
   Top-20 SVG PCC and RMSE also improve.
2. **With a few measured cores, correcting the soft output beats the authors' own correction.** Running the authors' released RRC (`tools/one_shot_rrc.py`, their code and defaults) on the soft score instead of the export ("RRC-soft") beats their default in **27/27 runs** across all three combinations: +12% (Breast5K), +12% (ATERA), +15% (pan-cancer), every block-bootstrap 95% CI above 0. RRC-soft was listed in the frozen protocol (`config_v4.json`) as an ablation; it is promoted here *after* seeing the results, because it was consistently best.
3. **HAC, our pre-declared primary hurdle-aware correction, is weaker than RRC-soft.** It beats the authors' RRC in 20/27 runs and passes the pre-declared rule on Breast5K and ATERA, but **fails it on the pan-cancer model** (CI > 0 in 1/3 primary repeats). HAC's specific two-part design adds an edge only on ATERA.
4. **Uncertainty (P3): passed.** HAC makes GHIST+'s P(expressed) well calibrated (calibration error 19–24× lower), and gives 90%-coverage magnitude intervals for expressing cells that are 20–48% narrower. Plain per-gene conformal intervals turned out to be degenerate on zero-inflated data.
5. **Failed or retired.**
   - Image features from an independent pathology encoder (Phikon, P4b) made HAC slightly worse (0/9 runs better, 7/9 worse).
   - The support-aware gate (P2) failed its pre-declared test.
   - Three biology priors (atlas co-expression, niche, cell-type magnitude) gave no gain in development probes.
   - My earlier v3 adapter lost to the authors' RRC in 18/18 runs.
6. **Diagnostic.** The authors' RRC tunes on a random 20% of core cells. That score exceeds the leave-one-core-out score in all 9 slide × config combinations, a mild optimism.

## What is new compared with the authors' code (checked against the official repo)

| Idea | In the authors' code? | Here |
|---|---|---|
| Ridge residual correction from measured cores | **yes** (`tools/one_shot_rrc.py`) | baseline |
| Core/TMA punch selection | **yes** (VQ, `utils/tma_select.py`) | not claimed |
| Soft hurdle expectation instead of the top-K gated export | no; export and evaluation use the hard gate | **new (zero-shot)** |
| Correcting occurrence and magnitude separately from cores (HAC) | no; RRC corrects the gated point prediction | **new** |
| Uncertainty intervals | no | v3 tool (conformal); prior art TISSUE (Sun et al. 2024) |

"Not in their code" is not proof of novelty. A literature check is still needed before any publication claim.

## Methods

```
GHIST+ released outputs per cell × gene:  q = occurrence probability,  μ = signed magnitude,  export = topK(q) · ReLU(μ)

GHIST+-Soft  (zero-shot)   S = q · ReLU(μ)
RRC-authors  (baseline)    authors' one_shot_rrc on the export (their functions and defaults)
RRC-soft     (ablation)    the same authors' RRC applied to S
HAC          (ours)        q̂ = clip(q + ridge([S, logit q]) , 0.001, 0.999)    occurrence correction, all genes as inputs
                           m̂ = max(a_g + b_g · ReLU(μ), 0.01)                  per-gene magnitude fitted on expressing core cells
                           output = q̂ · m̂                                      strictly positive, so no clipping ties
```

Everything is fitted on core cells only, with ridge α chosen by leave-one-core-out CV. It runs in seconds on CPU. Protocol: `config_v4.json`, frozen after the Breast5K zero-shot comparison and before any ATERA computation. Cores and evaluation are identical to v3: disk cores placed without RNA, a 200 µm buffer, and evaluation on the withheld rest of the slide. Code: `hurdle.py`, `p1_run.py`, `p1_report.py`.

## Results (`outputs_v4/summary.md`, `outputs_v4/figures/p1_scoreboard.png`)

**Zero-shot, all cells (no RNA):**

| slide | median PCC export → soft | top-20 SVG | top-50 SVG | mean RMSE |
|---|---|---|---|---|
| Breast5K (dev) | 0.0704 → **0.0882** | 0.383 → **0.436** | 0.243 → **0.307** | 0.334 → **0.201** |
| ATERA (confirmation) | 0.1085 → **0.1463** | 0.471 → **0.591** | 0.395 → **0.455** | 0.374 → **0.284** |

**With cores, primary config (4 × 1 mm², mean of 3 repeats).** Columns: GHIST+ / Soft / authors' RRC / RRC-soft / **HAC**.

| metric | Breast5K (dev) | ATERA (confirmation) |
|---|---|---|
| median PCC, all genes | 0.071 / 0.088 / 0.102 / 0.112 / **0.110** | 0.106 / 0.143 / 0.169 / 0.189 / **0.188** |
| median PCC, top-20 SVGs | 0.380 / 0.434 / 0.476 / 0.489 / **0.489** | 0.458 / 0.586 / 0.684 / 0.690 / **0.682** |
| median Spearman | 0.068 / 0.065 / 0.070 / 0.078 / **0.077** | 0.098 / 0.098 / 0.111 / 0.130 / **0.133** |
| mean RMSE ↓ | 0.335 / 0.201 / 0.099 / 0.098 / **0.097** | 0.375 / 0.284 / 0.221 / 0.219 / **0.221** |
| RMSE, expressing cells ↓ | 0.664 / 0.639 / 0.705 / 0.707 / **0.712** | 0.737 / 0.733 / 0.738 / 0.739 / **0.742** |
| SSIM, top-50 SVGs | 0.376 / 0.279 / 0.632 / 0.657 / **0.699** | 0.666 / 0.647 / 0.785 / 0.789 / **0.774** |
| CMD, all genes ↓ | 0.445 / 0.577 / 0.421 / 0.425 / **0.402** | 0.408 / 0.507 / 0.264 / 0.291 / **0.319** |

**Pre-declared decision, HAC vs authors' RRC:** **PASS** on both slides. Mean better, CI above 0 in 3/3 repeats on each, and guardrails within 2%. One caveat: on Breast5K, HAC's SSIM is 9.5% *better* than RRC's, and on ATERA it is 1.4% worse.

## Biology (`outputs_v4/biology_p1.txt`, `outputs_v4/<slide>/zero_shot_per_gene.csv`)

**Why the gate is a biological problem.** GHIST+'s hurdle head models two things: whether a gene is "on" in a cell (occurrence, q), and how much (magnitude). The export turns q into a yes/no decision: only the top-K cells per gene are "on".

Tissue expression is graded. Cell states form continua (activation, partial EMT, infiltration gradients near tumour borders). Xenium detection is also probabilistic, so a measured zero is often a missed transcript rather than true absence. The soft output keeps the graded probability; the export turns it into a cliff.

**What the data shows (zero-shot, all cells):**
- **Nearly every gene improves:** 216 of 258 genes on Breast5K (8 worse) and 358 of 413 on ATERA (16 worse).
- **Common genes gain most in absolute PCC.** By quartile of fraction of cells expressing (Breast5K), the gain is +0.044 for the most common quarter and +0.006 for the rarest (ρ = +0.64; ATERA +0.40). In relative terms the rarest genes gain most (+73% vs +22%).
- **The rarest genes remain essentially unpredictable** (PCC about 0.02–0.03). Genes expressed in under 1% of cells, such as the chemokines CCL20 and CCL8, are not rescued, which suggests a limit of morphology-based prediction for transient or rare signalling states.
- **All 8 cell-type marker programmes improve on both slides (16/16):**
  - malignant +0.081 (Breast5K)
  - epithelial +0.041 / +0.017
  - plasma +0.039 / +0.042
  - macrophage +0.033 / +0.032
  - B +0.030 / +0.011
  - myeloid only +0.006 / +0.004
- **Largest single-gene gains:** luminal / ER-pathway genes (CCND1 0.44 → 0.53, GATA3 0.47 → 0.56); epithelial identity (CDH1 0.41 → 0.59, AGR3 0.43 → 0.59, TFAP2A 0.27 → 0.41 on ATERA); T-cell / lymphoid genes (TRBC1 0.17 → 0.26, RORC 0.16 → 0.33, RHOH 0.27 → 0.40); vascular / smooth muscle (MYH11 0.32 → 0.43, AQP1 0.36 → 0.47).
- **Genes that got worse** include CD36, C11orf96, TSPAN8 and ITGB1. For these, GHIST+'s occurrence probability carries little or misleading signal, so the gate was hiding noise.

**HAC vs the authors' RRC, by programme** (per-gene ΔPCC, primary config):
- **ATERA immune programmes gain most:** plasma +0.022, myeloid +0.018, B +0.017, macrophage +0.017, T +0.015.
- **Abundant epithelial genes:** RRC is as good (−0.007).
- **Breast5K:** +0.005 to +0.016 for every programme except B (−0.006).
- **Interpretation:** correcting *whether* a gene is on (occurrence) matters most for markers of scattered immune cells, which a linear shift of the gated point prediction (RRC) handles poorly.

These are correlations with measured RNA, not causal biology.

## Where it is worse (report these)
- **Gene–gene structure on ATERA.** CMD is worse for HAC (0.319) than for the authors' RRC (0.264).
- **Top SVGs and SSIM on ATERA.** HAC does not beat RRC there: top-20 SVG PCC 0.682 vs 0.684, SSIM 0.774 vs 0.785. HAC's gain is concentrated in the many weakly expressed genes, which drive the median.
- **The zero-shot soft score alone** gives worse SSIM on Breast5K (0.376 → 0.279) and worse CMD on both slides. It ranks cells better per gene but blurs spatial maps and co-expression structure. Use it for per-gene ranking, not for maps.
- **Expressing cells.** All core corrections, including the authors' RRC, have higher error on expressing cells than GHIST+ on Breast5K, and they compress the spread of expression (spread ratio about 0.2).
- **Breast5K single-core config.** HAC only ties RRC there, and RRC-soft is best.

## Limits
- **Two slides.** BreastCancer2's export has no hurdle outputs, so P1 confirmation there needs the checkpoint re-run on the DGX.
- **Development slide.** Breast5K was used for development, and the zero-shot score was chosen among 4 candidates on it. ATERA is the clean test.
- **Simplified expectation.** The soft score is q·ReLU(μ), not the exact truncated-normal mean, because the model's σ is not exported.
- **Within-slide only.** Core corrections use measured RNA from the same slide. The zero-shot fix does not.
- **ATERA units.** The coordinate scale is assumed (see v3), which affects only the mm² labels.

## P2: support-aware trust gate. Result: FAIL (`config_v5.json`, `outputs_v5/summary.md`)
**Idea.** Shrink HAC's cell-specific correction for cells unlike any core cell, so that a correction learned from cores missing some tissue states is not applied to those states. **Stress test:** cores placed only in tumour-rich or only in stroma-rich regions, chosen without RNA.

| | Breast5K (dev) | ATERA (confirmation) |
|---|---|---|
| Biased cores make HAC worse than zero-shot on unfamiliar cells? | yes (tumour-only 0.089 vs 0.095) | **no** (HAC 0.190 vs zero-shot 0.171) |
| Gate fixes those cells (CI > 0) | 2 of 3 repeats per biased setting | n/a |
| Cost on all cells, unbiased cores | **−0.011 median PCC** (allowed −0.002) | −0.001 |

The gate does not meet the pre-declared bar. HAC's adjustments for unusual cells carry real signal, and on ATERA there was no extrapolation damage to prevent. The gate design needed two dated amendments after smoke tests on the development slide: one fixed mixing outputs on different scales, the other fixed gating the global recalibration. Both are recorded in the config.

**Side finding:** with tumour-only cores on ATERA, the authors' RRC drops *below* zero-shot on unfamiliar tissue (0.150 vs 0.171), while HAC stays above it (0.190). HAC is more robust to badly placed cores than RRC.

## P3: uncertainty. Result: PASS on the hurdle-aware rule (`config_v7.json`, `outputs_v7/summary.md`)
**First finding: plain per-gene conformal intervals are degenerate here.** They reach the nominal 91–92% coverage, but cover **0% of expressing entries** for HAC and RRC, because more than 90% of measured values are zero and the radius learns to cover zeros. Any headline "calibrated intervals" claim on this kind of data needs this check.

**Hurdle-aware uncertainty** (dated amendment after one development-slide smoke run, before any full or ATERA run) reports two things:

| unbiased cores | Breast5K (dev) | ATERA (confirmation) |
|---|---|---|
| Calibration error of P(expressed): raw GHIST+ q → HAC q̂ (ECE ↓) | 0.120 → **0.005** | 0.058 → **0.003** |
| Brier score ↓ | 0.0273 → **0.0104** | 0.0599 → **0.0451** |
| "If expressed, how much": 90% interval coverage on expressing entries | **0.904** | **0.899** |
| Interval width: HAC vs raw GHIST+ magnitude | **0.40** vs 0.77 (−48%) | **0.80** vs 0.99 (−20%) |

With tumour-only cores, coverage stays at 0.89–0.91 and ECE at 0.006–0.007. The same procedure on GHIST+'s raw outputs is miscalibrated: a stated 30% chance of expression does not come true 30% of the time.

## P4b (image features, Phikon on the M3 GPU): FAIL (`config_v6.json`, `outputs_v6/breast5k/`)
UNI2-h is gated, so we used **Phikon** (owkin/phikon, ViT-B/16, open, Owkin non-commercial licence), an encoder independent of the UNI2-h that GHIST+ uses internally.
- **H&E source and alignment:** 10x post-Xenium H&E with the inverted 10x alignment matrix, verified by a haematoxylin-contrast QA plus a visual overlay.
- **Features:** 112 µm crops at 0.5 µm/px; embeddings for all 564,723 Breast5K cells in 4 h on the M3 GPU (`extract_features.py`).
- **Test:** HAC with and without 64 image PCs, on the same cores as P1.

| Config | HAC | HAC + image | CI(image − HAC) > 0 / < 0 | genes ↑ / ↓ |
|---|---|---|---|---|
| 4 × 1 mm² diverse (primary) | **0.110** | 0.104 | 0/3 / 3/3 | 26 / 87 |
| 4 × 1 mm² random | **0.105** | 0.102 | 0/3 / 1/3 | 19 / 67 |
| 1 × 3.5 mm² | **0.101** | 0.100 | 0/3 / 3/3 | 10 / 67 |

**Reading.** Image features only improved co-expression structure slightly (CMD 0.402 → 0.392). The likely reason is that GHIST+ already captures most of the morphology signal, and about 27k core cells cannot support 64 extra features without overfitting. Development-slide evidence only: no other slide has both hurdle outputs and a public H&E.

## P5 (biology priors, development-slide probes only): no gain
Three biology-motivated additions were probed on Breast5K (primary config, seeds 0–1) before any protocol was frozen. None justified a confirmatory run (`probes/`):

| Idea | Result |
|---|---|
| Project the soft output onto co-expression programmes from the Wu et al. 2021 breast scRNA atlas (GSE176078, 100k cells; k = 5–80 PCs) | **worse**: median PCC 0.088 → 0.066–0.087; CMD 0.58 → 0.59–0.70. Dissociated-cell correlations are much stronger than Xenium's in-tissue correlations |
| Atlas cell-type profiles weighted by GHIST+'s cell-type probabilities | worse than the soft output (0.059 vs 0.089) |
| HAC + atlas programme scores (+ neighbourhood / niche averages) | no gain (0.1069 → 0.1064; with niche 0.1013) |
| HAC + cell-type-specific magnitude and occurrence | +0.0003 to +0.0006; expressing-cell RMSE unchanged (0.711) |

**Reading:** HAC already uses most of the information in GHIST+'s outputs. Expressing-cell RMSE is structurally higher for any expected-value output (P(on) × amount); the conditional amount m̂ is reported separately (P3). New information, such as image features, is the next lever.

## Next (`PLAN_v4.md`; free compute: `kaggle/README.md`; DGX runbook `dgx/README_DGX.md`)
- **P4a (CPU, run first):** P1–P3 on a second released checkpoint (pan-cancer model, same Breast5K slide).
- **P4b (GPU):** UNI2-h image-aware HAC. It needs the 10x H&E and alignment links in `dgx/he_sources.json`; freeze `config_v6.json` before running it.
- **Not packaged:** GHIST+ re-inference and the `train_tma_select.py` head-to-head. Both need the authors' unreleased preprocessed inputs.

## Reproduce
```sh
cd reproduction/ghist_adapt            # environment: requirements.txt; export OPENBLAS_NUM_THREADS=2 when running in parallel
python slides.py breast5k atera        # matched caches
python hurdle.py breast5k atera        # hurdle-output caches (q, μ) from the released .npz
python p1_run.py breast5k atera        # P0 + P1 (about 20 min per slide on an M3)
python p1_report.py                    # outputs_v4/summary.md + figure
python run_v3.py breast2 atera --headline-only --out outputs_rrc   # v3 vs authors' RRC head-to-head
```

## Archive: v3 (GHIST+-Adapt)
`config_v3.json`, `run_v3.py`, `report_v3.py`, `strata_v3.py`, `outputs/`, `outputs_rrc/`.
- **What the v3 runs measured:** GHIST+-Adapt beat the frozen GHIST+ export on 9–11 of 12 metrics on 3 slides (`outputs/summary.md`). Its cell-type and biology analyses are in `outputs/strata.md` and the `programme_weight_map_*.png` figures.
- **Why it is retired:** it lost to the authors' RRC in all 18 head-to-heads (`outputs_rrc/`).
- **Grid status:** the v3 learning-curve grid was paused when the method was superseded. Breast5K has 37 of 60 configs, ATERA 26 of 36, and BreastCancer2 is complete.
