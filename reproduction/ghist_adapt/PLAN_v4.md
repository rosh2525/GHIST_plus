# Plan v4: genuinely improving on GHIST+ (and on the authors' own core correction)

> **Status (2026-11-09).** Each package was run under a protocol frozen beforehand; amendments are dated in its config.
> - **P0 done.** The authors' RRC beats my old v3 adapter in 18/18 runs. Its random-cell tuning is mildly optimistic.
> - **P4b FAIL** (`config_v6.json`). Phikon image features, extracted for all 565k cells on the M3 GPU, made HAC slightly worse (CI < 0 in 7/9 runs, > 0 in 0/9).
> - **P4a done on Kaggle (pan-cancer checkpoint).** The soft output replicates (+31%). RRC-soft beats the authors' RRC in 27/27 runs across all three model/slide combinations. HAC fails the pre-declared rule on this model (1/3).
> - **P1 PASS on Breast5K/ATERA** (`config_v4.json`). New zero-shot soft hurdle expectation: +25% / +35% median PCC. HAC beats the authors' RRC on both slides.
> - **P2 FAIL** (`config_v5.json`). The trust gate helps locally on Breast5K with biased cores but costs 0.011 all-cell PCC. On ATERA, HAC does not hurt unfamiliar tissue, so there is nothing to gate. Side finding: with biased cores the authors' RRC drops below zero-shot on unfamiliar tissue and HAC does not.
> - **P3 PASS on the amended (hurdle-aware) rule** (`config_v7.json`). P(expressed) calibration error falls 19–24× (ECE 0.120→0.005 Breast5K, 0.058→0.003 ATERA). Conditional magnitude intervals reach 90.4% / 89.9% coverage on expressing cells and are 48% / 20% narrower than around raw GHIST+. The original marginal per-gene intervals turned out to be degenerate on zero-inflated data (they cover 0% of expressing entries); that is reported as a finding.
> - **P5 (biology priors, development probes): no gain.** Atlas co-expression projection (Wu et al. 2021), niche features and cell-type magnitude did not improve HAC on Breast5K; not taken to confirmation (`probes/`, README).
> - **No DGX available. P4 moves to free Kaggle compute** (`kaggle/README.md`, `kaggle/ghist_adapt_kaggle.ipynb`).
> - **P4 ready for the DGX** (`dgx/README_DGX.md`). P4a is a second-checkpoint confirmation (pan-cancer model, Breast5K), CPU only. P4b is the UNI2-h image-aware HAC; it needs the 10x H&E links filled in `dgx/he_sources.json`, and `config_v6.json` must be frozen before it runs.

Written 2026-11-09, after finding that v3's core correction duplicates the authors' released
`tools/one_shot_rrc.py` (shrunk ridge residual correction, "RRC"). In head-to-head runs RRC was slightly
better than v3. From now on **the baseline to beat is the authors' RRC on identical cores**, not frozen GHIST+.

## 0. Novelty check (done before planning)

Searched the official code (`model/`, `utils/`, `tools/`, `dataio/`, `train*.py`):

| Idea | In the authors' code? | Status for us |
|---|---|---|
| Ridge residual correction from measured cores | **yes**: `tools/one_shot_rrc.py` | baseline only |
| Core / TMA-punch selection | **yes**: VQ selector `utils/tma_select.py` | not claimed |
| Training with core cells of the target slide | **yes**: `train_tma_select.py` | baseline if runnable on DGX |
| Correction using the hurdle outputs (occurrence prob. + magnitude) | no; RRC uses point predictions only | **candidate (P1)** |
| Tissue-support / out-of-distribution gating of the correction | no (the only "trust" is ECRM neighbour weighting) | **candidate (P2)** |
| Prediction intervals / conformal calibration | no | **candidate (P3)**; prior art TISSUE (Sun et al. 2024), cite it |
| Core correction using UNI2-h image features | no; RRC sees GHIST+ outputs only | **candidate (P4)** |

A literature check is still required before any novelty claim. "Not in their code" is not "new".

## 1. Limitations targeted

| Paper limitation | Addressed by |
|---|---|
| B1: core-to-WSI anchoring is a global post-hoc correction applied to every cell, even tissue states the cores never saw | P2 |
| B2: magnitudes miscalibrated on new slides; RRC/v3 compress spread and clip to 0 (ties destroy ranking of sparse genes) | P1 |
| B3: no per-cell uncertainty | P3 |
| B5/B7: adapting needs a 72–120 h retrain; frozen encoder | P4 |
| B8: sparse genes ranked poorly (they dominate the median PCC) | P1, P4 |
| A1/A2: one test slide, top-SVG headline, no CIs | evaluation design (all genes + SVG groups, 3 slides, block-bootstrap CIs, frozen protocol) |

## 2. Work packages

**P0: baselines and protocol (Mac, about 1 day)**
- Finish the head-to-head GHIST+ vs authors' RRC (their functions, default config) vs v3 on identical cores, all 3 slides.
- Diagnostic: the authors' RRC picks α and blend on a random 20% of ROI cells. Measure how optimistic that is versus leave-one-core-out.
- Freeze `config_v4.json` before any P1–P3 run: cores, buffer, slides, metrics, decision rule.

**P1: hurdle-aware core correction (Mac, 3–4 days). Main candidate for better accuracy.**
GHIST+ predicts each gene as P(expressed) × magnitude, but RRC corrects only the final number and clips at 0.
- Correct the two parts separately from core cells:
  - occurrence: multigene logistic recalibration from [occurrence logits, predictions, cell-type probabilities]
  - magnitude: ridge on expressing core cells
- Output E[y] = π̂ · m̂. It is strictly positive, so no clipping ties.
- Hypothesis: better ranking of sparse genes (median PCC) and better magnitude on expressing cells, where RRC/v3 lose.
- Data: Breast5K and ATERA exports contain the hurdle outputs. The BreastCancer2 CSV does not, so P1 confirmation there needs the checkpoint re-run (P4) or uses ATERA as the confirmation slide.

**P2: support-aware (trust-gated) correction (Mac, 2–3 days)**
- For each cell, measure how far its GHIST+ state (prediction PCs + cell-type probabilities) is from the core cells' states (kNN distance).
- Shrink the correction λ(distance) → 0 outside the cores' support.
- Learn λ from leave-one-core-out: how the correction's benefit decays with distance.
- **Stress test** (the point of P2): cores placed only in one tissue state (e.g. tumour-rich), evaluated on the others. Hypothesis: RRC degrades unseen states and P2 does not. Also report the unbiased-core setting, where P2 should be neutral.

**P3: spatially valid uncertainty (Mac, about 2 days)**
- Conformal intervals around the corrected prediction from leave-one-core-out residuals.
- Mondrian conditioning on predicted cell type and on P2's support distance, so intervals widen away from core-like tissue.
- Report coverage and width on the withheld slide, per cell type and per distance band; compare with the same procedure on RRC and on raw GHIST+.
- This is a new capability, not an accuracy gain.

**P4: image-aware core adaptation (DGX, 1–2 weeks). Highest potential, highest setup cost.**
- Get UNI2-h access and the public H&E for the slides. Reconstruct the released checkpoint.
- Extract per-cell UNI2-h embeddings at the existing cell coordinates.
- Train a light head on [UNI2-h embedding, GHIST+ hurdle outputs] from core cells.
- Compare against RRC, P1 and, if runnable, the authors' `train_tma_select.py` with the same cores.
- Also regenerates BreastCancer2 predictions with hurdle outputs, for P1 confirmation.

## 3. Rules (so results are credible)
- Cores are chosen without RNA. Hyperparameters, gates, λ and interval radii are chosen by leave-one-core-out on core cells only. Evaluation is on the withheld slide with a 200 µm buffer.
- An extension counts as an improvement only if it beats **RRC** on the primary metric (median per-gene PCC, all genes) on at least 2 of 3 slides, with the block-bootstrap 95% CI of the difference excluding 0, and without degrading the guardrails (RMSE, error on expressing cells, SSIM) by more than 2%. Every result is reported either way.
- Breast5K is the development slide. Confirmation happens on slides not used to tune the extension.

## 4. Timeline
| Week | Work | Output |
|---|---|---|
| 1 | P0, P1, P2, P3 on the Mac | head-to-head table; P1/P2/P3 results vs RRC |
| 2–3 | P4 on DGX (needs UNI2-h access; request now) | image-aware adapter vs RRC vs TMA-select |
| 4 | write-up | report and slides |

**Minimum complete project without the DGX:** P0 + P1 + P2 + P3.

## 5. Expected outcomes (honest)
- P1: the best chance of beating RRC on accuracy on the Mac; not guaranteed.
- P2: likely neutral when cores cover the slide; useful when they don't. Its value is safety under extrapolation, which the paper's core-to-WSI pitch implicitly assumes away.
- P3: certain to deliver calibrated intervals; it is a capability, not an accuracy gain.
- P4: most likely to give the largest gain; depends on access and setup.
- If nothing beats RRC: the project becomes "an independent, leak-free multi-slide evaluation of GHIST+'s core correction, with uncertainty and failure analysis". That is still a legitimate BTP, but not a new method.
