# Review of Gemini Flash feedback
Date: 31 October 2026
Scope: mathematical and implementation review of the planning specification; no model training or biological experiment.

## Decision
Adopt the useful safeguards. Do not accept every explanation or proposed fix.

| Feedback | Assessment | Action |
|---|---|---|
| Zero-containing intervals provide weak occurrence supervision | Correct limitation, not an incorrect probability formula | Separate magnitude and occurrence flags; add a magnitude-only pilot and stratified occurrence/sensitivity checks |
| Gradient vanishes whenever mu<u | Incorrect as stated | Record exact gradient and dependence on sigma/tail probability |
| Cap occurrence q using SPLIT weights | Not justified by available evidence | Do not adopt; lineage mixture weights are not calibrated per-gene occurrence probabilities |
| Detach the conditional tail | Does not establish identification | Do not adopt as a default fix |
| Small intervals create large loss offsets | Correct for raw interval mass | Add development-fixed minimum-width/numerical checks and optional average-density diagnostic |
| Width normalization fixes sigma incentives | Incorrect for fixed endpoints | It adds a parameter-independent constant; gradients do not change |
| Partial-label composition can erase immune cells | Correct and material | Use fully trusted composition units or disable aggregate composition supervision consistently; threshold alone is insufficient |
| Keep auxiliary masks out of predicted ECRM probabilities | Correct | Add a fixed-input/weights forward-invariance acceptance test |
| Match raw/alternative target scales | Correct | Validate source units; no independent renormalization of corrected totals during interval construction |
| L5/L6 role swap holds out a scanner image group | Correct with qualifications | Predeclare secondary split; no claim of unseen scanner hardware/site/ancestry |

## Occurrence calculation
For interval [0,u], define T=P(Z>u | Z>0), q=sigmoid(s).
L=-log(1-q*T).
dL/ds=q*(1-q)*T/(1-q*T).

The derivative is positive for interior q and positive tail probability. It can be negligible when T is small, but mu<u alone does not imply zero gradient.

Analytic numerical checks (not biological data):
- q=0.99, mu=1, u=2, sigma=0.2: derivative approximately 2.84e-9.
- Same q,mu,u with sigma=2: derivative approximately 0.007913.

High occurrence probability alone does not guarantee a positive exported value: the GHIST+ cohort gate also depends on other cells' rankings and the positive magnitude. Nevertheless, weak occurrence supervision is a meaningful risk and must be evaluated separately.

## Interval width
For a narrow positive interval of fixed width Delta,
-log P([l,u]) = -log f(y) - log Delta + o(1).

Adding log Delta produces average-density loss and the positive point-density limit, but changes neither mu nor sigma gradients because Delta is fixed. Finite-difference checks for widths 0.2,0.02,0.002 confirmed identical gradients with and without this constant.

Mixed exact-density and interval-probability contributions are not inherently invalid; standard censored likelihoods contain both. Our ranges are heuristic weak labels, however, not verified censoring bounds.

Do not apply continuous-width normalization to intervals containing the discrete zero atom. Choose minimum widths/numerical tolerances on development data and freeze them. The review's 0.10 is an example, not an established optimum.

## Composition and graph correction
The current training implementation aggregates composition at batch level and maintains running class statistics. Fixing a tile target alone is not enough. Predicted aggregate, label histogram, cell-count weights and running statistics must use the same eligible population.

An 80% coverage threshold can still remove every cell of a 10% minority lineage. Downweighting this biased target does not make it unbiased. Initially use fully trusted tiles for aggregate composition losses; if too few exist, disable those auxiliary losses consistently while leaving image/graph nodes and forward computation intact.

Relevant inspected code:
- train.py around 2210–2263: composition aggregation and running class statistics.
- model/framework.py around 894–903: predicted cell-type probabilities passed to ECRM.
- dataio/dataset_input.py around 400–405: expr_scale * log1p applied to supplied table values.

## Data correction
The loader's log1p transform does not prove upstream CSVs contain raw counts. Record upstream units and construction provenance. Use a common raw-count basis or identical declared factors for observed and alternative targets. Library-normalized relative expression remains legitimate for separately defined evaluation, but is not the same absolute-count quantity.

Secondary acquisition-group split:
training L1/L3/L6; validation L5; test L7/L8.
Training and evaluation use different listed scanner-image groups. Validation and test still share group2, and donor/acquisition effects remain confounded with only two groups. This tests sensitivity to those acquisition groups, not broad cross-scanner generalization.

## Sources
- [Current plan](GHIST_PLUS_IMPLEMENTATION_PLAN.md)
- [Pinned training code](https://github.com/SydneyBioX/GHIST_plus/blob/c5c895174c19903cc27d2661934b9b726bdb9306/train.py)
- [Hurdle implementation](https://github.com/SydneyBioX/GHIST_plus/blob/c5c895174c19903cc27d2661934b9b726bdb9306/model/hurdle_distribution.py)
- [SciPy mixed censored likelihood documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.fit.html)
- [SciPy zero-width and interval conventions](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.CensoredData.html)
