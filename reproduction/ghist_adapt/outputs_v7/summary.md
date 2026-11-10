# P3: uncertainty (config_v7.json, incl. amendment 1)

## Breast5K (dev)

**Marginal per-gene intervals (original rule):** median coverage (all entries / expressing entries), mean width

| setting | GHIST+ | GHIST+-Soft | RRC-authors | HAC | HAC-mondrian |
|---|---|---|---|---|---|
| unbiased | 0.918 / 0.595, w 0.628 | 0.897 / 0.008, w 0.421 | 0.913 / 0.000, w 0.122 | 0.914 / 0.000, w 0.105 | 0.915 / 0.000, w 0.087 |
| tumour_biased | 0.940 / 0.150, w 0.604 | 0.888 / 0.008, w 0.438 | 0.915 / 0.000, w 0.136 | 0.915 / 0.000, w 0.118 | 0.910 / 0.000, w 0.094 |

**Hurdle-aware uncertainty (amendment 1):**

| setting | ECE raw q → q+c → HAC q̂ | median Brier raw → HAC | magnitude interval coverage (expressing) / width: HAC vs raw relu(μ) |
|---|---|---|---|
| unbiased | 0.120 → 0.033 → **0.005** | 0.0273 → **0.0104** | **0.904 / 0.399** vs 0.896 / 0.772 |
| tumour_biased | 0.119 → 0.033 → **0.007** | 0.0272 → **0.0107** | **0.910 / 0.427** vs 0.905 / 0.803 |

## ATERA (confirmation)

**Marginal per-gene intervals (original rule):** median coverage (all entries / expressing entries), mean width

| setting | GHIST+ | GHIST+-Soft | RRC-authors | HAC | HAC-mondrian |
|---|---|---|---|---|---|
| unbiased | 0.912 / 0.751, w 0.702 | 0.901 / 0.010, w 0.523 | 0.910 / 0.000, w 0.364 | 0.906 / 0.000, w 0.343 | 0.907 / 0.045, w 0.318 |
| tumour_biased | 0.918 / 0.024, w 0.609 | 0.857 / 0.007, w 0.505 | 0.870 / 0.000, w 0.333 | 0.870 / 0.000, w 0.316 | 0.881 / 0.022, w 0.297 |

**Hurdle-aware uncertainty (amendment 1):**

| setting | ECE raw q → q+c → HAC q̂ | median Brier raw → HAC | magnitude interval coverage (expressing) / width: HAC vs raw relu(μ) |
|---|---|---|---|
| unbiased | 0.058 → 0.020 → **0.003** | 0.0599 → **0.0451** | **0.899 / 0.798** vs 0.901 / 0.993 |
| tumour_biased | 0.060 → 0.033 → **0.006** | 0.0634 → **0.0477** | **0.891 / 0.733** vs 0.886 / 0.930 |

## Pre-declared decisions

- Breast5K (dev), hurdle-aware rule: ECE 0.120→0.005, Brier 0.0273→0.0104, magnitude coverage 0.904 → **PASS**
- Breast5K (dev), original marginal rule: HAC-mondrian coverage 0.915, width 0.087 → PASS (but expressing-entry coverage 0.000: marginal intervals are not informative for expressing cells)
- ATERA (confirmation), hurdle-aware rule: ECE 0.058→0.003, Brier 0.0599→0.0451, magnitude coverage 0.899 → **PASS**
- ATERA (confirmation), original marginal rule: HAC-mondrian coverage 0.907, width 0.318 → PASS (but expressing-entry coverage 0.045: marginal intervals are not informative for expressing cells)
