# Where GHIST+-Adapt helps and where it hurts (exploratory)

Primary config (4 cores x 1 mm², diverse), mean of 3 repeats. Δ = Adapt − GHIST+, signed so that **positive = better**.

## Breast5K — strata by GHIST+ H&E-predicted type

| cell type | cells | profile PCC GHIST+ → Adapt | gene PCC within type | RMSE | RMSE (measured>0) | verdict |
|---|---|---|---|---|---|---|
| B | 92,291 | 0.148 → 0.254 | 0.015 → 0.024 | 0.223 → 0.072 | 0.701 → 0.717 (DM 0.706) | use Adapt |
| Epithelial | 19,708 | 0.190 → 0.246 | 0.034 → 0.044 | 0.288 → 0.106 | 0.705 → 0.719 (DM 0.695) | use Adapt |
| Macrophage | 51,645 | 0.102 → 0.194 | 0.027 → 0.039 | 0.301 → 0.081 | 0.702 → 0.725 (DM 0.704) | use Adapt |
| Malignant | 126,770 | 0.297 → 0.402 | 0.022 → 0.040 | 0.283 → 0.114 | 0.708 → 0.715 (DM 0.700) | use Adapt |
| Myeloid | 21,902 | 0.112 → 0.211 | 0.036 → 0.049 | 0.292 → 0.081 | 0.688 → 0.720 (DM 0.695) | use Adapt |
| Plasma | 45,922 | 0.100 → 0.190 | 0.024 → 0.037 | 0.306 → 0.080 | 0.693 → 0.720 (DM 0.701) | use Adapt |
| Stromal | 66,915 | 0.161 → 0.256 | 0.034 → 0.047 | 0.361 → 0.094 | 0.702 → 0.718 (DM 0.692) | use Adapt |
| T | 75,032 | 0.142 → 0.244 | 0.016 → 0.025 | 0.244 → 0.070 | 0.693 → 0.715 (DM 0.704) | use Adapt |
| Unassigned | 17,098 | 0.080 → 0.161 | 0.034 → 0.050 | 0.209 → 0.065 | 0.693 → 0.720 (DM 0.701) | use Adapt |

Most improved genes: C11orf96 (-0.00→0.32), CACNG4 (0.12→0.41), ECM1 (0.10→0.37), CTTN (0.31→0.52), STC2 (0.21→0.40)

Most worsened genes: MYH11 (0.33→0.31), SOX17 (0.16→0.15), SPIB (0.03→0.03), TAC1 (0.04→0.03), EPHA4 (0.04→0.04)

## ATERA — strata by GHIST+ H&E-predicted type

| cell type | cells | profile PCC GHIST+ → Adapt | gene PCC within type | RMSE | RMSE (measured>0) | verdict |
|---|---|---|---|---|---|---|
| B | 4,317 | 0.198 → 0.303 | 0.053 → 0.061 | 0.339 → 0.203 | 0.744 → 0.749 (DM 0.727) | use Adapt |
| Epithelial | 32,068 | 0.273 → 0.624 | 0.036 → 0.062 | 0.320 → 0.218 | 0.768 → 0.748 (DM 0.731) | use Adapt |
| Macrophage | 8,173 | 0.197 → 0.330 | 0.055 → 0.080 | 0.342 → 0.196 | 0.763 → 0.776 (DM 0.753) | use Adapt |
| Malignant | 41,652 | 0.407 → 0.768 | 0.017 → 0.039 | 0.333 → 0.220 | 0.752 → 0.733 (DM 0.723) | use Adapt |
| Myeloid | 2,498 | 0.180 → 0.290 | 0.062 → 0.074 | 0.323 → 0.197 | 0.749 → 0.770 (DM 0.748) | use Adapt |
| Plasma | 2,592 | 0.208 → 0.329 | 0.079 → 0.092 | 0.351 → 0.216 | 0.754 → 0.775 (DM 0.744) | use Adapt |
| Stromal | 26,550 | 0.275 → 0.416 | 0.048 → 0.056 | 0.342 → 0.190 | 0.742 → 0.756 (DM 0.736) | use Adapt |
| T | 4,066 | 0.197 → 0.303 | 0.045 → 0.060 | 0.328 → 0.197 | 0.745 → 0.746 (DM 0.726) | use Adapt |
| Unassigned | 1,730 | 0.263 → 0.480 | 0.074 → 0.115 | 0.311 → 0.201 | 0.765 → 0.769 (DM 0.740) | use Adapt |

Most improved genes: DST (0.03→0.52), C15orf48 (0.07→0.48), CD36 (-0.00→0.38), CTTN (0.25→0.63), RUNX1 (0.19→0.55)

Most worsened genes: MS4A1 (0.14→0.12), GPX3 (0.02→0.01), ADH1B (0.07→0.06), BAALC (0.04→0.03), CD83 (0.08→0.07)

## BreastCancer2 — strata by measured scClassify label

| cell type | cells | profile PCC GHIST+ → Adapt | gene PCC within type | RMSE | RMSE (measured>0) | verdict |
|---|---|---|---|---|---|---|
| B | 1,051 | 0.280 → 0.340 | 0.016 → 0.041 | 0.546 → 0.260 | 0.804 → 0.681 (DM 0.697) | use Adapt |
| Epithelial | 4,778 | 0.495 → 0.619 | 0.051 → 0.075 | 0.671 → 0.335 | 0.839 → 0.698 (DM 0.705) | use Adapt |
| Macrophage | 8,881 | 0.330 → 0.465 | 0.072 → 0.110 | 0.690 → 0.314 | 0.886 → 0.709 (DM 0.701) | use Adapt |
| Malignant | 7,502 | 0.509 → 0.700 | 0.086 → 0.147 | 0.769 → 0.361 | 0.854 → 0.697 (DM 0.706) | use Adapt |
| Myeloid | 1,383 | 0.233 → 0.369 | 0.033 → 0.064 | 0.727 → 0.297 | 0.856 → 0.686 (DM 0.683) | use Adapt |
| Plasma | 1,130 | 0.368 → 0.446 | 0.028 → 0.044 | 0.645 → 0.285 | 0.821 → 0.687 (DM 0.689) | use Adapt |
| Stromal | 5,461 | 0.399 → 0.509 | 0.081 → 0.119 | 0.748 → 0.309 | 0.885 → 0.700 (DM 0.686) | use Adapt |
| T | 23,018 | 0.317 → 0.469 | 0.025 → 0.053 | 0.551 → 0.259 | 0.825 → 0.694 (DM 0.707) | use Adapt |
| Unassigned | 4,314 | 0.406 → 0.491 | 0.069 → 0.091 | 0.732 → 0.325 | 0.874 → 0.694 (DM 0.684) | use Adapt |

Most improved genes: ALDH1A3 (-0.01→0.51), AKR1C1 (0.05→0.43), ESR1 (0.19→0.50), KRT23 (0.07→0.37), EGFR (0.15→0.45)

Most worsened genes: CD19 (0.05→0.05), TAC1 (0.09→0.08), TPSAB1 (0.06→0.06), FOXC2 (0.08→0.07), ADH1B (0.32→0.32)

## Marker-gene programmes (per-gene PCC change, whole evaluation set)

| slide | programme | genes ↑ / ↓ | mean ΔPCC | per gene |
|---|---|---|---|---|
| Breast5K | T | 3 / 0 | +0.019 | CD3E +0.017, CD3G +0.000, CD247 +0.019, TRAC +0.039 |
| Breast5K | B | 1 / 0 | +0.003 | BANK1 +0.000, CD19 +0.000, CD79A +0.000, CD79B +0.000, MS4A1 +0.013 |
| Breast5K | Plasma | 2 / 0 | +0.014 | MZB1 +0.007, PRDM1 +0.000, TENT5C +0.050, TNFRSF17 +0.000 |
| Breast5K | Macrophage | 3 / 0 | +0.039 | C1QC +0.111, CD68 +0.052, CD163 +0.000, MMP12 +0.000, MRC1 +0.030 |
| Breast5K | Myeloid | 3 / 0 | +0.012 | CD14 +0.001, CLEC9A +0.009, FCER1A +0.000, IL3RA +0.002, ITGAX +0.059, LILRA4 -0.001, MNDA +0.015 |
| Breast5K | Stromal | 4 / 1 | +0.036 | CCDC80 +0.060, DPT +0.000, LRRC15 +0.000, LUM +0.084, MYH11 -0.026, PDGFRA -0.000, PDGFRB +0.044, POSTN +0.122 |
| Breast5K | Epithelial | 5 / 0 | +0.084 | AGR3 +0.021, CDH1 +0.044, EPCAM +0.096, FOXA1 +0.099, GATA3 +0.159 |
| Breast5K | Malignant | 4 / 0 | +0.102 | APOBEC3B +0.004, CCND1 +0.194, ERBB2 +0.107, MKI67 +0.092, TOP2A +0.114 |
| ATERA | T | 3 / 0 | +0.039 | CD3E +0.034, CD3G +0.000, CD247 +0.009, TRAC +0.111 |
| ATERA | B | 2 / 1 | +0.011 | BANK1 +0.042, CD19 -0.003, CD79A +0.028, CD79B +0.000, MS4A1 -0.013 |
| ATERA | Plasma | 3 / 0 | +0.057 | MZB1 +0.066, PRDM1 +0.063, TENT5C +0.100, TNFRSF17 +0.000 |
| ATERA | Macrophage | 6 / 0 | +0.021 | APOC1 +0.003, C1QA +0.049, C1QC +0.022, CD68 +0.010, CD163 +0.014, MMP12 +0.008, MRC1 +0.038 |
| ATERA | Myeloid | 5 / 0 | +0.012 | CD14 +0.018, CLEC9A +0.039, FCER1A +0.008, IL3RA +0.009, ITGAX +0.000, LILRA4 +0.004, MNDA +0.015, S100A8 +0.000 |
| ATERA | Stromal | 6 / 0 | +0.026 | ACTA2 +0.004, CCDC80 +0.004, DPT +0.003, LRRC15 +0.001, LUM +0.048, MMP2 +0.045, MYH11 +0.006, PDGFRA +0.028, PDGFRB +0.012, POSTN +0.105 |
| ATERA | Epithelial | 7 / 0 | +0.099 | AGR3 +0.319, CDH1 +0.360, EPCAM +0.056, FOXA1 +0.082, GATA3 +0.058, KRT5 +0.084, KRT7 +0.000, KRT8 +0.023, KRT14 +0.000, PIGR +0.004 |
| ATERA | Malignant | 6 / 0 | +0.042 | APOBEC3B +0.070, CCND1 +0.084, ERBB2 +0.000, LYPD3 +0.008, MKI67 +0.099, S100A14 +0.000, TACSTD2 +0.010, TOP2A +0.065 |
| BreastCancer2 | T | 4 / 0 | +0.091 | CD3E +0.073, CD3G +0.113, CD247 +0.115, TRAC +0.061 |
| BreastCancer2 | B | 2 / 0 | +0.017 | BANK1 +0.000, CD19 -0.004, CD79A +0.065, CD79B +0.003, MS4A1 +0.022 |
| BreastCancer2 | Plasma | 3 / 0 | +0.045 | MZB1 +0.056, PRDM1 +0.061, TENT5C +0.061, TNFRSF17 +0.000 |
| BreastCancer2 | Macrophage | 7 / 0 | +0.048 | APOC1 +0.027, C1QA +0.046, C1QC +0.047, CD68 +0.046, CD163 +0.078, MMP12 +0.010, MRC1 +0.079 |
| BreastCancer2 | Myeloid | 6 / 0 | +0.050 | CD14 +0.048, CLEC9A +0.000, FCER1A +0.000, IL3RA +0.046, ITGAX +0.032, LILRA4 +0.070, MNDA +0.018, S100A8 +0.189 |
| BreastCancer2 | Stromal | 9 / 0 | +0.056 | ACTA2 +0.134, CCDC80 +0.044, DPT +0.011, LRRC15 +0.000, LUM +0.047, MMP2 +0.061, MYH11 +0.028, PDGFRA +0.049, PDGFRB +0.046, POSTN +0.138 |
| BreastCancer2 | Epithelial | 10 / 0 | +0.104 | AGR3 +0.144, CDH1 +0.104, EPCAM +0.115, FOXA1 +0.099, GATA3 +0.136, KRT5 +0.094, KRT7 +0.129, KRT8 +0.124, KRT14 +0.025, PIGR +0.071 |
| BreastCancer2 | Malignant | 7 / 0 | +0.074 | APOBEC3B +0.138, CCND1 +0.058, ERBB2 +0.137, LYPD3 +0.158, MKI67 +0.019, S100A14 +0.000, TACSTD2 +0.073, TOP2A +0.011 |
