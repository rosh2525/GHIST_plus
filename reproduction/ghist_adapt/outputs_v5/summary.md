# P2: support-aware trust gate (config_v5.json, incl. amendments 1–2)

Median per-gene PCC (all genes) on evaluation cells, mean of 3 repeats. OUT = cells outside the cores' support (λ < 1).

## Breast5K (dev)

| setting | out-of-support cells | OUT: Soft / RRC / HAC / HAC-trust | CI(trust−HAC) > 0 on OUT | ALL: HAC / HAC-trust | HAC worse than Soft on OUT? |
|---|---|---|---|---|---|
| unbiased | 7% | 0.0960 / 0.0970 / 0.0974 / 0.0999 | 1/3 | 0.1097 / 0.0989 | no |
| tumour_biased | 6% | 0.0950 / 0.0925 / 0.0886 / 0.0946 | 2/3 | 0.0971 / 0.0928 | yes |
| stroma_biased | 9% | 0.0951 / 0.0968 / 0.0917 / 0.0973 | 2/3 | 0.1020 / 0.0937 | yes |

## ATERA (confirmation)

| setting | out-of-support cells | OUT: Soft / RRC / HAC / HAC-trust | CI(trust−HAC) > 0 on OUT | ALL: HAC / HAC-trust | HAC worse than Soft on OUT? |
|---|---|---|---|---|---|
| unbiased | 4% | 0.1689 / 0.1714 / 0.1981 / 0.1922 | 0/3 | 0.1880 / 0.1870 | no |
| tumour_biased | 10% | 0.1712 / 0.1503 / 0.1904 / 0.1817 | 0/3 | 0.1881 / 0.1873 | no |
| stroma_biased | 2% | 0.1596 / 0.1551 / 0.1831 / 0.1834 | 0/3 | 0.1857 / 0.1862 | no |

## Pre-declared decision

- Breast5K (dev): problem exists in a biased setting = True; gate fixes it where it exists (CI>0 in ≥2/3) = True; unbiased all-cell cost +0.0109 (allowed ≤ 0.002) → FAIL
- ATERA (confirmation): problem exists in a biased setting = False; gate fixes it where it exists (CI>0 in ≥2/3) = n/a; unbiased all-cell cost +0.0010 (allowed ≤ 0.002) → FAIL

**Overall: P2 does not meet the bar on these slides.**

