# Benchmarks (v0.1.0, 2026-11-03)

Two benchmarks with known answers, both run with the `tumor` panel (layers C + T). Scripts are in `benchmarks/`;
result tables in `benchmarks/results/`.

## 1. PRAD in-silico mixtures (hg19): does the method find tumour DNA, and how little?

**Design.** Fragments from a prostate cancer plasma sample (PRAD_14230, ichorCNA tumour fraction 46.5%) mixed into
a healthy person's plasma (HD45) at 0, 0.5, 1, 2, 5, 10 and 20% tumour fraction; 10 M and 30 M fragments; 3
replicates each. A second healthy person (HD46) at 0% shows person-to-person variation. Detection = all 3 replicates
beyond every unspiked sample (HD45 and HD46) at that depth. All files hg19 (native build of C and T).

![PRAD dose response](img/prad_dose_response.png)

| Depth | Mode | Spearman, T_PRAD vs tumour fraction | Rank of T_PRAD among 23 TCGA tracks | Detection limit |
|---|---|---|---|---|
| 10 M | sets | −0.63 | 1 | 20% |
| 10 M | sample | −0.88 | 1 | **2%** |
| 30 M | sets | −0.77 | 1 | 5% |
| 30 M | sample | −0.94 | 1 | **2%** |

- **T_PRAD is the most tumour-fraction-responsive TCGA track in every setting.**
- **Direction: more tumour DNA lowers coverage at the tumour's specific open regions** (negative score). Read
  tumour tracks this way: lower = more tumour-derived signal.
- `sample` mode (all bins) is more sensitive than `sets` mode (top-1,000 bins), especially at low depth.
- Caveat: the null has only two healthy people. The earlier benchmark with four people showed that person-to-person
  variation, not depth, sets the floor, so real-world limits will be higher than 2%.

## 2. Cabat pan-cancer cohort (hg38): does each cancer light up its own track?

**Design.** 139 cases (pancreatic 28, lung 23, breast 22, prostate 21, colorectal 18, ovary 12, liver 8, stomach 7)
and 139 matched controls, QC-pass. Each cancer type vs all controls, adjusted for batch and sex; own TCGA track as
primary (pancreas / ovary have no TCGA track, so the ENCODE tissue was used). ichorCNA tumour fractions are mostly
below ~3%.

![Cabat matrix](img/benchmark_matrix.png)

| Cancer | Own track | sets t (p) | sample t (p) | differential rho (p) | Spearman with ichorCNA TF (sets) |
|---|---|---|---|---|---|
| Breast | T_BRCA | 1.79 (0.07) | **3.18 (0.002), rank 1/23** | **0.023 (<0.001), rank 1/23** | −0.29 |
| Stomach | T_STAD | 1.27 (0.21) | **2.61 (0.012), rank 3/23** | 0.017 (0.03) | **−0.88 (p 0.009, n 7)** |
| Liver | T_LIHC | −1.09 (0.28) | −2.49 (0.017) | −0.043 (0.01) | 0.26 |
| Prostate | T_PRAD | 1.38 (0.16) | 1.07 (0.28) | 0.013 (0.23) | −0.22 |
| Colorectal | T_COAD | 1.77 (0.08) | 0.21 (0.83) | 0.001 (0.84) | −0.34 |
| Lung | T_LUSC / T_LUAD | 2.02 (0.04) / 0.40 | −1.23 / −1.13 | n.s. | 0.09 / **−0.42 (p 0.04)** |
| Pancreatic | C_pancreas | −2.30 (0.02) | −0.83 | n.s. | −0.08 |
| Ovary | C_ovary | −2.18 (0.03) | −3.93 (<0.001) | −0.086 (<0.001) | −0.56 (p 0.06) |

- **No diagonal.** Most cancer types share the same pattern of TCGA tracks (see the matrix rows), so the main
  case-vs-control difference in this cohort is common to all cancers, not tumour-type specific. At ichorCNA
  fractions of ~1–3%, tumour-specific DNA is near the detection limit found in benchmark 1.
- **Breast** is the one clear own-track hit (rank 1 of 23 in sample and differential modes).
- **Within cases, the own-track score falls as ichorCNA tumour fraction rises** in 6 of 8 cancer types, the same
  direction as the mixture benchmark (stomach and lung LUAD individually significant). This is the expected
  tumour signal, but weak at these fractions.
- The shared pattern needs follow-up before biological interpretation: cases and controls may differ in collection
  era or processing, which no covariate here captures.

## 3. Null calibration (Cabat controls, random labels)

| Mode | Runs | Mean share of tracks at p < 0.05 | Max in one run |
|---|---|---|---|
| sets | 10 | 5.2% | 24.8% |
| sample | 10 | 5.0% | 16.8% |
| differential | 3 | 0.3% | 1.0% |

On average the per-sample modes are calibrated (~5%). Single runs can show many tracks at once because tracks are
correlated (one chance shift moves related tissues together); use the q-values and the primary-track correction.
Differential mode is conservative.

## Take-aways

1. The method recovers tumour-type-specific signal when tumour DNA is ≥ 2–5% (`sample` mode most sensitive).
2. Tumour tracks go **down** with more tumour DNA.
3. In low-fraction real cohorts, expect shared, non-specific case/control shifts; pre-register tracks and check
   dose-response with an independent tumour-fraction estimate (ichorCNA) before interpreting.
