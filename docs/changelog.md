# Changelog

## 0.4.0 (2026-11-19)
- **`markers` mode** (third mode): 10 kb marker-set coverage (mean over each set's top-N bins of `--markers`) with GC-decile
  matched random sets, as in the C1D1 and Ahmed NU21C01 analyses (cfdna_chromatin.enrichment `markers`).
- `sample` mode: `--sample-reference all`, and the per-sample report (paired sign-flip or unpaired AUC contrasts, PCA,
  heatmaps, strip plots; `--group-col`, `--groups`, `--contrasts`, `--patient-col`, `--label-col`, `--connect`, `--pca-centre`),
  also used by `markers`.
- `blood_cns_glioma` panel and `--specificity combined-z`; panels can set `specificity` and a `tracks` subset.
  `sample --panel blood_cns_glioma --sample-reference all` is the C1D1 / Ahmed 50 kb tissue enrichment.
- Samples can give a precomputed per-sample coverage csv (`coverage` column) instead of `frag_path`.
- Regression (`tests/regression_c1d1_ahmed.py`, both cohorts vs the original outputs): per-sample tissue scores within 1e-4
  (r = 1.000), marker scores identical, contrast effects / p / q identical (one Ahmed sign-flip p differs by one flip of 512).
  Random-set `random_p` differs by Monte Carlo noise only (r 0.96; the grid lists bins in numeric chromosome order, so the
  same seed draws different random sets); 7 of 8 and 2 of 3–4 sets called at random_p < 0.05 agree (Ahmed, C1D1).

## 0.3.0 (2026-11-13)
- `--gc-correct yes|no` (default yes): GC LOWESS can be switched off, for every mode.
- `--detrend-bins N`: rolling-median detrend window (default unchanged, ~2 Mb).
- `--diff-stat pooled` (default `t`): differential mode on pooled groups, d = detrend(log2 mean CPM cases + 1) - detrend(log2 mean
  CPM reference + 1), one two-level factor, no covariates; label-permutation nulls re-pool and re-select. Reproduces the
  cfdna-chromatin-xai GBM progression / pseudoprogression figure with `--gc-correct no --detrend-bins 101 --specificity combined`.
- `--specificity combined` (default `layer`): one specificity table across the panel's layers (openness minus the row mean), as in
  cfdna-chromatin-xai.
- `--loo`: leave-one-out re-pool + re-score of the top bins (pooled statistic) -> `enrichment_top_bins_loo.csv`.
- Top-bins table gains `z_closed_form` (closed-form random-bin z, descriptive); new figure `top_bins_activity_<effect>.png` labels
  each panel by coverage direction and the cfDNA reading (more coverage = less active).

## 0.2.0 (2026-11-05)
- **`sets` mode removed.** Two modes remain, `sample` (per sample) and `differential` (group level); both use every bin
  and agree closely. Why: [method, section 5](method.md#5-why-there-is-no-sets-mode). `--set-size` removed; a run
  without groups now stops with a message instead of falling back to `sets`.
- Benchmark scripts and write-up updated to the two modes.

## 0.1.2 (2026-11-05)
- Ranked bar plots of every track for each effect in all three modes (two panels like the cfdna-chromatin-xai plot, coloured by track type), `--bar-top`, and a `plot` command to redraw them for a finished run.
- Benchmark runs re-plotted with the new figures.

## 0.1.1 (2026-11-03)
- hg19 reference (layers C and T; clean grid from the hg19 ENCODE blacklist v2, UCSC gaps, N > 10%, CRG 50-mer mappability < 0.80).
- `count` command and `--cache`: per-sample counts computed once (parallel) and reused.
- MAPQ column detected from a fragment-file header (FinaleDB files: column 5).
- Benchmarks: PRAD mixtures (hg19), Cabat pan-cancer (hg38), null calibration.

## 0.1.0 (2026-11-03)
- First standalone release, generalised from `pregnancy_tissue_kit` v1.1.
- Layers C (ENCODE DNase), T (TCGA ATAC), W (maternal–fetal interface snATAC) on hg38, 50 kb and 10 kb.
- Panels `prenatal`, `tumor`, `glioma`, `pan`.
- Modes `sets`, `sample`, `differential`; factorial designs with covariates and stratified permutations.
- Random bin sets matched on GC × mappability; top-bin p-values from label permutations.
