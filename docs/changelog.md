# Changelog

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
