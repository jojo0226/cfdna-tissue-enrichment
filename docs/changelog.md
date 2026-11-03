# Changelog

## 0.1.0 (2026-11-03)
- First standalone release, generalised from `pregnancy_tissue_kit` v1.1.
- Layers C (ENCODE DNase), T (TCGA ATAC), W (maternal–fetal interface snATAC) on hg38, 50 kb and 10 kb.
- Panels `prenatal`, `tumor`, `glioma`, `pan`.
- Modes `sets`, `sample`, `differential`; factorial designs with covariates and stratified permutations.
- Random bin sets matched on GC × mappability; top-bin p-values from label permutations.
