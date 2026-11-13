# Outputs and interpretation

## Files (in `--out`)

| file | content |
|---|---|
| `qc.csv`, `qc.png` | fragments read / used / in grid, median length, short fraction (100–150 bp), zero bins, low-depth flag |
| `design_matrix.csv` | the model: intercept, factor indicators, covariates |
| `sample/scores_E.csv` | samples × tracks: E vs the reference group |
| `sample/effect_tests.csv` | per track × effect: beta, t, `p_ols`, `p_perm`, `q_perm`, `q_perm_primary` |
| `sample/pairwise_vs_reference.csv` | each group vs the reference group: n, difference, AUC, `p_mannwhitney`, `p_perm`, q values |
| `*/scores_heatmap.png`, `*/primary_by_group.png` | per-sample heatmap; primary tracks by group |
| `*/ranked_<effect>.png` | every track ranked, two panels (lower / higher coverage in cases), coloured by type (immune/blood, tissue, TCGA, derived); sample: group-effect t on E, differential: all-bins rho; * permutation p < 0.05, ** and q < 0.10. Redraw with `cfdna-tissue-enrich plot --out RUN [--top N]` |
| `differential/enrichment_all_bins.csv` (+ `.png`) | per effect × track: Spearman rho, `p_perm`, `q_perm`, `q_perm_primary` |
| `differential/enrichment_top_bins.csv` | per effect × direction × track: mean specificity of top bins, `p_perm`, q values; `z_vs_matched_random` with GC / mappability balance and `z_closed_form` (both descriptive) |
| `differential/bin_effects.csv.gz` | per effect × bin: moderated t and p (`--diff-stat pooled`: pooled difference `d`) |
| `differential/top_bins_activity_<effect>.png` | top-bin tracks by `z_closed_form`, one panel per coverage direction, titled with the cfDNA reading (more coverage = less active) |
| `differential/enrichment_top_bins_loo.csv` | `--loo`: per dropped sample × direction × track, mean specificity and `z_closed_form` |
| `run_manifest.json` | version, panel definition, layer provenance, design, every parameter |

## Reading the results

- **Scores are relative.** They compare coverage at a track's specific bins with other bins in the same samples. Compare groups; do not read a value as a fraction of cfDNA.
- **Check the direction with a known contrast first.** Whether more DNA from a cell type raises or lowers coverage at its specific bins depends on assay and depth. Examples: pregnant vs non-pregnant (placental tracks), high vs low tumour fraction (tumour tracks), early vs late gestation.
- **Agreement across modes is the strongest evidence.** `sample` is per-sample; `differential` is group level. A track moving the same way in both is more convincing than one significant in a single mode.
- **Small groups.** With < 8 per group, report effect sizes and AUC rather than p-values alone. Differential mode needs residual df ≥ 8; the log warns when lower.
- **Covariates.** Adjust for what moves the tissue mix: gestational age in pregnancy; age, sex and batch in most cohorts.

## Limits

- hg38 only so far (a T2T build is planned).
- The C and T layers come from an hg19 atlas moved to hg38 by bin midpoint (10 kb); W was placed on hg38 natively.
- W windows are not peaks; at 50 kb their signal is diluted, which is why everything is based on specificity.
- Seven W cell types have mixed fetal / maternal origin (dNK, T, B, FB, pvSMC, cDC, Ery); they are kept as tracks but excluded from the origin axis.
- Differential mode's moderated t uses a fixed prior (4 df, median variance), a simple stand-in for limma.
