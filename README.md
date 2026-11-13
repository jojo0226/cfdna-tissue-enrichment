# cfdna-tissue-enrichment

Which tissues and cell types does a group difference in plasma cfDNA coverage come from? This tool answers
that with **one standardised method** against interchangeable **chromatin-accessibility atlases**:
ENCODE tissues, TCGA tumours, the maternal–fetal interface, and any atlas you add.

```
cfDNA fragment files ──► clean hg38 bins ──► GC-corrected, detrended coverage
                                                   │
          atlas layers (C, T, W, …) ──► track specificity per bin
                                                   │
                              ┌────────────────────┴────────────────────┐
                            sample                                differential
                  per sample, all bins vs a              per bin group effect, then
                  reference group                         tissue enrichment of it
```

- **Same method for every project.** Only the *panel* changes: which atlas layers and which pre-registered tracks.
- **Two readouts.** A per-sample score (`sample`) and a group-level enrichment (`differential`). Both use every bin, with no bin selection, and agree closely (track rankings r ≈ 0.91–0.97 across the benchmarks).
- **Group designs** from a single comparison to factorial (e.g. term/preterm × SGA/non-SGA) with covariates.
- **Clean by construction.** ENCODE blacklist, gaps, high-N and low-mappability bins are removed before anything is counted.
- **Comparable with cfdna-chromatin-xai (0.3.0).** Options switch off GC correction, change the detrend window, pool groups
  instead of fitting a per-bin t, and use one specificity table across layers; with these settings the tool reproduces the
  cfdna-chromatin-xai GBM progression / pseudoprogression figure (r = 0.98), while still reporting label-permutation p-values.

| Panel | Layers | Use |
|---|---|---|
| `prenatal` | ENCODE tissues (C) + maternal–fetal interface snATAC (W) | pregnancy outcomes, placental signal |
| `tumor` | ENCODE tissues (C) + TCGA tumour ATAC (T) | any cancer; set `--primary` for your tumour type |
| `glioma` | C + T, glioma / CNS tracks pre-registered | GBM / LGG projects |
| `pan` | C + T + W | exploration, nothing pre-registered |

## Quick start

```bash
conda env create -f environment.yml && conda activate cfdna-tissue-enrichment && pip install -e .
bash tests/smoke_test.sh
cfdna-tissue-enrich run --panel prenatal --samples samples.csv --out results/ \
    --factor birth:term:preterm --factor sga:nonSGA:SGA --covariates gestational_age_weeks
```

### Reproducing a cfdna-chromatin-xai figure (0.3.0)

```bash
cfdna-tissue-enrich run --panel glioma --samples samples.csv --out results_xai/ --modes differential \
    --group group --reference-level pseudo --case progression \
    --gc-correct no --detrend-bins 101 --diff-stat pooled --specificity combined --loo
```

| option | default | effect |
|---|---|---|
| `--gc-correct yes\|no` | `yes` | per-sample GC LOWESS (all modes) |
| `--detrend-bins N` | 41 (50 kb) / 205 (10 kb) | rolling-median detrend window, ~2 Mb by default |
| `--diff-stat t\|pooled` | `t` | differential per-bin statistic: moderated t, or the difference of pooled groups (one two-level factor, no covariates) |
| `--specificity layer\|combined` | `layer` | specificity within each layer, or across all panel layers (cfdna-chromatin-xai) |
| `--loo` | off | pooled statistic: leave-one-out re-pool + re-score of the top bins |

Defaults reproduce 0.2.0 exactly. `z_closed_form` (the cfdna-chromatin-xai `spec_z`) is reported for top bins but is
descriptive: top bins are selected from the data, so judge significance by `p_perm`. In the GBM benchmark the pattern
reproduced at z up to ±29 had `p_perm` ≥ 0.11 and largely disappeared with GC correction. Direction: in cfDNA, open chromatin
loses coverage, so top bins with **more coverage in cases** point to cell types **more active in the reference group**; the
`top_bins_activity` figure labels both.

Documentation: `docs/` (published as a GitHub Pages site by `.github/workflows/docs.yml`; preview locally with
`pip install mkdocs-material && mkdocs serve`).

| Page | |
|---|---|
| [Install and quick start](docs/quickstart.md) | inputs, commands, options |
| [Method](docs/method.md) | the standardised method |
| [Panels and layers](docs/panels.md) | what each atlas and panel contains |
| [Outputs and interpretation](docs/outputs.md) | files and how to read them |
| [Benchmarks](docs/benchmarks.md) | PRAD mixtures, pan-cancer cohort, null calibration |
| [Add a new atlas](docs/add_atlas.md) | how to plug in another reference |

## Repository layout

```
cfdna_tissue_enrichment/   package (cli, reference, counts, stats, modes, figures)
panels/                    panel definitions (JSON)
reference/hg38/            grid + layers C, T, W at 50 kb and 10 kb, layers.json, tracks_W.csv
reference/hg19/            grid + layers C, T at 50 kb and 10 kb (native hg19 build of the atlas)
benchmarks/                benchmark scripts and result tables
examples/                  8 small fragment files (BH01 subsamples; arbitrary labels) + samples.csv
tests/smoke_test.sh        runs every panel on the examples
build/build_layers.py      how the reference layers were built (maintainers)
docs/, mkdocs.yml          documentation site
```

## Citation

Cite the atlas sources you use: ENCODE Project Consortium (C); Corces et al., *Science* 2018 (T);
Wang et al., *Nature* 653:167–179 (2026), with pseudobulk matrices prepared by Paz Polak (W).
Example data: BH01 (Snyder et al., *Cell* 2016).

License: to be chosen by the lab before the repository is made public.

Developed in the EpiFluid Lab (Northwestern University), from `cfdna-chromatin-xai`.
