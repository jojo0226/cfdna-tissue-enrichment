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

| Panel | Layers | Use |
|---|---|---|
| `prenatal` | ENCODE tissues (C) + maternal–fetal interface snATAC (W) | pregnancy outcomes, placental signal |
| `tumor` | ENCODE tissues (C) + TCGA tumour ATAC (T) | any cancer; set `--primary` for your tumour type |
| `glioma` | C + T, glioma / CNS tracks pre-registered | GBM / LGG projects |
| `pan` | C + T + W | exploration, nothing pre-registered |

Start with [Install and quick start](quickstart.md). The method is described in [Method](method.md).

Developed in the EpiFluid Lab (Northwestern University), from `cfdna-chromatin-xai`.
