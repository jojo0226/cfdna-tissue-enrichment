# Install and quick start

## Install

```bash
git clone <repo-url> cfdna-tissue-enrichment && cd cfdna-tissue-enrichment
conda env create -f environment.yml && conda activate cfdna-tissue-enrichment
pip install -e .                      # provides the `cfdna-tissue-enrich` command
```

Python ≥ 3.9 with numpy, pandas, scipy, statsmodels and matplotlib. No internet or cluster paths are needed at run time.
The hg38 reference (`reference/hg38/`, ~80 MB) ships with the repository.

## Smoke test

```bash
bash tests/smoke_test.sh              # every panel on the example data + a reproducibility check
```

## Inputs

**Sample sheet (CSV)**, one row per sample:

| column | required | meaning |
|---|---|---|
| `sample_id` | yes | unique name |
| `frag_path` | one of these two | fragment file; absolute or relative to the sheet |
| `coverage` | | per-sample coverage csv (`key` + one `*_cpm` column, e.g. from `fragbed6_to_coverage_mapq.py`); CPM used as given |
| group columns | for group analyses | one column per two-level factor, e.g. `birth` = term / preterm, `sga` = nonSGA / SGA |
| covariates | no | e.g. `gestational_age_weeks`, `age`, `sex`, `batch` |

**Fragment files**: one fragment per line, tab-separated, plain or `.gz`, **aligned to hg38**:
`chrom  start  end  [MAPQ]  [...]`. MAPQ is read from column 4 when numeric (`--mapq-col N` to change, `0` for none).
FinaleDB fragment files work directly. Recommended depth: ≥ 5 M fragments at 50 kb, ≥ 20 M at 10 kb.

## Run

```bash
# pregnancy, 2 x 2 design, reference group = term non-SGA
cfdna-tissue-enrich run --panel prenatal --samples samples.csv --out results/ \
    --factor birth:term:preterm --factor sga:nonSGA:SGA --covariates gestational_age_weeks,bmi,batch:cat

# glioma vs controls
cfdna-tissue-enrich run --panel glioma --samples samples.csv --out results/ \
    --group diagnosis --reference-level control --case GBM --covariates age,sex

# another cancer: tumor panel with your own pre-registered tracks
cfdna-tissue-enrich run --panel tumor --samples samples.csv --out results/ \
    --group diagnosis --reference-level control --case PRAD --primary T_PRAD,C_neutrophil

# per-sample tissue enrichment and marker coverage, paired timepoints (the Ahmed NU21C01 analysis)
cfdna-tissue-enrich run --panel blood_cns_glioma --modes sample --sample-reference all --samples s50.csv --out tissue_50kb/ \
    --groups Pre,Surgery,On-treatment --contrasts Surgery:Pre,On-treatment:Surgery,On-treatment:Pre --patient-col patient --label-col cycle --connect
cfdna-tissue-enrich run --modes markers --bin 10kb --markers immune_marker_bins_10kb_immune_GBM_LGG.csv.gz --samples s10.csv --out markers_10kb/ \
    --groups Pre,Surgery,On-treatment --contrasts Surgery:Pre,On-treatment:Surgery,On-treatment:Pre --patient-col patient --label-col cycle --connect

cfdna-tissue-enrich panels                     # list panels
cfdna-tissue-enrich tracks --panel prenatal    # list a panel's tracks; * = primary
```

## Main options

| option | default | meaning |
|---|---|---|
| `--panel` | `pan` | panel name or a panel JSON path |
| `--modes` | `sample,differential` | any of `sample`, `differential`, `markers`; differential needs groups (`--factor` or `--group`), sample too unless `--sample-reference all` |
| `--sample-reference` | `group` | `all`: sample-mode z against every sample (C1D1 / Ahmed) (0.4.0) |
| `--markers`, `--marker-top` | —, 1000 | markers mode: marker atlas and bins per set; use `--bin 10kb` (0.4.0) |
| `--group-col`, `--groups`, `--contrasts` | `group`, all, adjacent | per-sample report: group column, groups kept (in order), `case:ref` contrasts (0.4.0) |
| `--patient-col`, `--label-col`, `--connect` | `patient` | paired sign-flip tests when patients repeat (`none` = unpaired); plot labels; PCA arrows (0.4.0) |
| `--pca-centre` | sample `no`, markers `both` | report PCA on raw and / or sample-centred scores (0.4.0) |
| `--factor` | — | `COLUMN:REFERENCE:CASE`, repeatable |
| `--group`, `--reference-level`, `--case` | — | shortcut for one factor |
| `--interaction` | off | factor × factor terms |
| `--covariates` | — | comma-separated; numeric with ≥ 6 values = linear, else categorical; force with `:num` / `:cat` |
| `--primary` | panel's list | pre-registered tracks (own BH correction) |
| `--bin` | `50kb` | `50kb` or `10kb` |
| `--top-n` | 2000 | differential mode, bins per direction |
| `--n-perm`, `--n-perm-bins` | 2000, 200 | permutations for score tests / differential mode |
| `--seed` | 17 | same inputs + seed = identical output |
| `--gc-correct` | `yes` | `no` skips the per-sample GC LOWESS (0.3.0) |
| `--detrend-bins` | 41 / 205 | rolling-median detrend window in bins (0.3.0) |
| `--diff-stat` | `t` | `pooled` = difference of pooled groups, one two-level factor, no covariates (0.3.0) |
| `--specificity` | panel's, else `layer` | `combined` = one table across the panel's layers, as cfdna-chromatin-xai (0.3.0); `combined-z` = per-bin z across all panel tracks, as C1D1 / Ahmed (0.4.0) |
| `--loo` | off | pooled statistic: leave-one-out re-pool + re-score of the top bins (0.3.0) |
