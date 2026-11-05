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
| `frag_path` | yes | fragment file; absolute or relative to the sheet |
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

cfdna-tissue-enrich panels                     # list panels
cfdna-tissue-enrich tracks --panel prenatal    # list a panel's tracks; * = primary
```

## Main options

| option | default | meaning |
|---|---|---|
| `--panel` | `pan` | panel name or a panel JSON path |
| `--modes` | `sample,differential` | either or both; both need groups (`--factor` or `--group`) |
| `--factor` | — | `COLUMN:REFERENCE:CASE`, repeatable |
| `--group`, `--reference-level`, `--case` | — | shortcut for one factor |
| `--interaction` | off | factor × factor terms |
| `--covariates` | — | comma-separated; numeric with ≥ 6 values = linear, else categorical; force with `:num` / `:cat` |
| `--primary` | panel's list | pre-registered tracks (own BH correction) |
| `--bin` | `50kb` | `50kb` or `10kb` |
| `--top-n` | 2000 | differential mode, bins per direction |
| `--n-perm`, `--n-perm-bins` | 2000, 200 | permutations for score tests / differential mode |
| `--seed` | 17 | same inputs + seed = identical output |
