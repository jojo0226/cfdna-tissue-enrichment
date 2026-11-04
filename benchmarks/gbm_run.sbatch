#!/bin/bash
#SBATCH -A b1198
#SBATCH -p b1198
#SBATCH -t 24:00:00
#SBATCH --mem 48G
#SBATCH -c 4
#SBATCH -J cfte_gbm
#SBATCH -o /projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/logs/gbm_%j.log
# GBM benchmarks (hg38) + PRAD differential (hg19). Counts are cached; all modes.
set -eu
G=/projects/b1198/epifluidlab/jocelyn/software/cfdna-tissue-enrichment; B=/projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks; PY=/home/pfu2606/.conda/envs/meth/bin/python
O=$B/gbm; cd $G
$PY benchmarks/gbm_prepare.py $O
$PY -m cfdna_tissue_enrichment count --samples $O/gbm58_samples.csv --cache $O/cache
$PY -m cfdna_tissue_enrichment count --samples $O/prog_samples.csv --cache $O/cache
IMM="C_neutrophil;C_CD14-positive monocyte;C_B cell;C_natural killer cell;C_CD4-positive, alpha-beta T cell;C_CD8-positive, alpha-beta T cell"
run() { $PY -m cfdna_tissue_enrichment run "$@" 2>&1 | grep -v Warning | tail -2; }
run --panel glioma --samples $O/sheet_GBM_vs_schwannoma.csv --cache $O/cache --out $O/runs/GBM_vs_schwannoma --group diagnosis --reference-level Schwannoma --case GBM --covariates sex:cat,age
run --panel tumor --samples $O/sheet_Breast_met_vs_schwannoma.csv --cache $O/cache --out $O/runs/Breast_met_vs_schwannoma --group diagnosis --reference-level Schwannoma --case Breast_met --covariates sex:cat,age --primary "T_BRCA;$IMM"
run --panel glioma --samples $O/prog_samples.csv --cache $O/cache --out $O/runs/progression_vs_pseudo --group group --reference-level Pseudoprogression --case Progression
# PRAD mixtures: differential mode, spiked >= 5% vs unspiked, per depth
P=$B/prad_mix
for dep in 10 30; do
  $PY -c "import pandas as pd; d=pd.read_csv('$P/samples_mixtures.csv'); d=d[(d.depth_M==$dep)&((d.tf_pct>=5)|(d.tf_pct==0))]; d['hi']=(d.tf_pct>=5).map({True:'yes',False:'no'}); d.to_csv('$P/samples_diff_d${dep}M.csv', index=False)"
  run --panel tumor --build hg19 --samples $P/samples_diff_d${dep}M.csv --cache $P/cache --mapq-col 0 --out $P/runs/diff_d${dep}M --modes differential --group hi --reference-level no --case yes --primary T_PRAD
done
$PY $G/benchmarks/gbm_summarise.py $B
