#!/bin/bash
#SBATCH -A b1198
#SBATCH -p b1198
#SBATCH -t 24:00:00
#SBATCH --mem 48G
#SBATCH -c 4
#SBATCH -J cfte_cabat_run
#SBATCH -o /projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/logs/cabat_analyse_%j.log
# Per cancer type: cases vs all controls, tumor panel, all three modes; then controls-only nulls (sample).
set -eu
G=/projects/b1198/epifluidlab/jocelyn/software/cfdna-tissue-enrichment; O=/projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/cabat; PY=/home/pfu2606/.conda/envs/meth/bin/python
declare -A TR=([Prostate]=T_PRAD [Colorectal]=T_COAD [Liver]=T_LIHC,C_liver [Lung]=T_LUAD,T_LUSC,Tgroup_lung [Breast]=T_BRCA [Stomach]=T_STAD [Pancreatic]=C_pancreas [Ovary]=C_ovary)
cd $G
for ct in Prostate Colorectal Liver Lung Breast Stomach Pancreatic Ovary; do
  echo "== $ct $(date +%T)"
  $PY -m cfdna_tissue_enrichment run --panel tumor --samples $O/sheet_$ct.csv --cache $O/cache --out $O/runs/$ct \
      --group status --reference-level control --case case --covariates batch:cat,sex:cat --primary ${TR[$ct]},C_neutrophil 2>&1 | grep -v Warning | tail -3
done
for k in 0 1 2 3 4 5 6 7 8 9; do
  $PY -m cfdna_tissue_enrichment run --panel tumor --samples $O/null_$k.csv --cache $O/cache --out $O/null/run_$k --modes sample \
      --group fake --reference-level A --case B --covariates batch:cat,sex:cat 2>&1 | grep -v Warning | tail -1
done
for k in 0 1 2; do
  $PY -m cfdna_tissue_enrichment run --panel tumor --samples $O/null_$k.csv --cache $O/cache --out $O/null/diff_$k --modes differential \
      --group fake --reference-level A --case B --covariates batch:cat,sex:cat 2>&1 | grep -v Warning | tail -1
done
$PY $G/benchmarks/cabat_summarise.py $O
