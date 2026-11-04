#!/bin/bash
#SBATCH -A b1198
#SBATCH -p b1198
#SBATCH -t 12:00:00
#SBATCH --mem 96G
#SBATCH -c 2
#SBATCH -J cfte_prad
#SBATCH -o /projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/logs/prad_%j.log
# PRAD in-silico mixtures on hg19: build mixtures into the count cache, then run the tumor panel (sets + sample;
# reference group = unspiked samples, both healthy people) per depth, and summarise.
set -eu
G=/projects/b1198/epifluidlab/jocelyn/software/cfdna-tissue-enrichment; O=/projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/prad_mix; PY=/home/pfu2606/.conda/envs/meth/bin/python
[ -s $O/samples_mixtures.csv ] || $PY -u $G/benchmarks/make_mixtures.py $O/pools $O $G
cd $G
for dep in 10 30; do
  $PY -c "import pandas as pd; d=pd.read_csv('$O/samples_mixtures.csv'); d[d.depth_M==$dep].to_csv('$O/samples_d${dep}M.csv', index=False)"
  $PY -m cfdna_tissue_enrichment run --panel tumor --build hg19 --samples $O/samples_d${dep}M.csv --cache $O/cache --mapq-col 0 --out $O/runs/d${dep}M \
      --modes sets,sample --group spiked --reference-level no --case yes --primary T_PRAD,C_neutrophil 2>&1 | grep -v Warning | tail -2
done
$PY $G/benchmarks/prad_summarise.py $O
