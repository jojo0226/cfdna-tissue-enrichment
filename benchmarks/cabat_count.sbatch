#!/bin/bash
#SBATCH -A b1198
#SBATCH -p b1198
#SBATCH -t 4:00:00
#SBATCH --mem 12G
#SBATCH -c 1
#SBATCH -J cfte_cabat_count
#SBATCH -o /projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/logs/cabat_count_%a.log
set -eu
G=/projects/b1198/epifluidlab/jocelyn/software/cfdna-tissue-enrichment; O=/projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/cabat
cd $G && /home/pfu2606/.conda/envs/meth/bin/python -m cfdna_tissue_enrichment count --samples $O/cabat_samples.csv --cache $O/cache --build hg38 --bin 50kb --index $SLURM_ARRAY_TASK_ID
