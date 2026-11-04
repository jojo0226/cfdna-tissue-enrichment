#!/bin/bash
#SBATCH -A b1198
#SBATCH -p b1198
#SBATCH -t 6:00:00
#SBATCH --mem 32G
#SBATCH -c 2
#SBATCH -J cfte_hg19
#SBATCH -o /projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/logs/build_hg19_%j.log
set -eu
G=/projects/b1198/epifluidlab/jocelyn/software/cfdna-tissue-enrichment; PY=/home/pfu2606/.conda/envs/meth/bin/python
$PY -u $G/build/build_layers_hg19.py --out $G/reference/hg19
cd $G && export PATH=$(dirname $PY):$PATH && bash tests/smoke_test.sh
