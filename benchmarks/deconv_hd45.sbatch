#!/bin/bash
#SBATCH -A e30682
#SBATCH -p short
#SBATCH -t 04:00:00
#SBATCH --mem 48G
#SBATCH -c 2
#SBATCH -J cfte_deconv
#SBATCH -o /projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/logs/deconv_hd45_%j.log
# Coverage-based deconvolution prototype on HD45 (WGBS truth available) and HD46 (second healthy person, no truth).
set -eu
G=/projects/b1198/epifluidlab/jocelyn/software/cfdna-tissue-enrichment; P=/projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/prad_mix/pools
cd $G; /home/pfu2606/.conda/envs/meth/bin/python -u benchmarks/deconv_hd45.py /projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/deconv_hd45 $G \
  HD45=$P/HD45.frag.gz HD46=$P/HD46.frag.gz 2>&1 | grep -v Warning
