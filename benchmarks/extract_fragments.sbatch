#!/bin/bash
#SBATCH -A b1198
#SBATCH -p b1198
#SBATCH -t 24:00:00
#SBATCH --mem 16G
#SBATCH -c 8
#SBATCH -J cfte_extract
#SBATCH --array=0-2
#SBATCH -o /projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/logs/extract_%a_%j.log
# Subsampled fragment pools (hg19) from the deep PRAD / healthy BAMs (read-only source). Properly paired, MAPQ >= 30,
# primary, not duplicate; one line per fragment (from the leftmost mate): chrom start end mapq.
set -eu
O=/projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/prad_mix/pools; mkdir -p $O
D=/projects/b1198/epifluidlab/yaping/projects/ccinference/broad_data/deep_wgs_bam
SAMT=/gpfs/projects/b1198/epifluidlab/conda_envs/d087dc3f3e18fae57be60f688793b5c8_/bin/samtools
NAMES=(PRAD HD45 HD46); BAMS=(14230_1_WGS_all.merged.calmd.bam HD_45_WGS.calmd.bam HD_46_WGS.calmd.bam); FRAC=(17.10 17.50 17.35)   # seed 17, fraction .10 / .50 / .35
i=$SLURM_ARRAY_TASK_ID; n=${NAMES[$i]}
$SAMT view -@ 7 -s ${FRAC[$i]} -f 3 -F 3852 -q 30 $D/${BAMS[$i]} | awk -v OFS='\t' '$9 > 0 && $9 <= 1000 {print $3, $4 - 1, $4 - 1 + $9, $5}' | gzip -1 > $O/$n.frag.gz.tmp
mv $O/$n.frag.gz.tmp $O/$n.frag.gz; echo "$n $(zcat $O/$n.frag.gz | wc -l) fragments"
