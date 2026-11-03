#!/bin/bash
#SBATCH -A b1198
#SBATCH -p b1198
#SBATCH -t 4:00:00
#SBATCH --mem 32G
#SBATCH -c 2
#SBATCH -J cfte_build
#SBATCH -o /projects/b1198/epifluidlab/jocelyn/software/cfdna-tissue-enrichment/build/build_and_test_%j.log
# Maintainer job: build reference layers, make the example data, run the smoke test.
set -eu
G=/projects/b1198/epifluidlab/jocelyn/software/cfdna-tissue-enrichment; PY=/home/pfu2606/.conda/envs/meth/bin/python
cd /projects/b1198/epifluidlab/jocelyn/software/pregnancy_tissue_kit_build          # holds the extracted Wang h5ad
$PY -u $G/build/build_layers.py --out $G/reference/hg38
# examples: 8 x ~300 k fragments, random subsamples of the kit's BH01 example files (public healthy sample)
mkdir -p $G/examples/frag
for k in 1 2 3 4 5 6 7 8; do
  zcat /projects/b1198/epifluidlab/jocelyn/software/pregnancy_tissue_kit_build/pregnancy_tissue_kit_v1/example/frag/example_$k.bed.gz | awk -v s=$k 'BEGIN{srand(s)} rand()<0.3' | gzip -9 > $G/examples/frag/example_$k.bed.gz
done
sed 's#frag/#frag/#' /projects/b1198/epifluidlab/jocelyn/software/pregnancy_tissue_kit_build/pregnancy_tissue_kit_v1/example/samples.csv > $G/examples/samples.csv
du -sh $G/reference $G/examples
cd $G && export PATH=$(dirname $PY):$PATH && bash tests/smoke_test.sh
