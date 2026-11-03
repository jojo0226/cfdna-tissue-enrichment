#!/bin/bash
# Runs every panel on the example data and checks the run is reproducible. ~3-5 minutes on one CPU.
set -euo pipefail
cd "$(dirname "$0")/.."
T=${TMPDIR:-/tmp}/cfte_smoke; rm -rf $T; mkdir -p $T
python -m cfdna_tissue_enrichment panels
for p in prenatal tumor glioma pan; do
  python -m cfdna_tissue_enrichment run --panel $p --samples examples/samples.csv --out $T/$p --mapq-col 0 \
    --factor birth:term:preterm --factor sga:nonSGA:SGA --covariates gestational_age_weeks --n-perm 500 --n-perm-bins 50 > $T/$p.log 2>&1 && echo "OK  $p" || { echo "FAIL $p"; tail -20 $T/$p.log; exit 1; }
done
python -m cfdna_tissue_enrichment run --panel prenatal --samples examples/samples.csv --out $T/prenatal_rerun --mapq-col 0 \
  --factor birth:term:preterm --factor sga:nonSGA:SGA --covariates gestational_age_weeks --n-perm 500 --n-perm-bins 50 > /dev/null 2>&1
for f in sets/scores_corrected.csv sample/scores_E.csv differential/enrichment_all_bins.csv; do cmp $T/prenatal/$f $T/prenatal_rerun/$f && echo "reproducible $f"; done
