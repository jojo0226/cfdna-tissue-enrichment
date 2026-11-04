#!/usr/bin/env python
"""GBM benchmark sheets (hg38, FinaleDB fragments): 58-sample cohort (GBM, schwannoma, metastatic breast, pituitary
adenoma; one library batch) and the 10-sample progression vs pseudo-progression set. usage: gbm_prepare.py OUTDIR"""
import glob, os, re, sys
import pandas as pd
G = "/projects/b1198/epifluidlab/jocelyn/project/gbm"; out = sys.argv[1]; os.makedirs(out, exist_ok=True)
m = pd.read_csv(f"{G}/all_gbm_samples/sample_all_58_meta.txt", sep="\t")
fr = {re.sub(r"_S\d+\.hg38\.frag\.bed\.gz$", "", os.path.basename(f)): f for f in glob.glob(f"{G}/all_gbm_samples/frag/*.hg38.frag.bed.gz")}
s = pd.DataFrame({"sample_id": m.NuSeq_ID, "frag_path": m.NuSeq_ID.map(fr), "diagnosis": m.Diagnosis.str.replace(r"\s*\(Metastatic\)", "_met", regex=True).str.replace(" ", "_"),
                  "sex": m.Sex, "age": pd.to_numeric(m["Age at time of collection"], errors="coerce"), "pair_id": m.Pair_ID})
s = s[s.frag_path.notna()]; s.to_csv(f"{out}/gbm58_samples.csv", index=False); print(s.diagnosis.value_counts().to_dict(), "missing frag:", m.NuSeq_ID.map(fr).isna().sum())
for case in ("GBM", "Breast_met"):
    s[s.diagnosis.isin([case, "Schwannoma"])].to_csv(f"{out}/sheet_{case}_vs_schwannoma.csv", index=False)
p = pd.read_csv(f"{G}/finaledb_wgs_gbm_prog/gbm_metadata.csv", encoding="utf-8-sig")
fp = {re.sub(r"_S\d+\.hg38\.frag\.bed\.gz$", "", os.path.basename(f)): f for f in glob.glob(f"{G}/finaledb_wgs_gbm_prog/frag/*.hg38.frag.bed.gz")}
q = pd.DataFrame({"sample_id": p.NISeq_ID, "frag_path": p.NISeq_ID.map(fp), "group": p.Group, "sex": p.Sex, "age": p["Age at time of collection"]})
q = q[q.frag_path.notna()]; q.to_csv(f"{out}/prog_samples.csv", index=False); print(q.group.value_counts().to_dict())
