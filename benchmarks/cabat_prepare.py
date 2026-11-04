#!/usr/bin/env python
"""Cabat pan-cancer benchmark sheets (hg38): QC-pass samples, matched case/control design, ichorCNA tumour fraction.
Writes cabat_samples.csv (all), per-cancer sheets (that cancer's cases + all controls) and null sheets (controls only,
random fake labels, 10 repeats)."""
import os, sys
import numpy as np, pandas as pd
W = "/projects/b1198/epifluidlab/jocelyn/project/cfdna_wgs_pancancer_validation_20220225"; out = sys.argv[1]; os.makedirs(out, exist_ok=True)
t = pd.read_csv(f"{W}/samples_cabat.tsv", sep="\t"); t = t[t["pass"].astype(bool)]
ich = pd.read_csv(f"{W}/ichorcna/ichor_summary.csv")[["sample_id", "tumour_fraction"]]
t = t.merge(ich, on="sample_id", how="left")
s = pd.DataFrame({"sample_id": t.sample_id, "frag_path": t.frag_hg38, "status": t.case_control, "cancer_type": t.cancer_type.fillna("control"),
                  "batch": t.batch, "sex": t.sex, "pair_id": t.pair_id, "ichor_tf": t.tumour_fraction})
s.to_csv(f"{out}/cabat_samples.csv", index=False); print(len(s), s.status.value_counts().to_dict(), s.cancer_type.value_counts().to_dict())
for ct in sorted(s.cancer_type.unique()):
    if ct == "control": continue
    s[(s.cancer_type == ct) | (s.status == "control")].to_csv(f"{out}/sheet_{ct}.csv", index=False)
ctl = s[s.status == "control"].reset_index(drop=True); rng = np.random.default_rng(17)
for k in range(10):
    c = ctl.copy(); c["fake"] = np.where(rng.permutation(len(c)) < len(c) // 2, "A", "B"); c.to_csv(f"{out}/null_{k}.csv", index=False)
