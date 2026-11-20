#!/usr/bin/env python
"""Regression: the C1D1 and Ahmed NU21C01 per-sample analyses (originally cfdna_chromatin.enrichment `tissue` / `markers`)
rerun with `run --modes sample` (blood_cns_glioma panel, --sample-reference all) and `run --modes markers`, compared with
the original scores.csv / contrast_tests.csv. Usage: regression_c1d1_ahmed.py OUTDIR"""
import os, subprocess, sys
import numpy as np, pandas as pd
G = "/projects/b1198/epifluidlab/jocelyn/project/gbm"; T = sys.argv[1]; TOOL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MK_AH = f"{G}/Ahmed_9.22.26/reference/immune_markers_gbm/immune_marker_bins_10kb_immune_GBM_LGG.csv.gz"
CO = dict(c1d1=dict(sheet=f"{G}/c1d1_treatment/enrich/samples.csv", cov=f"{G}/c1d1_treatment/coverage_", orig=f"{G}/c1d1_treatment/enrich",
                    args=["--groups", "Responder,Progressor", "--patient-col", "none"], markers=MK_AH),  # same file as in its run.json, since moved to reference/
          ahmed=dict(sheet=f"{G}/Ahmed_9.22.26/by_timepoint/enrich/samples.csv", cov=f"{G}/Ahmed_9.22.26/features/coverage_", orig=f"{G}/Ahmed_9.22.26/by_timepoint/enrich",
                     args=["--groups", "Pre,Surgery,On-treatment", "--contrasts", "Surgery:Pre,On-treatment:Surgery,On-treatment:Pre", "--label-col", "cycle", "--connect"], markers=MK_AH))
NICE = {"CD4-positive, alpha-beta T cell": "CD4+ T", "CD8-positive, alpha-beta T cell": "CD8+ T", "CD14-positive monocyte": "CD14+ monocyte",
        "common myeloid progenitor, CD34-positive": "CD34+ progenitor", "natural killer cell": "NK cell"}
nice = lambda t: "TCGA_" + t[2:] if t.startswith("T_") else NICE.get(t[2:], t[2:])
res = []
for name, c in CO.items():
    for mode, bin_, sub, extra in (("sample", "50kb", "tissue_50kb", ["--panel", "blood_cns_glioma", "--sample-reference", "all"]),
                                   ("markers", "10kb", "markers_10kb", ["--markers", c["markers"]])):
        m = pd.read_csv(c["sheet"], dtype=str, keep_default_na=False); m["coverage"] = [f"{c['cov']}{bin_}/{s}_coverage.csv" for s in m.sample_id]
        os.makedirs(f"{T}/{name}", exist_ok=True); sh = f"{T}/{name}/samples_{bin_}.csv"; m.to_csv(sh, index=False); out = f"{T}/{name}/{sub}"
        cmd = [sys.executable, "-W", "ignore", "-m", "cfdna_tissue_enrichment", "run", "--samples", sh, "--out", out, "--modes", mode, "--bin", bin_] + extra + c["args"]
        p = subprocess.run(cmd, cwd=TOOL, capture_output=True, text=True); print("==", name, mode, "exit", p.returncode, flush=True)
        if p.returncode: print(p.stdout[-1500:], p.stderr[-3000:]); continue
        O = pd.read_csv(f"{c['orig']}/{sub}/scores.csv").set_index("sample_id"); rep = f"{out}/{mode}/report"
        N = pd.read_csv(f"{rep}/scores.csv").set_index("sample_id")
        if mode == "sample": N = N.rename(columns={k: nice(k) for k in N.columns})
        feats = [f for f in O.columns if f != "group"]; miss = [f for f in feats if f not in N]
        d = (N.loc[O.index, [f for f in feats if f in N]] - O[[f for f in feats if f in N]]).abs()
        res.append(dict(cohort=name, mode=mode, table="scores", n=len(O), features=len(feats), missing=len(miss), max_abs_diff=d.values.max(),
                        r=np.corrcoef(N.loc[O.index, feats].values.ravel(), O[feats].values.ravel())[0, 1] if not miss else np.nan))
        OT = pd.read_csv(f"{c['orig']}/{sub}/contrast_tests.csv"); NT = pd.read_csv(f"{rep}/contrast_tests.csv")
        if mode == "sample": NT["feature"] = NT.feature.map(nice)
        J = OT.merge(NT, on=["contrast", "feature"], suffixes=("_o", "_n"))
        for col in ("effect", "p", "q") + (("random_p",) if "random_p" in OT else ()):
            res.append(dict(cohort=name, mode=mode, table=f"contrast {col}", n=len(J), features=len(OT), missing=len(OT) - len(J),
                            max_abs_diff=(J[f"{col}_o"] - J[f"{col}_n"]).abs().max(), r=J[[f"{col}_o", f"{col}_n"]].corr().iloc[0, 1]))
pd.set_option("display.width", 200); print(pd.DataFrame(res).round(6).to_string(index=False))
