#!/usr/bin/env python
"""Summarise the Cabat benchmark: does each cancer type's own track stand out, does it track ichorCNA tumour
fraction, and are the nulls calibrated. usage: cabat_summarise.py CABAT_DIR"""
import glob, os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
O = sys.argv[1]; R = os.path.join(O, "runs"); S = pd.read_csv(os.path.join(O, "cabat_samples.csv")).set_index("sample_id")
TARGET = {"Prostate": ["T_PRAD"], "Colorectal": ["T_COAD"], "Liver": ["T_LIHC", "C_liver"], "Lung": ["T_LUAD", "T_LUSC", "Tgroup_lung"],
          "Breast": ["T_BRCA"], "Stomach": ["T_STAD"], "Pancreatic": ["C_pancreas"], "Ovary": ["C_ovary"]}
EFF = "status:case_vs_control"; rows = []; mats = {}
for ct, tg in TARGET.items():
    d = os.path.join(R, ct)
    if not os.path.isdir(d): continue
    cases = S.index[(S.cancer_type == ct)]
    for mode, f, col in (("sets", "sets/effect_tests.csv", "t"), ("sample", "sample/effect_tests.csv", "t"), ("differential", "differential/enrichment_all_bins.csv", "rho")):
        p = os.path.join(d, f)
        if not os.path.exists(p): continue
        e = pd.read_csv(p); e = e[e.effect == EFF]; key = "score" if "score" in e else "track"; e = e.set_index(key)
        tt = e[e.index.str.startswith("T_")][col]
        mats.setdefault(mode, {})[ct] = tt
        sc = os.path.join(d, "sets/scores_corrected.csv" if mode == "sets" else "sample/scores_E.csv")
        for t_ in tg:
            if t_ not in e.index: continue
            v = e.loc[t_, col]; pool = e[e.index.str.startswith(t_.split("_")[0] + "_")][col]
            r = dict(cancer=ct, n_cases=len(cases), mode=mode, track=t_, stat=v, p_perm=e.loc[t_, "p_perm"],
                     rank_in_layer=int((pool > v).sum() + 1), n_in_layer=len(pool), rank_abs=int((pool.abs() > abs(v)).sum() + 1))
            if mode != "differential" and os.path.exists(sc):
                X = pd.read_csv(sc, index_col=0); c = [s for s in cases if s in X.index and np.isfinite(S.loc[s, "ichor_tf"])]
                if len(c) >= 5:
                    rho, pr = spearmanr(X.loc[c, t_], S.loc[c, "ichor_tf"]); r.update(tf_spearman=rho, tf_p=pr, tf_n=len(c))
            rows.append(r)
T = pd.DataFrame(rows); T.to_csv(os.path.join(O, "benchmark_summary.csv"), index=False)
pd.set_option("display.width", 220); print(T.round(4).to_string(index=False))
# matrix figure: TCGA tracks (columns) x cancer types (rows), effect t (sets, sample) or rho (differential)
fig, axs = plt.subplots(1, len(mats), figsize=(7 * len(mats), 4.6), facecolor="white", squeeze=False)
for ax, (mode, m) in zip(axs[0], mats.items()):
    M = pd.DataFrame(m).T; M = M[sorted(M.columns)]; v = np.nanmax(np.abs(M.to_numpy())) or 1
    im = ax.imshow(M.to_numpy(), cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
    ax.set_xticks(range(M.shape[1])); ax.set_xticklabels([c[2:] for c in M.columns], rotation=90, fontsize=7); ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(M.index, fontsize=8)
    for i, ct in enumerate(M.index):
        for t_ in TARGET.get(ct, []):
            if t_ in M.columns: j = list(M.columns).index(t_); ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False, lw=1.6, ec="k"))
    fig.colorbar(im, ax=ax, fraction=.03); ax.set_title(f"{mode}: cases vs controls ({'rho' if mode == 'differential' else 't'}); box = own TCGA track", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(O, "benchmark_matrix.png"), dpi=150)
# nulls
nl = []
for mode, pat in (("sets", "null/run_*/sets/effect_tests.csv"), ("sample", "null/run_*/sample/effect_tests.csv"), ("differential", "null/diff_*/differential/enrichment_all_bins.csv")):
    for f in glob.glob(os.path.join(O, pat)):
        e = pd.read_csv(f); nl.append(dict(mode=mode, run=f.split("/")[-3], n_tests=len(e), frac_p_lt_005=(e.p_perm < .05).mean(), min_q=e.q_perm.min()))
N = pd.DataFrame(nl); N.to_csv(os.path.join(O, "null_calibration.csv"), index=False)
if len(N): print(N.groupby("mode")[["frac_p_lt_005", "min_q"]].agg(["mean", "max"]).round(3).to_string())
