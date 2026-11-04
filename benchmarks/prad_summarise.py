#!/usr/bin/env python
"""PRAD mixture benchmark: T_PRAD score vs spiked tumour fraction; detection = replicate score above the highest
unspiked sample (HD45 and HD46 resamples) at the same depth. usage: prad_summarise.py PRAD_DIR"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
O = sys.argv[1]; S = pd.read_csv(os.path.join(O, "samples_mixtures.csv")).set_index("sample_id"); rows = []; fig, axs = plt.subplots(2, 2, figsize=(11, 8), facecolor="white")
for i, dep in enumerate((10, 30)):
    for j, (mode, f) in enumerate((("sets", "sets/scores_corrected.csv"), ("sample", "sample/scores_E.csv"))):
        p = os.path.join(O, "runs", f"d{dep}M", f)
        if not os.path.exists(p): continue
        X = pd.read_csv(p, index_col=0); m = S.loc[X.index]; ax = axs[i, j]
        for t_ in [c for c in X.columns if c.startswith("T_")]:
            v = X[t_]; null = v[m.spiked == "no"]; thr = null.max(); sgn = np.sign(spearmanr(v[m.background == "HD45"], m.tf_pct[m.background == "HD45"])[0] or 1)
            rho = spearmanr(v[m.background == "HD45"], m.tf_pct[m.background == "HD45"])[0]
            det = {tf: float(((v[(m.tf_pct == tf)] * sgn) > (null * sgn).max()).mean()) for tf in sorted(m.tf_pct.unique()) if tf > 0}
            lod = min([tf for tf, fr in det.items() if fr == 1 and all(det[x] == 1 for x in det if x >= tf)], default=np.nan)
            rows.append(dict(depth_M=dep, mode=mode, track=t_, spearman_with_tf=rho, lod_pct=lod, **{f"detect_{k:g}pct": vv for k, vv in det.items()}))
        for t_, col in (("T_PRAD", "#d9622b"), ("C_neutrophil", "#3a78c2")):
            if t_ not in X: continue
            for bg, mk in (("HD45", "o"), ("HD46", "s")):
                s = m.background == bg; ax.scatter(m.tf_pct[s] + np.random.default_rng(1).uniform(-.1, .1, s.sum()), X.loc[s, t_], c=col, marker=mk, s=18, label=f"{t_} ({bg})")
        ax.set_xscale("symlog", linthresh=0.5); ax.set_xlabel("spiked tumour fraction (%)"); ax.set_title(f"{mode}, {dep} M fragments", fontsize=9); ax.legend(fontsize=6)
fig.tight_layout(); fig.savefig(os.path.join(O, "prad_dose_response.png"), dpi=150)
T = pd.DataFrame(rows); T.to_csv(os.path.join(O, "prad_summary.csv"), index=False)
pd.set_option("display.width", 220)
for (dep, mode), d in T.groupby(["depth_M", "mode"]):
    d = d.sort_values("spearman_with_tf", key=np.abs, ascending=False); print(f"\n== {mode}, {dep} M  (T_PRAD rank by |rho| among 23 TCGA tracks: {list(d.track).index('T_PRAD') + 1})")
    print(d.head(6).round(3).to_string(index=False))
