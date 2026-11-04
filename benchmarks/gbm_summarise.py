#!/usr/bin/env python
"""Summarise GBM benchmarks and the PRAD differential run: primary-track results per mode. usage: gbm_summarise.py BENCH_DIR"""
import os, sys
import pandas as pd
B = sys.argv[1]; pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 40)
runs = [("gbm/runs/GBM_vs_schwannoma", "diagnosis:GBM_vs_Schwannoma"), ("gbm/runs/Breast_met_vs_schwannoma", "diagnosis:Breast_met_vs_Schwannoma"),
        ("gbm/runs/progression_vs_pseudo", "group:Progression_vs_Pseudoprogression"), ("prad_mix/runs/diff_d10M", "hi:yes_vs_no"), ("prad_mix/runs/diff_d30M", "hi:yes_vs_no")]
out = []
for d, eff in runs:
    for mode, f, col in (("sets", "sets/effect_tests.csv", "t"), ("sample", "sample/effect_tests.csv", "t"), ("differential", "differential/enrichment_all_bins.csv", "rho"),
                         ("diff_top_down", "differential/enrichment_top_bins.csv", "mean_specificity")):
        p = os.path.join(B, d, f)
        if not os.path.exists(p): continue
        e = pd.read_csv(p); e = e[e.effect == eff]
        if mode == "diff_top_down": e = e[e.direction == "down"]
        k = "score" if "score" in e else "track"
        T = e[e[k].str.startswith("T_")].copy(); T["rank_T"] = T[col].rank(ascending=True if col != "rho" else True)
        for _, r in e[e.primary].iterrows():
            rk = int(T.set_index(k)["rank_T"].get(r[k], -1)) if r[k].startswith("T_") else ""
            out.append(dict(run=d.split("/")[-1], mode=mode, track=r[k], stat=round(r[col], 4), p_perm=r.p_perm, q_primary=round(r.get("q_perm_primary", float("nan")), 3), rank_low_among_T=rk))
O = pd.DataFrame(out); O.to_csv(os.path.join(B, "gbm_prad_diff_summary.csv"), index=False); print(O.to_string(index=False))
