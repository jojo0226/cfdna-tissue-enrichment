import pandas as pd, os, warnings; warnings.filterwarnings("ignore")
R = "/projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/cabat/runs"
blood = {"C_neutrophil": "Neutrophil", "C_CD14-positive monocyte": "Monocyte", "C_B cell": "B", "C_natural killer cell": "NK",
         "C_CD4-positive, alpha-beta T cell": "CD4 T", "C_CD8-positive, alpha-beta T cell": "CD8 T", "C_common myeloid progenitor, CD34-positive": "CMP",
         "C_spleen": "Spleen", "C_thymus": "Thymus", "C_liver": "Liver"}
rows = []
for ct in sorted(os.listdir(R)):
    for mode, f, col in (("sets", "sets/effect_tests.csv", "t"), ("sample", "sample/effect_tests.csv", "t"), ("differential", "differential/enrichment_all_bins.csv", "rho")):
        p = os.path.join(R, ct, f)
        if not os.path.exists(p): continue
        e = pd.read_csv(p); e = e[e.effect == "status:case_vs_control"]; e = e.set_index("score" if "score" in e else "track")
        for b, nm in blood.items():
            if b in e.index: rows.append(dict(cancer=ct, mode=mode, track=nm, stat=e.loc[b, col], p=e.loc[b, "p_perm"], q=e.loc[b, "q_perm"]))
T = pd.DataFrame(rows); T.to_csv("/projects/b1198/epifluidlab/jocelyn/project/cfte_benchmarks/cabat/immune_tracks.csv", index=False)
pd.set_option("display.width", 250)
for mode in ("sets", "sample", "differential"):
    d = T[T["mode"] == mode]
    M = d.pivot(index="cancer", columns="track", values="stat"); P = d.pivot(index="cancer", columns="track", values="p")
    S = M.round(3).astype(str) + P.apply(lambda c: c.map(lambda v: "*" if v < 0.05 else ""))
    print("\n== " + mode + " (" + ("rho" if mode == "differential" else "t") + "; * = permutation p < 0.05)"); print(S[list(blood.values())].to_string())
