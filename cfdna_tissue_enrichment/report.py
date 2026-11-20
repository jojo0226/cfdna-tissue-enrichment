"""Per-sample report for the sample and markers modes (the cfdna_chromatin.enrichment figure and test set, as used for the
C1D1 and Ahmed NU21C01 analyses): contrasts (paired sign-flip or unpaired AUC; contrasts.py), PCA (raw and / or
sample-centred) with loadings, clustered heatmap, effect plots with the random-set band (markers), strip plots of the
strongest features, PCA of per-patient change, lineage scores and set overlap (markers)."""
import os, re
import numpy as np, pandas as pd
from . import contrasts as CT, sample_figures as FG, markers as MK

IMM_P = [r"\bcd4-positive\b", r"\bcd8-positive\b", r"\bb cell\b", r"\bnatural killer cell\b", r"\bnk cell\b", r"\bneutrophil\b",
         r"\bmonocyte\b", r"\bmacrophage\b", r"\bgranulocyte\b", r"\bdendritic cell\b", r"\bmyeloid progenitor\b", r"\blymphocyte\b",
         r"\bleukocyte\b", r"\bpbmc\b", r"\bthymus\b", r"\bspleen\b", r"\bpeyer's patch\b", r"\bgerminal center\b"]
CNS_P = [r"\bbrain\b", r"\bcerebell", r"\bhippocamp", r"\bneuron", r"\bastrocy", r"\bmicroglia", r"\boligodend", r"\bopc\b",
         r"\bspinal cord\b", r"\bfrontal cortex\b", r"\boccipital\b", r"\btemporal lobe\b", r"\bparietal\b"]
NOT_CNS = [r"renal cortex"]
NICE = {"CD4-positive, alpha-beta T cell": "CD4+ T", "CD8-positive, alpha-beta T cell": "CD8+ T", "CD14-positive monocyte": "CD14+ monocyte",
        "common myeloid progenitor, CD34-positive": "CD34+ progenitor", "natural killer cell": "NK cell"}
TISSUE_ORDER = ["blood / immune", "CNS", "glioma (TCGA)", "other tumour (TCGA)", "other tissue"]


def track_category(t):
    """Display category of a tissue track (C_<ENCODE>, T_<TCGA>; other layers and derived tracks -> 'other tissue')."""
    if t.startswith("T_"):
        return "glioma (TCGA)" if t[2:] in ("GBM", "LGG") else "other tumour (TCGA)"
    x = (t[2:] if t.startswith("C_") else t).lower().replace("_", " ")
    if any(re.search(p, x) for p in IMM_P): return "blood / immune"
    if not any(re.search(p, x) for p in NOT_CNS) and any(re.search(p, x) for p in CNS_P): return "CNS"
    return "other tissue"


def track_nice(t):
    if t.startswith("T_"): return "TCGA_" + t[2:]
    x = t[2:] if t.startswith("C_") else t
    return NICE.get(x, x)


def order(feats, cat, levels):
    return sorted(feats, key=lambda f: (levels.index(cat[f]) if cat[f] in levels else 99, f))


def run(S, m, o, kind, cat, nice, ylab, out, rand=None, J=None, log=print):
    """S = samples x features (rows aligned with m). o: group_col, groups (ordered levels), patient_col, label_col,
    contrasts, connect, top_features, pca_centre."""
    os.makedirs(out, exist_ok=True); groups = o.groups
    st = FG.Style(m, o.group_col, groups, o.patient_col, o.label_col)
    cons = [CT.Contrast(m, c, r, o.group_col, o.patient_col) for c, r in CT.parse_contrasts(o.contrasts, groups)]
    S.assign(sample_id=m.sample_id.values, group=m[o.group_col].values).to_csv(os.path.join(out, "scores.csv"), index=False)
    tabs = []
    for con in cons:
        T = CT.table(S, con, cat)
        band = CT.random_band(con, rand, T.set_index("feature").effect) if rand else None
        if band is not None:
            T = T.merge(band, on="feature", how="left"); T["random_q"] = CT.bh(T.random_p.to_numpy())
        tabs.append(T)
        FG.fig_effect(T, con, cat, nice, out, f"{'Tissue enrichment' if kind == 'sample' else 'Marker-set coverage'}: {con.label()}", band)
        for _, r in T.sort_values("p").head(o.top_features).iterrows():
            FG.fig_strip(st, S[r.feature].to_numpy(), r.feature, nice, con, r, out, ylab)
        FG.fig_pca_change(st, S, con, out, f"PCA of per-patient change ({'tissue' if kind == 'sample' else 'markers'})")
    if tabs: pd.concat(tabs).to_csv(os.path.join(out, "contrast_tests.csv"), index=False)
    pcs = []
    for centred in {"no": [False], "yes": [True], "both": [False, True]}[o.pca_centre]:
        P, Vt, ev = FG.pca(S, centre_rows=centred); P, Vt = FG.orient(P, Vt, cons[0] if cons else None)
        lines = []
        for con in cons:
            r1, r2 = con.test(P[:, 0]), con.test(P[:, 1])
            lines.append(f"{con.label()}: PC1 p {r1['p']:.2g}, PC2 p {r2['p']:.2g}" + (f" ({r1['n_up']}/{r1['n']}, {r2['n_up']}/{r2['n']} up)" if con.kind == "paired" else ""))
            pcs += [dict(version="sample-centred" if centred else "raw", contrast=con.label(), pc=f"PC{j + 1}", var=ev[j], **con.test(P[:, j])) for j in range(2)]
        sfx = "_centred" if centred else ""; tag = " sample-centred\n(each sample's mean over features removed)" if centred else ""
        what = f"{S.shape[1]} {'atlas tracks' if kind == 'sample' else 'marker sets'}"
        FG.fig_pca(st, P, ev, "\n".join(lines), out, f"PCA of per-sample {'tissue enrichment' if kind == 'sample' else 'marker coverage'} ({what}){tag}",
                   connect=o.connect, fname=f"01_pca{sfx}.png")
        FG.fig_loadings(Vt, ev, list(S.columns), cat, nice, out, tag=", sample-centred" if centred else "", suffix=sfx)
        k = min(3, len(Vt))
        pd.DataFrame(Vt[:k].T, index=S.columns, columns=[f"PC{j + 1}" for j in range(k)]).assign(category=[cat.get(f, "") for f in S.columns]).to_csv(os.path.join(out, f"pca_loadings{sfx}.csv"))
        pd.DataFrame(P[:, :k], columns=[f"PC{j + 1}" for j in range(k)]).assign(sample_id=m.sample_id.values, group=m[o.group_col].values).to_csv(os.path.join(out, f"pca_scores{sfx}.csv"), index=False)
    pd.DataFrame(pcs).to_csv(os.path.join(out, "pca_tests.csv"), index=False)
    FG.fig_heatmap(st, S, cat, nice, out, f"Per-sample {'tissue enrichment' if kind == 'sample' else 'marker coverage'}, samples and features clustered")
    if kind == "markers":
        L, lcat = FG.lineage_scores(S, cat)
        n_by = {k: sum(1 for f in S.columns if (f if f.startswith("TCGA_") else cat.get(f)) == k) for k in L.columns}
        FG.fig_heatmap(st, L, lcat, lambda k: f"{k} ({n_by[k]} set{'s' if n_by[k] > 1 else ''})", out,
                       "Per-sample marker coverage by lineage (mean of each lineage's sets), samples clustered", fname="04b_heatmap_lineage.png", cluster_rows=False)
        L.assign(sample_id=m.sample_id.values, group=m[o.group_col].values).to_csv(os.path.join(out, "lineage_scores.csv"), index=False)
    if J is not None:
        FG.fig_overlap(J, nice, out)
    if tabs:
        log(pd.concat(tabs).sort_values("p")[["contrast", "feature", "category", "effect", "p", "q"] + [c for c in ("n_up", "random_p") if c in tabs[0]]].head(20).round(4).to_string(index=False))


def markers_categories(feats):
    return {f: MK.lineage(f) for f in feats}


MARKER_ORDER = list(MK.LINEAGE) + ["tumour (TCGA)", "other tumour (TCGA)", "other"]
