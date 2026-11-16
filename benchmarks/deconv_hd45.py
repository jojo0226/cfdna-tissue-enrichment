"""Prototype: coverage-based cell-type deconvolution of one cfDNA sample, checked on HD45 against WGBS truth.

Model (per bin i, normalised coverage y_i from counts.preprocess):
    -y_i = c + b * mappability_i + g * mean_openness_i + sum_k beta_k * a_k(i) + e_i,   beta_k >= 0
a_k = openness of component k (ENCODE DNase layer C, quantile-normalised across tracks; components = single cell
types or means of related tissues). The intercept, mappability and the across-component mean openness (shared open
chromatin) are projected out first, so beta comes from cell-type-specific openness only. beta_k / sum(beta) is a
relative contribution: it assumes every source lowers coverage at its open bins by the same amount per unit of DNA.
Bins = union of the top-N specific bins of every component. CIs: bootstrap over bins.

usage: python deconv_hd45.py OUT_DIR REPO SAMPLE=FRAG [SAMPLE=FRAG ...]"""
import os, sys
import numpy as np, pandas as pd
from scipy.optimize import nnls

out, repo = sys.argv[1], sys.argv[2]; samples = dict(s.split("=", 1) for s in sys.argv[3:]); sys.path.insert(0, repo)
from cfdna_tissue_enrichment import reference as RF, counts as CN

BUILD, N_BOOT, SEED = "hg19", 200, 17
TRUTH = "/projects/b1198/epifluidlab/jocelyn/project/finaleme_hd45/truth/HD45_WGBS.deconv.hg38.tsv"
FINALEME = "/projects/b1198/epifluidlab/jocelyn/project/finaleme_hd45/hg19/HD45_FinaleMe.deconv.hg19.tsv"

# components: name -> ENCODE layer-C tracks (averaged); 'category' links to the WGBS (Loyfer) cell types below
COMP = {
    "B": ["B cell"], "Monocyte": ["CD14-positive monocyte"], "T": ["CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell"],
    "NK": ["natural killer cell"], "Neutrophil": ["neutrophil"], "Progenitor_CD34": ["common myeloid progenitor, CD34-positive"],
    "Liver": ["liver"], "Lung": ["lung", "left lung", "right lung", "upper lobe of left lung"],
    "Heart": ["heart", "heart left ventricle", "heart right ventricle", "left cardiac atrium", "right atrium auricular region"],
    "Skeletal_muscle": ["muscle of arm", "muscle of back", "muscle of leg", "muscle of trunk", "psoas muscle"],
    "Kidney": ["kidney", "left kidney", "right kidney"], "GI": ["stomach", "small intestine", "sigmoid colon", "transverse colon", "large intestine"],
    "Brain": ["brain", "cerebellum", "frontal cortex", "occipital lobe"], "Vascular": ["endothelial cell of umbilical vein", "dermis blood vessel endothelial cell", "coronary artery", "tibial artery"],
    "Adipose": ["omental fat pad"], "Skin": ["skin of body"], "Pancreas": ["pancreas"],
}
# WGBS (Loyfer atlas) cell type -> component used for comparison; unlisted types go to 'Other'
TRUTH_MAP = {"Blood-B": "B", "Blood-Mono+Macro": "Monocyte", "Blood-T": "T", "Blood-NK": "NK", "Blood-Granul": "Neutrophil", "Eryth-prog": "Progenitor_CD34",
             "Liver-Hep": "Liver", "Lung-Ep-Alveo": "Lung", "Lung-Ep-Bron": "Lung", "Heart-Cardio": "Heart", "Heart-Fibro": "Heart", "Skeletal-Musc": "Skeletal_muscle",
             "Kidney-Ep": "Kidney", "Gastric-Ep": "GI", "Colon-Ep": "GI", "Small-Int-Ep": "GI", "Neuron": "Brain", "Oligodend": "Brain", "Endothel": "Vascular",
             "Adipocytes": "Adipose", "Epid-Kerat": "Skin", "Pancreas-Acinar": "Pancreas", "Pancreas-Alpha": "Pancreas", "Pancreas-Beta": "Pancreas",
             "Pancreas-Delta": "Pancreas", "Pancreas-Duct": "Pancreas"}
BLOOD = ["B", "Monocyte", "T", "NK", "Neutrophil", "Progenitor_CD34"]


def components(A):
    Q = A.rank(axis=0, pct=True)                                         # quantile-normalise tracks (comparable scales)
    miss = [t for v in COMP.values() for t in v if t not in Q]
    if miss: sys.exit(f"tracks missing from layer C: {miss}")
    return pd.DataFrame({k: Q[v].mean(1) for k, v in COMP.items()})


def residualise(M, Z):
    B, *_ = np.linalg.lstsq(Z, M, rcond=None); return M - Z @ B


def solve(y, Xr, Z):
    yr = residualise(y[:, None], Z)[:, 0]; beta, _ = nnls(Xr, -yr); return beta


def read_deconv(path):
    d = pd.read_csv(path, sep="\t"); s = d.groupby(d.cell_type.map(TRUTH_MAP).fillna("Other")).proportion.sum()
    return s.reindex(list(COMP) + ["Other"]).fillna(0)


os.makedirs(out, exist_ok=True); rows, res = [], {}
for binsz, top in (("10kb", 5000), ("50kb", 1000)):
    rdir = RF.ref_dir(BUILD, os.path.join(repo, "reference")); g = RF.load_grid(rdir, binsz); lut = RF.lookup_tables(g, RF.BIN[binsz])
    layers, _ = RF.load_layers(rdir, binsz, ["C"], g.key.to_numpy()); X = components(layers["C"])
    Zs = X.sub(X.mean(1), axis=0); spec = Zs.div(Zs.std(0), axis=1)     # component specificity, for picking informative bins
    keep = np.zeros(len(g), bool)
    for k in X:
        v = spec[k].to_numpy(); ok = np.isfinite(v); keep[np.where(ok)[0][np.argsort(-v[ok])][:top]] = True
    for sid, frag in samples.items():
        cache = os.path.join(out, f"{sid}.{binsz}.npz")
        if os.path.exists(cache): cts = np.load(cache)["counts"]
        else:
            cts, hist, n_all, n_pass, _ = CN.count_fragments(frag, lut, RF.BIN[binsz], len(g), 30, 50, 500, 4); np.savez_compressed(cache, counts=cts)
            print(sid, binsz, f"{n_pass:,} fragments used", flush=True)
        y = CN.preprocess(cts, g, binsz)
        ok = keep & np.isfinite(y) & np.isfinite(X.to_numpy()).all(1) & g.mappability.notna().to_numpy()
        Xm = X.to_numpy()[ok]; Z = np.column_stack([np.ones(ok.sum()), g.mappability.to_numpy()[ok], Xm.mean(1)])
        Xr = residualise(Xm, Z); yv = y[ok]; beta = solve(yv, Xr, Z); rng = np.random.default_rng(SEED); boot = []
        for _ in range(N_BOOT):
            i = rng.integers(0, len(yv), len(yv)); boot.append(solve(yv[i], residualise(Xm[i], Z[i]), Z[i]))
        boot = np.array(boot); tot = boot.sum(1, keepdims=True); frac = np.divide(boot, tot, out=np.zeros_like(boot), where=tot > 0)
        p = beta / beta.sum() if beta.sum() > 0 else beta
        r = pd.DataFrame({"sample": sid, "bin": binsz, "component": list(X.columns), "beta": beta, "proportion": p,
                          "ci_lower": np.percentile(frac, 2.5, 0), "ci_upper": np.percentile(frac, 97.5, 0)})
        rows.append(r); res[(sid, binsz)] = r.set_index("component").proportion
        print(f"{sid} {binsz}: {ok.sum():,} bins, sum beta {beta.sum():.3f}", flush=True)

R = pd.concat(rows, ignore_index=True); R.to_csv(os.path.join(out, "deconv_coverage.csv"), index=False, float_format="%.4f")
truth = read_deconv(TRUTH); fme = read_deconv(FINALEME)
cmp_ = pd.DataFrame({"WGBS_truth": truth, "FinaleMe": fme})
for (sid, binsz), s in res.items(): cmp_[f"coverage_{sid}_{binsz}"] = s.reindex(cmp_.index).fillna(0)
cmp_.loc["Blood_total"] = cmp_.loc[BLOOD].sum()
cmp_.to_csv(os.path.join(out, "comparison_HD45.csv"), float_format="%.4f")


def agree(a, b):  # within-blood shares (renormalised to blood total), Pearson + mean absolute error
    a = a[BLOOD] / a[BLOOD].sum(); b = b[BLOOD] / b[BLOOD].sum(); return np.corrcoef(a, b)[0, 1], (a - b).abs().mean()


print("\n" + (cmp_ * 100).round(1).to_string())
print("\nwithin-blood agreement with WGBS truth (Pearson r, mean abs error in % of blood):")
for c in cmp_.columns[1:]:
    if "HD46" in c: continue
    r_, m_ = agree(cmp_[c], cmp_.WGBS_truth); print(f"  {c:<26} r {r_:.2f}  MAE {100 * m_:.1f}")
