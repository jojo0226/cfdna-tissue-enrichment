"""Reference loading: genome grid, atlas layers, panels, track specificity.

A *layer* is one per-bin openness matrix from one source and assay (e.g. C = ENCODE DNase, T = TCGA ATAC,
W = maternal-fetal-interface snATAC). Layers live in reference/<build>/layer_<ID>_<bin>.csv.gz with a `key`
column plus one column per track; reference/<build>/layers.json describes them.

A *panel* (panels/<name>.json) picks layers and adds derived tracks (groups, axes) and pre-registered primary
tracks. Every panel uses the same method; only the reference changes.

Specificity (the one definition used everywhere): within a layer, each bin's openness is standardised across
that layer's tracks (z), then each track's specificity = its z minus the mean z of the layer's other tracks.
Groups = mean specificity of member tracks; axes = mean(plus members) - mean(minus members)."""
import json, os, sys
import numpy as np, pandas as pd

PKG = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(PKG)
BIN = {"50kb": 50000, "10kb": 10000}


def ref_dir(build, root=None):
    d = os.path.join(root or os.environ.get("CFTE_REFERENCE", os.path.join(ROOT, "reference")), build)
    if not os.path.isdir(d):
        sys.exit(f"reference folder not found: {d} (set CFTE_REFERENCE or pass --reference)")
    return d


def list_panels():
    pd_ = os.path.join(ROOT, "panels")
    return sorted(f[:-5] for f in os.listdir(pd_) if f.endswith(".json"))


def load_panel(name_or_path):
    p = name_or_path if name_or_path.endswith(".json") else os.path.join(ROOT, "panels", f"{name_or_path}.json")
    if not os.path.exists(p):
        sys.exit(f"panel '{name_or_path}' not found; available: {list_panels()}")
    return json.load(open(p))


def load_grid(rdir, binsz):
    g = pd.read_csv(os.path.join(rdir, f"grid_{binsz}.csv.gz"))
    g = g.assign(_c=g.chrom.str.replace("chr", "").astype(int)).sort_values(["_c", "start"]).drop(columns="_c").reset_index(drop=True)
    return g


def lookup_tables(g, bs):
    lut = {}
    for c, d in g.groupby("chrom"):
        a = np.full(int(d.start.max() // bs) + 1, -1, np.int64); a[(d.start // bs).to_numpy()] = d.index.to_numpy(); lut[c] = a
    return lut


def load_layers(rdir, binsz, layer_ids, keys):
    meta = json.load(open(os.path.join(rdir, "layers.json")))
    out = {}
    for L in layer_ids:
        if L not in meta["layers"]:
            sys.exit(f"layer '{L}' not in {rdir}/layers.json (available: {list(meta['layers'])})")
        A = pd.read_csv(os.path.join(rdir, f"layer_{L}_{binsz}.csv.gz")).set_index("key").reindex(keys)
        out[L] = A.apply(pd.to_numeric, errors="coerce")
    return out, meta


def specificity(layers, panel):
    cols = {}
    for L, A in layers.items():
        Z = A.sub(A.mean(1), axis=0).div(A.std(1).replace(0, np.nan), axis=0); n = Z.shape[1]
        S = Z - (Z.sum(1).values[:, None] - Z.values) / (n - 1)
        for c in S.columns: cols[f"{L}_{c}"] = S[c].to_numpy()
    for name, g in panel.get("groups", {}).items():
        mem = [f"{g['layer']}_{m}" for m in g["members"] if f"{g['layer']}_{m}" in cols]
        if mem: cols[name] = np.nanmean(np.column_stack([cols[m] for m in mem]), 1)
    for name, ax in panel.get("axes", {}).items():
        lp, lm = ax.get("plus_layer", ax.get("layer")), ax.get("minus_layer", ax.get("layer"))
        P = [f"{lp}_{m}" for m in ax["plus"] if f"{lp}_{m}" in cols]
        M = [f"{lm}_{m}" for m in ax["minus"] if f"{lm}_{m}" in cols]
        if P and M:
            v = np.nanmean(np.column_stack([cols[m] for m in P]), 1) - np.nanmean(np.column_stack([cols[m] for m in M]), 1)
            cols[name] = v
            if ax.get("mirror"): cols[ax["mirror"]] = -v
    keys = next(iter(layers.values())).index
    return pd.DataFrame(cols, index=keys)


def primary_tracks(panel, override=None):
    if override:
        sep = ";" if ";" in override else ","                            # use ';' when track names contain commas
        return [t.strip() for t in override.split(sep) if t.strip()]
    return list(panel.get("primary", []))
