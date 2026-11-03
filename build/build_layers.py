#!/usr/bin/env python
"""Build the reference layers (maintainers only; users download / clone the result).

Output: reference/<build>/{grid_<bin>.csv.gz, layer_<ID>_<bin>.csv.gz, tracks_W.csv, layers.json} for bin = 50kb, 10kb.

Grid: clean hg38 autosomal bins from cfdna-chromatin-xai ms_xai (ENCODE blacklist ENCFF356LFX, gaps, N > 10% and
45-mer mappability < 0.80 removed), with GC and mappability per bin.

Layers (one source + assay each; values are per-bin openness, any scale -- specificity re-standardises per bin):
  C  ENCODE DNase-seq, 72 tissues / primary cells + neutrophil (cfdna-chromatin-xai ff_openness atlas, hg38;
     combined-openness z per track). TCGA columns are split off into T.
  T  TCGA ATAC-seq, 23 cancer types (same atlas file, columns C_TCGA_*).
  W  Wang et al. 2026 maternal-fetal interface pseudobulk snATAC (5 kb windows, hg19 file; each window id numbers
     the source hg38 5 kb tiling, so windows are placed on hg38 natively by midpoint). Per cell type: CPM -> log1p
     -> z across bins. Erythroid dropped (67 nuclei, 41% window coverage).

Add a new layer: write layer_<ID>_<bin>.csv.gz (key + one column per track, on the same grid) and an entry in
layers.json; then reference it from a panel. Nothing else changes."""
import argparse, collections, json, os, sys
import numpy as np, pandas as pd

R = "/projects/b1198/epifluidlab/jocelyn/software/cfdna-chromatin-xai"
D = dict(grid50=f"{R}/ms_xai/matrices/ms_50kb_reference_clean_final.bed", gc50=f"{R}/ms_xai/matrices/ms_50kb_bins_gc.csv", map50=f"{R}/ms_xai/matrices/ms_50kb_mappability.csv",
         grid10=f"{R}/ms_xai/matrices/ms_10kb_reference_clean_final.bed", gc10=f"{R}/ms_xai/matrices/ms_10kb_bins_gc.csv", map10=f"{R}/ms_xai/matrices/ms_10kb_mappability.csv",
         c50=f"{R}/analysis/fetal_fraction/ff_openness_atlas_hg38_50kb_tissue73_tcga23.csv.gz", c10=f"{R}/analysis/fetal_fraction/ff_openness_atlas_hg38_10kb_tissue73_tcga23.csv.gz",
         wang="wang_mfi_atlas_hg19/wang_mfi_win5kb_hg19.h5ad", chain="/projects/b1198/epifluidlab/jocelyn/project/ad/resources/hg19ToHg38.over.chain.gz")


def grid(bed, gc, mp):
    g = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "key"])
    g = g.merge(pd.read_csv(gc)[["key", "gc"]], on="key", how="left").merge(pd.read_csv(mp)[["key", "mappability"]], on="key", how="left")
    g = g[g.chrom.isin([f"chr{i}" for i in range(1, 23)]) & g.gc.notna()]
    g = g.assign(_c=g.chrom.str.replace("chr", "").astype(int)).sort_values(["_c", "start"]).drop(columns="_c").reset_index(drop=True)
    return g


def wang(path, chain):
    import h5py
    f = h5py.File(path, "r"); X = f["X"][:].astype(np.float64)
    dec = lambda a: [x.decode() if isinstance(x, bytes) else x for x in a]
    ob = pd.DataFrame({"track": dec(f["obs/cell_type"][:]), "origin": np.array(dec(f["obs/origin/categories"][:]))[f["obs/origin/codes"][:]],
                       "origin_purity": f["obs/origin_purity"][:], "low_n": f["obs/low_n_flag"][:], "n_nuclei": f["obs/n_cells_total"][:],
                       "total_fragments": f["obs/total_fragments"][:].astype(float)})
    n = np.array([int(x[1:]) for x in dec(f["var/_index"][:])]); cats = dec(f["var/chrom/categories"][:]); chrom = np.array(cats)[f["var/chrom/codes"][:]]; s19 = f["var/start"][:]
    sys.path.insert(0, os.path.join(R, "pylib")); from pyliftover import LiftOver
    lo = LiftOver(chain); rng = np.random.default_rng(0); base, agree = {}, {}
    for c in cats:
        m = np.where(chrom == c)[0]; cnt = collections.Counter()
        for i in rng.choice(m, min(600, len(m)), replace=False):
            l = lo.convert_coordinate(c, int(s19[i]))
            if l and l[0][0] == c: cnt[int(n[i]) - int(round(l[0][1] / 5000))] += 1
        b, k = cnt.most_common(1)[0]; base[c] = b; agree[c] = round(k / sum(cnt.values()), 3)
    s38 = np.array([(n[i] - base[chrom[i]]) * 5000 for i in range(len(n))])
    keep = (ob.track != "Ery").to_numpy()
    return X[keep], ob[keep].reset_index(drop=True), chrom, s38, agree


def main():
    ap = argparse.ArgumentParser(); [ap.add_argument(f"--{k}", default=v) for k, v in D.items()]; ap.add_argument("--out", required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    X, ob, chrom, s38, agree = wang(a.wang, a.chain)
    ob["clean_origin"] = (ob.origin_purity >= 0.90) & ~ob.low_n & ob.origin.isin(["Fetal", "Maternal"]); ob.to_csv(os.path.join(a.out, "tracks_W.csv"), index=False)
    stats = {}
    for tag, bs, bed, gc, mp, cp in (("50kb", 50000, a.grid50, a.gc50, a.map50, a.c50), ("10kb", 10000, a.grid10, a.gc10, a.map10, a.c10)):
        g = grid(bed, gc, mp); keys = g.key.to_numpy(); g.to_csv(os.path.join(a.out, f"grid_{tag}.csv.gz"), index=False, float_format="%.5g")
        C = pd.read_csv(cp).set_index("key").reindex(keys)
        cc = [c for c in C.columns if c.startswith("C_") and not c.startswith("C_TCGA_")]; tc = [c for c in C.columns if c.startswith("C_TCGA_")]
        C[cc].rename(columns=lambda c: c[2:]).reset_index().to_csv(os.path.join(a.out, f"layer_C_{tag}.csv.gz"), index=False, float_format="%.4g")
        C[tc].rename(columns=lambda c: c[len("C_TCGA_"):]).reset_index().to_csv(os.path.join(a.out, f"layer_T_{tag}.csv.gz"), index=False, float_format="%.4g")
        pos = {k: i for i, k in enumerate(g.chrom + ":" + g.start.astype(str))}
        rows = np.array([pos.get(f"{c}:{m}", -1) for c, m in zip(chrom, (s38 + 2500) // bs * bs)]); ok = rows >= 0
        M = np.stack([np.bincount(rows[ok], weights=X[j, ok], minlength=len(g)) for j in range(len(X))]); hit = np.bincount(rows[ok], minlength=len(g)) > 0
        W = np.log1p(M / ob.total_fragments.to_numpy()[:, None] * 1e6).T; W[~hit] = np.nan; W = (W - np.nanmean(W, 0)) / np.nanstd(W, 0)
        pd.concat([g[["key"]], pd.DataFrame(W, columns=ob.track)], axis=1).to_csv(os.path.join(a.out, f"layer_W_{tag}.csv.gz"), index=False, float_format="%.4g")
        stats[tag] = dict(bins=len(g), C_tracks=len(cc), T_tracks=len(tc), W_tracks=int(len(ob)), W_windows_placed=round(float(ok.mean()), 4), W_bins_covered=round(float(hit.mean()), 4))
        print(tag, stats[tag], flush=True)
    layers = {
        "C": dict(name="ENCODE DNase tissues and cell types", assay="DNase-seq (bulk)", n_tracks=stats["50kb"]["C_tracks"],
                  source="ENCODE portal DNase-seq narrowPeak, 72 biosamples (tissues, blood cells, endothelium, trophoblast) + one neutrophil track; cfdna-chromatin-xai ff_openness atlas",
                  values="combined-openness z per track", provenance="hg19 build lifted to hg38 by bin midpoint (10 kb); 50 kb hg38 provenance not documented", citation="ENCODE Project Consortium, Nature 2020"),
        "T": dict(name="TCGA tumour ATAC", assay="ATAC-seq (bulk tumour)", n_tracks=stats["50kb"]["T_tracks"],
                  source="TCGA-ATAC cancer-type-specific peak calls (Corces et al., Science 2018), 23 cancer types; same atlas file as C",
                  values="openness z per track", provenance="as C", citation="Corces et al., Science 362:eaav1898 (2018)"),
        "W": dict(name="Maternal-fetal interface snATAC", assay="snATAC-seq (pseudobulk)", n_tracks=stats["50kb"]["W_tracks"],
                  source="Wang et al. 2026, 32 cell types (erythroid dropped); 5 kb pseudobulk window counts prepared by Paz Polak",
                  values="CPM -> log1p -> z per cell type", provenance="windows placed on hg38 natively from their hg38 window numbers", hg38_offset_agreement=agree,
                  track_meta="tracks_W.csv", citation="Wang et al., Nature 653:167-179 (2026), doi:10.1038/s41586-026-10316-x")}
    json.dump(dict(build="hg38", bins=stats, grid="hg38 autosomes; ENCODE blacklist ENCFF356LFX, gaps, N>10%, 45-mer mappability<0.80 removed", layers=layers),
              open(os.path.join(a.out, "layers.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
