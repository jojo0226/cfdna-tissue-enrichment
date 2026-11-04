#!/usr/bin/env python
"""Build the hg19 reference (maintainers only): clean grid + layers C and T at 50 kb and 10 kb.

Grid: hg19 chr1-22 tiled in 50 kb / 10 kb bins. A bin is removed if it overlaps the hg19 ENCODE blacklist v2
(Amemiya et al. 2019) or a UCSC gap by any amount, has > 10% N, or mean CRG 50-mer mappability < 0.80
(hg38 uses a 45-mer track; the threshold is the same). GC = G+C over non-N bases.
Layers: C and T from the cfdna-chromatin-xai hg19 openness atlas (the build the atlas was originally made on, so
no liftover). W is not built for hg19.
Output: reference/hg19/{grid_<bin>.csv.gz, layer_C_<bin>.csv.gz, layer_T_<bin>.csv.gz, layers.json}."""
import argparse, gzip, json, os
import numpy as np, pandas as pd

R = "/projects/b1198/epifluidlab/jocelyn/software/cfdna-chromatin-xai/analysis/fetal_fraction"
H = "/projects/b1198/epifluidlab/jocelyn/project/ref/hg19"
D = dict(fasta=f"{H}/hg19.fa", blacklist=f"{H}/hg19-blacklist.v2.bed", gap=f"{H}/hg19_gap.txt.gz", mappability=f"{H}/wgEncodeCrgMapabilityAlign50mer.bigWig",
         c50=f"{R}/ff_openness_atlas_hg19_50kb_tissue73_tcga23.csv.gz", c10=f"{R}/ff_openness_atlas_hg19_10kb_tissue73_tcga23.csv.gz")
CHR = [f"chr{i}" for i in range(1, 23)]


def overlaps(starts, ends, iv):
    """Boolean per bin: overlaps any interval in iv (sorted DataFrame s, e)."""
    hit = np.zeros(len(starts), bool)
    if len(iv) == 0: return hit
    s, e = iv.s.to_numpy(), iv.e.to_numpy()
    lo = np.searchsorted(e, starts, side="right"); hi = np.searchsorted(s, ends, side="left")
    return hi > lo


def grid(bs, a):
    import pysam, pyBigWig
    fa = pysam.FastaFile(a.fasta); bw = pyBigWig.open(a.mappability)
    bl = pd.read_csv(a.blacklist, sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])
    gp = pd.read_csv(a.gap, sep="\t", header=None, usecols=[1, 2, 3], names=["c", "s", "e"])
    rows = []
    for c in CHR:
        L = fa.get_reference_length(c); st = np.arange(0, L, bs); en = np.minimum(st + bs, L)
        seq = np.frombuffer(fa.fetch(c).upper().encode(), dtype="S1")
        isn = seq == b"N"; isgc = (seq == b"G") | (seq == b"C")
        idx = np.arange(len(seq)) // bs; nN = np.bincount(idx, isn, len(st)); nGC = np.bincount(idx, isgc, len(st)); ln = (en - st).astype(float)
        mp = np.array([bw.stats(c, int(s), int(e), type="mean")[0] or 0.0 for s, e in zip(st, en)], float)
        b_ = overlaps(st, en, bl[bl.c == c].sort_values("s")); g_ = overlaps(st, en, gp[gp.c == c].sort_values("s"))
        d = pd.DataFrame({"chrom": c, "start": st, "end": en, "gc": nGC / np.maximum(ln - nN, 1), "N_fraction": nN / ln, "mappability": mp, "blacklist": b_, "gap": g_})
        rows.append(d); print(c, len(d), flush=True)
    g = pd.concat(rows, ignore_index=True); g["key"] = g.chrom + ":" + g.start.astype(str) + "-" + g.end.astype(str)
    keep = (g.N_fraction <= 0.10) & (g.mappability >= 0.80) & ~g.blacklist & ~g.gap & (g.end - g.start == bs)
    st = dict(total=len(g), N_gt_10pct=int((g.N_fraction > .10).sum()), map_lt_080=int((g.mappability < .80).sum()), blacklist=int(g.blacklist.sum()), gap=int(g.gap.sum()), kept=int(keep.sum()))
    return g[keep][["key", "chrom", "start", "end", "gc", "mappability"]].reset_index(drop=True), st


def main():
    ap = argparse.ArgumentParser(); [ap.add_argument(f"--{k}", default=v) for k, v in D.items()]; ap.add_argument("--out", required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); stats = {}
    for tag, bs, cp in (("50kb", 50000, a.c50), ("10kb", 10000, a.c10)):
        g, st = grid(bs, a); keys = g.key.to_numpy(); g.to_csv(os.path.join(a.out, f"grid_{tag}.csv.gz"), index=False, float_format="%.5g")
        C = pd.read_csv(cp).set_index("key").reindex(keys)
        cc = [c for c in C.columns if c.startswith("C_") and not c.startswith("C_TCGA_")]; tc = [c for c in C.columns if c.startswith("C_TCGA_")]
        C[cc].rename(columns=lambda c: c[2:]).reset_index().to_csv(os.path.join(a.out, f"layer_C_{tag}.csv.gz"), index=False, float_format="%.4g")
        C[tc].rename(columns=lambda c: c[len("C_TCGA_"):]).reset_index().to_csv(os.path.join(a.out, f"layer_T_{tag}.csv.gz"), index=False, float_format="%.4g")
        stats[tag] = dict(grid=st, bins=len(g), C_tracks=len(cc), T_tracks=len(tc), C_bins_covered=round(float(np.isfinite(C[cc].to_numpy(float)).all(1).mean()), 4))
        print(tag, stats[tag], flush=True)
    layers = {
        "C": dict(name="ENCODE DNase tissues and cell types", assay="DNase-seq (bulk)", n_tracks=stats["50kb"]["C_tracks"],
                  source="ENCODE portal DNase-seq narrowPeak (hg19), 72 biosamples + one neutrophil track; cfdna-chromatin-xai ff_openness atlas",
                  values="combined-openness z per track", provenance="native hg19 (the atlas's original build)", citation="ENCODE Project Consortium, Nature 2020"),
        "T": dict(name="TCGA tumour ATAC", assay="ATAC-seq (bulk tumour)", n_tracks=stats["50kb"]["T_tracks"],
                  source="TCGA-ATAC cancer-type-specific peak calls (Corces et al., Science 2018), 23 cancer types; same atlas file as C",
                  values="openness z per track", provenance="as C", citation="Corces et al., Science 362:eaav1898 (2018)")}
    json.dump(dict(build="hg19", bins=stats, grid="hg19 autosomes; ENCODE blacklist v2 (hg19), UCSC gaps, N>10%, CRG 50-mer mappability<0.80 removed", layers=layers),
              open(os.path.join(a.out, "layers.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
