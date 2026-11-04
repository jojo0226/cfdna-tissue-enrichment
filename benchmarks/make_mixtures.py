#!/usr/bin/env python
"""In-silico mixtures for the PRAD benchmark, written straight into the tool's count cache (hg19 grid).
Tumour fraction TF in the mixture = share of PRAD fragments x TF0 (PRAD neat ichorCNA TF0 = 0.4651).
Background HD45 at TF in {0, .5, 1, 2, 5, 10, 20%}; HD46 (a second person) at TF 0 only. Depth 10 M and 30 M
fragments; 3 replicates (different seeds). Placeholder fragment files are written so 'run' accepts the sheet;
counts come from the cache (run with the same --mapq-col 0 / length filters).
usage: make_mixtures.py POOLDIR OUTDIR REPO_ROOT"""
import os, sys, gzip
import numpy as np, pandas as pd
pool, out, repo = sys.argv[1:4]; sys.path.insert(0, repo)
from cfdna_tissue_enrichment import reference as RF
TF0 = 0.4651; BIN = "50kb"; bs = RF.BIN[BIN]
g = RF.load_grid(RF.ref_dir("hg19", os.path.join(repo, "reference")), BIN); lut = RF.lookup_tables(g, bs)


def load(name):
    d = pd.read_csv(os.path.join(pool, f"{name}.frag.gz"), sep="\t", header=None, names=["c", "s", "e", "q"], dtype={"c": str})
    d["c"] = np.where(d.c.str.startswith("chr"), d.c, "chr" + d.c); ln = (d.e - d.s).to_numpy()
    keep = (ln >= 50) & (ln <= 500); d = d[keep]; ln = ln[keep]; mid = ((d.s + d.e) // 2).to_numpy(); row = np.full(len(d), -1)
    for c, idx in d.groupby("c").indices.items():
        if c in lut:
            b = mid[idx] // bs; a = lut[c]; ok = b < len(a); r = np.full(len(idx), -1); r[ok] = a[b[ok]]; row[idx] = r
    print(name, len(d), "fragments,", (row >= 0).mean().round(3), "in grid", flush=True)
    return row, ln
P = {n: load(n) for n in ("PRAD", "HD45", "HD46")}
cache = os.path.join(out, "cache"); os.makedirs(cache, exist_ok=True); os.makedirs(os.path.join(out, "frag"), exist_ok=True); rows = []
tag = f"hg19.{BIN}.q30.l50-500.m0"
for depth in (10_000_000, 30_000_000):
    for bg, tfs in (("HD45", [0, .005, .01, .02, .05, .10, .20]), ("HD46", [0])):
        for tf in tfs:
            for rep in (1, 2, 3):
                rng = np.random.default_rng(hash((depth, bg, tf, rep)) % 2**32); q = tf / TF0; nP = int(round(depth * q)); nB = depth - nP
                rr, ll = [], []
                for name, n in (("PRAD", nP), (bg, nB)):
                    if n == 0: continue
                    r, l = P[name]; i = rng.choice(len(r), n, replace=n > len(r)); rr.append(r[i]); ll.append(l[i])
                r = np.concatenate(rr); l = np.concatenate(ll)
                counts = np.bincount(r[r >= 0], minlength=len(g)); hist = np.bincount(l, minlength=501)[:501]
                sid = f"mix_{bg}_tf{tf * 100:g}_d{depth // 1_000_000}M_r{rep}"
                fp = os.path.join(out, "frag", f"{sid}.placeholder.bed"); open(fp, "w").close()
                np.savez_compressed(os.path.join(cache, f"{sid}.{tag}.npz"), counts=counts.astype(np.int32), hist=hist, n_all=depth, n_pass=depth, has_q=False, frag_path=fp)
                rows.append(dict(sample_id=sid, frag_path=fp, background=bg, tf_pct=tf * 100, depth_M=depth // 1_000_000, rep=rep, spiked="yes" if tf > 0 else "no",
                                 person=bg))
                print(sid, flush=True)
pd.DataFrame(rows).to_csv(os.path.join(out, "samples_mixtures.csv"), index=False)
