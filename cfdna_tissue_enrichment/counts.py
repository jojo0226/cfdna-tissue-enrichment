"""Fragments -> per-bin counts -> normalised coverage (identical to cfdna-enrich prep.py).

Fragment file: tab-separated chrom, start, end [, MAPQ, ...], plain or .gz. Fragment midpoints are counted in the
bins of the clean grid; fragments whose midpoint falls in a removed bin (blacklist, gap, N > 10%, mappability
< 0.80, non-autosomal) are dropped."""
import gzip
import numpy as np, pandas as pd

DETREND = {"50kb": 41, "10kb": 205}                                      # ~2 Mb rolling-median window, in bins


def count_fragments(path, lut, bs, nbins, min_mapq=30, lmin=50, lmax=500, mapq_col=None, chunk=4_000_000):
    counts = np.zeros(nbins, np.int64); hist = np.zeros(lmax + 1, np.int64); n_all = n_pass = 0
    opener = gzip.open if path.endswith(".gz") else open
    header = None
    with opener(path, "rt") as fh:
        first = fh.readline()
        while first.startswith("#") or first.lower().startswith(("chrom", "track", "browser")):
            if first.lower().lstrip("#").startswith("chrom"): header = [h.strip().lower() for h in first.lstrip("#").rstrip("\n").split("\t")]
            first = fh.readline()
    f = first.rstrip("\n").split("\t"); has_q = False; qc_ = None

    def _num(x):
        try:
            v = float(x); return 0 <= v <= 255
        except ValueError:
            return False
    if mapq_col is not None:                                             # explicit (1-based; 0 = none)
        has_q = mapq_col > 0; qc_ = mapq_col - 1
    elif header and "mapq" in header:                                    # e.g. FinaleDB: chrom start end name mapq strand
        qc_ = header.index("mapq"); has_q = True
    elif len(f) >= 4 and _num(f[3]):
        qc_, has_q = 3, True
    elif len(f) >= 5 and f[3] in (".", "") and _num(f[4]):               # BED6 with MAPQ in the score column
        qc_, has_q = 4, True
    use = [0, 1, 2] + ([qc_] if has_q else [])
    for df in pd.read_csv(path, sep="\t", header=None, usecols=use, comment="#", chunksize=chunk, dtype={0: str}, low_memory=False, on_bad_lines="skip"):
        df = df[pd.to_numeric(df[1], errors="coerce").notna() & pd.to_numeric(df[2], errors="coerce").notna()]
        n_all += len(df)
        ch = df[0].astype(str); ch = ch.where(ch.str.startswith("chr"), "chr" + ch).to_numpy()
        s = df[1].to_numpy(np.int64); e = df[2].to_numpy(np.int64); ln = e - s
        ok = (ln >= lmin) & (ln <= lmax)
        if has_q: ok &= pd.to_numeric(df[qc_], errors="coerce").fillna(0).to_numpy() >= min_mapq
        mid = (s + e) // 2
        for c in np.unique(ch[ok]):
            if c not in lut: continue
            m = ok & (ch == c); b = mid[m] // bs; a = lut[c]; b = b[b < len(a)]; r = a[b]; r = r[r >= 0]
            counts += np.bincount(r, minlength=nbins)
        hist += np.bincount(ln[ok], minlength=lmax + 1)[: lmax + 1]; n_pass += int(ok.sum())
    return counts, hist, n_all, n_pass, has_q


def preprocess(counts, g, binsz):
    """log2(CPM + 1) -> LOWESS on GC (frac 0.10, 3 iterations, level kept) -> per-chromosome centred ~2 Mb
    rolling-median detrend."""
    from statsmodels.nonparametric.smoothers_lowess import lowess
    cpm = counts / max(counts.sum(), 1) * 1e6; y = np.log2(cpm + 1); gc = g.gc.to_numpy(float); ok = np.isfinite(gc) & np.isfinite(y)
    fit = np.full(len(y), np.nan); fit[ok] = lowess(endog=y[ok], exog=gc[ok], frac=0.10, it=3, return_sorted=False)
    y = y - fit + np.nanmean(y[ok]); out = np.empty_like(y); w = DETREND[binsz]
    for _, idx in g.groupby("chrom", sort=False).indices.items():
        out[idx] = y[idx] - pd.Series(y[idx]).rolling(w, center=True, min_periods=1).median().to_numpy()
    return out


def qc_row(sid, counts, hist, n_all, n_pass, has_q, min_fragments):
    L = np.arange(len(hist)); tot = max(hist.sum(), 1)
    return dict(sample_id=sid, fragments_read=n_all, fragments_used=n_pass, fragments_in_grid=int(counts.sum()), mapq_filter=has_q,
                median_length=float(L[np.searchsorted(np.cumsum(hist), tot / 2)]), short_fraction_100_150=float(hist[100:151].sum() / tot),
                mono_fraction_151_220=float(hist[151:221].sum() / tot), bins_zero=float((counts == 0).mean()),
                mean_fragments_per_bin=float(counts.mean()), flag_low_depth=n_pass < min_fragments)
