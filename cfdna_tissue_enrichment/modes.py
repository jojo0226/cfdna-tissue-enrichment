"""The two modes. Both take the normalised coverage matrix Xp (samples x grid bins) and the track specificity
SP (grid bins x tracks) from reference.specificity().

  sample        per sample: bin z against a reference group, Spearman with each track's specificity over all bins
  differential  per bin: group effect (moderated t, default) or pooled-group difference (--diff-stat pooled); per
                track: Spearman of the bin statistic with specificity over all bins, and specificity of the top-N bins
                up / down; label-permutation nulls (top bins re-selected per permutation); closed-form random-bin z
                (descriptive)"""
import numpy as np, pandas as pd


def closed_form_z(m, Sm, n):
    """z of a top-n bin set's mean specificity against random n-bin sets from all bins (closed form, finite-population
    corrected; the cfdna-chromatin-xai spec_z). Descriptive: top bins are selected, so it is anti-conservative."""
    N = len(Sm); return (m - Sm.mean(0)) / (Sm.std(0) / np.sqrt(n) * np.sqrt((N - n) / (N - 1)))
from scipy import stats as sst
from .stats import ols_t, strata, permuted_design, bh


def match_strata(g):
    gc = pd.qcut(g.gc.rank(method="first"), 10, labels=False).to_numpy()
    if "mappability" in g:
        return gc * 5 + pd.qcut(g.mappability.rank(method="first"), 5, labels=False).to_numpy()
    return gc


def matched_random(idx, st, n_random, rng):
    cnt = np.bincount(st[idx], minlength=st.max() + 1); parts = []
    for k in np.nonzero(cnt)[0]:
        pool = np.where(st == k)[0]; rep = len(pool) < cnt[k]
        parts.append(np.stack([rng.choice(pool, cnt[k], replace=rep) for _ in range(n_random)]))
    return np.concatenate(parts, axis=1)


# ------------------------------------------------------------------ sample
def mode_sample(Xp, samples, SP, ref_mask, loo_min=4):
    ok = np.isfinite(SP.to_numpy()).all(1) & np.isfinite(Xp).all(0); X = Xp[:, ok]; Sm = SP.to_numpy()[ok]
    sr = np.apply_along_axis(sst.rankdata, 0, Sm); sr = (sr - sr.mean(0)) / sr.std(0)
    ri = np.where(ref_mask)[0]; Z = np.empty_like(X)
    for i in range(len(X)):
        rows = ri[ri != i] if (ref_mask[i] and len(ri) >= loo_min) else ri
        med = np.median(X[rows], 0); mad = 1.4826 * np.median(np.abs(X[rows] - med), 0)
        mad = np.where(mad > 0, mad, np.median(mad[mad > 0])); Z[i] = (X[i] - med) / mad
    zr = np.apply_along_axis(sst.rankdata, 1, Z); zr = (zr - zr.mean(1, keepdims=True)) / zr.std(1, keepdims=True)
    return pd.DataFrame(zr @ sr / ok.sum(), index=samples, columns=SP.columns), int(ok.sum()), bool(len(ri) >= loo_min)


# ------------------------------------------------------------------ differential
def mode_differential(Xp, samples, SP, X, terms, primary, g, n_perm, top_n, n_random, seed):
    rows_idx = [samples.index(s) for s in X.index]
    ok = np.isfinite(SP.to_numpy()).all(1) & np.isfinite(Xp).all(0)
    Y = Xp[rows_idx][:, ok]; Sm = SP.to_numpy()[ok]; Xa = X.to_numpy(float); cols = list(X.columns)
    sr = np.apply_along_axis(sst.rankdata, 0, Sm); sr = (sr - sr.mean(0)) / sr.std(0)
    _, T, d = ols_t(Xa, Y, moderate=True); rng = np.random.default_rng(seed); tcols = [X.columns.get_loc(t) for t in terms]
    gk = g[ok].reset_index(drop=True); st = match_strata(gk); gcv = gk.gc.to_numpy(float); mpv = gk.mappability.to_numpy(float) if "mappability" in gk else None
    n = min(top_n, ok.sum() // 4)

    def rho(tv):
        r = sst.rankdata(np.nan_to_num(tv)); r = (r - r.mean()) / r.std(); return r @ sr / len(r)

    def topspec(tv):
        tv = np.nan_to_num(tv); up = np.argpartition(-tv, n)[:n]; dn = np.argpartition(tv, n)[:n]
        return np.nanmean(Sm[up], 0), np.nanmean(Sm[dn], 0), up, dn
    sf, top, bins_out = [], [], []
    for t, j in zip(terms, tcols):
        tv = T[j]; obs = rho(tv); m_up, m_dn, o_up, o_dn = topspec(tv); grp = strata(X, tcols, j)
        null = np.empty((n_perm, Sm.shape[1])); nu = np.empty_like(null); nd = np.empty_like(null)
        for b in range(n_perm):
            tp = ols_t(permuted_design(Xa, cols, j, grp, rng), Y, moderate=True)[1][j]; null[b] = rho(tp); nu[b], nd[b], _, _ = topspec(tp)
        r = pd.DataFrame({"effect": t, "track": SP.columns, "rho": obs, "p_perm": (np.abs(null) >= np.abs(obs)[None] - 1e-12).mean(0)})
        r["primary"] = r.track.isin(primary); r["q_perm"] = bh(r.p_perm); r.loc[r.primary, "q_perm_primary"] = bh(r.loc[r.primary, "p_perm"]); sf.append(r)
        bins_out.append(pd.DataFrame({"effect": t, "key": g.key.to_numpy()[ok], "t": tv, "p": 2 * sst.t.sf(np.abs(tv), d)}))
        for direction, m, o, nl in (("up", m_up, o_up, nu), ("down", m_dn, o_dn, nd)):
            Rr = matched_random(o, st, n_random, rng); rnd = np.nanmean(Sm[Rr], 1)
            tt = pd.DataFrame({"effect": t, "direction": direction, "n_bins": n, "track": SP.columns, "mean_specificity": m, "perm_null_mean": nl.mean(0),
                               "p_perm": (np.abs(nl - nl.mean(0)) >= np.abs(m - nl.mean(0))[None] - 1e-12).mean(0),
                               "z_vs_matched_random": (m - rnd.mean(0)) / rnd.std(0, ddof=1), "z_closed_form": closed_form_z(m, Sm, n), "gc_top": gcv[o].mean(), "gc_random": gcv[Rr].mean(),
                               "map_top": mpv[o].mean() if mpv is not None else np.nan, "map_random": mpv[Rr].mean() if mpv is not None else np.nan})
            tt["primary"] = tt.track.isin(primary); tt["q_perm"] = bh(tt.p_perm); tt.loc[tt.primary, "q_perm_primary"] = bh(tt.loc[tt.primary, "p_perm"]); top.append(tt)
    return pd.concat(sf, ignore_index=True), pd.concat(top, ignore_index=True), pd.concat(bins_out, ignore_index=True), int(ok.sum()), d


# ------------------------------------------------------------------ differential, pooled groups (cfdna-chromatin-xai)
def mode_differential_pooled(Lev, samples, SP, X, term, primary, g, n_perm, top_n, n_random, seed, window, loo=False):
    """One two-level factor. Lev = per-sample log2(CPM + 1) before the detrend (GC-corrected or not, per --gc-correct).
    Per bin d = detrend(log2(mean CPM of cases + 1)) - detrend(log2(mean CPM of reference + 1)): samples are pooled,
    within-group variance is not used. Per track: Spearman of d with specificity (all bins); mean specificity of the
    top-N bins by d ('up' = more coverage in cases) and by -d ('down'); label-permutation nulls (re-pool and re-select
    per permutation); GC x mappability matched random z and closed-form z (descriptive). loo=True: drop each sample,
    re-pool, re-score the top bins."""
    from .counts import detrend
    rows_idx = [samples.index(s_) for s_ in X.index]; case = X[term].to_numpy(float) == 1
    L = Lev[rows_idx]; lin = 2.0 ** L - 1; ok = np.isfinite(SP.to_numpy()).all(1) & np.isfinite(L).all(0)
    Sm = SP.to_numpy()[ok]; sr = np.apply_along_axis(sst.rankdata, 0, Sm); sr = (sr - sr.mean(0)) / sr.std(0)
    rng = np.random.default_rng(seed); n = min(top_n, ok.sum() // 4)
    gk = g[ok].reset_index(drop=True); st = match_strata(gk); gcv = gk.gc.to_numpy(float); mpv = gk.mappability.to_numpy(float) if "mappability" in gk else None

    def dstat(cm, keep=None):
        k = np.ones(len(cm), bool) if keep is None else keep
        a = detrend(np.log2(np.nanmean(lin[k & cm], 0) + 1), g, window); b = detrend(np.log2(np.nanmean(lin[k & ~cm], 0) + 1), g, window)
        return (a - b)[ok]

    def rho(tv):
        r = sst.rankdata(np.nan_to_num(tv)); r = (r - r.mean()) / r.std(); return r @ sr / len(r)

    def topspec(tv):
        tv = np.nan_to_num(tv); up = np.argpartition(-tv, n)[:n]; dn = np.argpartition(tv, n)[:n]
        return np.nanmean(Sm[up], 0), np.nanmean(Sm[dn], 0), up, dn
    d = dstat(case); obs = rho(d); m_up, m_dn, o_up, o_dn = topspec(d)
    null = np.empty((n_perm, Sm.shape[1])); nu = np.empty_like(null); nd = np.empty_like(null)
    for b in range(n_perm):
        dp = dstat(rng.permutation(case)); null[b] = rho(dp); nu[b], nd[b], _, _ = topspec(dp)
    sf = pd.DataFrame({"effect": term, "track": SP.columns, "rho": obs, "p_perm": (np.abs(null) >= np.abs(obs)[None] - 1e-12).mean(0)})
    sf["primary"] = sf.track.isin(primary); sf["q_perm"] = bh(sf.p_perm); sf.loc[sf.primary, "q_perm_primary"] = bh(sf.loc[sf.primary, "p_perm"])
    top = []
    for direction, m, o, nl in (("up", m_up, o_up, nu), ("down", m_dn, o_dn, nd)):
        Rr = matched_random(o, st, n_random, rng); rnd = np.nanmean(Sm[Rr], 1)
        tt = pd.DataFrame({"effect": term, "direction": direction, "n_bins": n, "track": SP.columns, "mean_specificity": m, "perm_null_mean": nl.mean(0),
                           "p_perm": (np.abs(nl - nl.mean(0)) >= np.abs(m - nl.mean(0))[None] - 1e-12).mean(0),
                           "z_vs_matched_random": (m - rnd.mean(0)) / rnd.std(0, ddof=1), "z_closed_form": closed_form_z(m, Sm, n), "gc_top": gcv[o].mean(), "gc_random": gcv[Rr].mean(),
                           "map_top": mpv[o].mean() if mpv is not None else np.nan, "map_random": mpv[Rr].mean() if mpv is not None else np.nan})
        tt["primary"] = tt.track.isin(primary); tt["q_perm"] = bh(tt.p_perm); tt.loc[tt.primary, "q_perm_primary"] = bh(tt.loc[tt.primary, "p_perm"]); top.append(tt)
    lo = []
    if loo:
        for i, s_ in enumerate(X.index):
            k = np.ones(len(case), bool); k[i] = False
            if (k & case).sum() < 1 or (k & ~case).sum() < 1: continue
            mu_, md_, _, _ = topspec(dstat(case, k))
            for direction, m in (("up", mu_), ("down", md_)):
                lo.append(pd.DataFrame({"effect": term, "direction": direction, "dropped": s_, "track": SP.columns, "mean_specificity": m, "z_closed_form": closed_form_z(m, Sm, n)}))
    bins_out = pd.DataFrame({"effect": term, "key": g.key.to_numpy()[ok], "d": d})
    return sf, pd.concat(top, ignore_index=True), bins_out, int(ok.sum()), (pd.concat(lo, ignore_index=True) if lo else None)
