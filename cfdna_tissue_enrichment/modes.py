"""The three modes. All take the normalised coverage matrix Xp (samples x grid bins) and the track specificity
SP (grid bins x tracks) from reference.specificity().

  sets          per sample: mean coverage over each track's top-N specific bins minus GC x mappability matched
                random bin sets (the 10-kb marker method of cfdna-enrich, with sets defined from the atlas)
  sample        per sample: bin z against a reference group, Spearman with each track's specificity over all bins
  differential  per bin: group effect (moderated t); per track: Spearman of bin t with specificity over all bins,
                and specificity of the top-N bins up / down; label-permutation nulls (top bins re-selected per
                permutation)"""
import numpy as np, pandas as pd
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


# ------------------------------------------------------------------ sets
def mode_sets(Xp, samples, sets, g, n_random, seed):
    st = match_strata(g); rng = np.random.default_rng(seed)
    R = {s: matched_random(idx, st, n_random, rng) for s, idx in sets.items()}
    raw, cor, zz = {}, {}, {}
    for i, sid in enumerate(samples):
        x = Xp[i]; raw[sid], cor[sid], zz[sid] = {}, {}, {}
        for s, idx in sets.items():
            obs = float(np.nanmean(x[idx])); rnd = np.nanmean(x[R[s]], 1)
            raw[sid][s] = obs; cor[sid][s] = obs - rnd.mean(); zz[sid][s] = (obs - rnd.mean()) / rnd.std(ddof=1)
    info = pd.DataFrame([dict(set=s, n_bins=len(idx), gc_set=g.gc.to_numpy()[idx].mean(), gc_random=g.gc.to_numpy()[R[s]].mean(),
                              map_set=g.mappability.to_numpy()[idx].mean() if "mappability" in g else np.nan,
                              map_random=g.mappability.to_numpy()[R[s]].mean() if "mappability" in g else np.nan) for s, idx in sets.items()])
    return (pd.DataFrame(d).T for d in (raw, cor, zz)), info


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
                               "z_vs_matched_random": (m - rnd.mean(0)) / rnd.std(0, ddof=1), "gc_top": gcv[o].mean(), "gc_random": gcv[Rr].mean(),
                               "map_top": mpv[o].mean() if mpv is not None else np.nan, "map_random": mpv[Rr].mean() if mpv is not None else np.nan})
            tt["primary"] = tt.track.isin(primary); tt["q_perm"] = bh(tt.p_perm); tt.loc[tt.primary, "q_perm_primary"] = bh(tt.loc[tt.primary, "p_perm"]); top.append(tt)
    return pd.concat(sf, ignore_index=True), pd.concat(top, ignore_index=True), pd.concat(bins_out, ignore_index=True), int(ok.sum()), d
