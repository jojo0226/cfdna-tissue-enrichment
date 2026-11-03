"""Group design and tests, shared by all modes.

Factors: --factor COLUMN:REFERENCE:CASE (two levels each). The reference group is every factor at its reference
level. Effects are estimated jointly with covariates (OLS); permutation p-values shuffle one factor's labels
within strata of the other factors."""
import itertools, sys
import numpy as np, pandas as pd
from scipy import stats as sst


def bh(p):
    p = np.asarray(p, float); out = np.full(len(p), np.nan); ok = np.isfinite(p); q = p[ok]; n = len(q)
    if n:
        o = np.argsort(q); r = np.empty(n); r[o] = np.minimum.accumulate((q[o] * n / np.arange(1, n + 1))[::-1])[::-1]; out[ok] = np.minimum(r, 1)
    return out


def parse_factors(specs, meta):
    F = []
    for s in specs:
        try:
            col, ref, case = s.split(":")
        except ValueError:
            sys.exit(f"--factor must be COLUMN:REFERENCE:CASE, got '{s}'")
        if col not in meta: sys.exit(f"--factor {s}: column '{col}' not in samples file")
        lv = set(meta[col].dropna().astype(str))
        if not {ref, case} <= lv: sys.exit(f"--factor {s}: levels {sorted(lv)} do not include both '{ref}' and '{case}'")
        F.append((col, ref, case))
    return F


def design_matrix(meta, F, covs, interaction):
    m = meta.copy(); keep = np.ones(len(m), bool)
    for col, ref, case in F: keep &= m[col].astype(str).isin([ref, case]).to_numpy()
    covs = [(c.split(":")[0], c.split(":")[1] if ":" in c else None) for c in covs]
    for c, _ in covs: keep &= m[c].notna().to_numpy()
    m = m[keep]; X = pd.DataFrame({"intercept": 1.0}, index=m.index); terms = []
    for col, ref, case in F:
        t = f"{col}:{case}_vs_{ref}"; X[t] = (m[col].astype(str) == case).astype(float); terms.append(t)
    if interaction and len(F) >= 2:
        for a, b in itertools.combinations(list(terms), 2): X[f"{a} x {b}"] = X[a] * X[b]; terms.append(f"{a} x {b}")
    for c, kind in covs:
        v = pd.to_numeric(m[c], errors="coerce")
        if kind == "num" or (kind is None and v.notna().all() and m[c].astype(str).nunique() >= 6):
            X[c] = v - v.mean()
        else:
            for lv in sorted(m[c].astype(str).unique())[1:]: X[f"{c}[{lv}]"] = (m[c].astype(str) == lv).astype(float)
    return X, terms, m


def cells(m, F):
    return m[[c for c, _, _ in F]].astype(str).agg("-".join, axis=1), "-".join(r for _, r, _ in F)


def ols_t(X, Y, moderate=False):
    """X n x p, Y n x k -> beta, t (p x k), residual df. moderate: shrink residual variances toward their median
    (prior 4 df), a simple stand-in for limma's empirical Bayes."""
    XtXi = np.linalg.pinv(X.T @ X); B = XtXi @ X.T @ Y; R = Y - X @ B; d = X.shape[0] - np.linalg.matrix_rank(X)
    s2 = (R ** 2).sum(0) / max(d, 1)
    if moderate:
        d0 = 4.0; s2 = (d0 * np.nanmedian(s2) + d * s2) / (d0 + d); d = d + d0
    se = np.sqrt(np.outer(np.diag(XtXi), s2)); return B, B / np.where(se > 0, se, np.nan), d


def strata(X, term_cols, j):
    others = [k for k in term_cols if k != j and " x " not in X.columns[k]]
    key = X.iloc[:, others].astype(int).astype(str).agg("".join, axis=1).to_numpy() if others else np.zeros(len(X), dtype="<U1")
    return [np.where(key == s)[0] for s in np.unique(key)]


def permuted_design(Xa, cols, j, groups, rng):
    Xp = Xa.copy(); v = Xa[:, j].copy()
    for i in groups: v[i] = v[rng.permutation(i)]
    Xp[:, j] = v
    for k, name in enumerate(cols):
        if " x " in name:
            a_, b_ = name.split(" x "); Xp[:, k] = Xp[:, cols.index(a_)] * Xp[:, cols.index(b_)]
    return Xp


def effect_tests(scores, X, terms, primary, n_perm, seed):
    Y = scores.loc[X.index].to_numpy(float); ok = np.isfinite(Y).all(0); Xa = X.to_numpy(float)
    B, T, d = ols_t(Xa, Y[:, ok]); rng = np.random.default_rng(seed); rows = []; tcols = [X.columns.get_loc(t) for t in terms]; cols = list(X.columns)
    for t, j in zip(terms, tcols):
        grp = strata(X, tcols, j); null = np.empty((n_perm, ok.sum()))
        for b in range(n_perm): null[b] = ols_t(permuted_design(Xa, cols, j, grp, rng), Y[:, ok])[1][j]
        r = pd.DataFrame({"score": scores.columns[ok], "effect": t, "beta": B[j], "t": T[j], "p_ols": 2 * sst.t.sf(np.abs(T[j]), d),
                          "p_perm": (np.abs(null) >= np.abs(T[j])[None] - 1e-12).mean(0)})
        r["primary"] = r.score.isin(primary); r["q_perm"] = bh(r.p_perm); r.loc[r.primary, "q_perm_primary"] = bh(r.loc[r.primary, "p_perm"])
        rows.append(r)
    return pd.concat(rows, ignore_index=True)


def pairwise_tests(scores, lab, ref, primary, n_perm, seed):
    rng = np.random.default_rng(seed); rows = []
    for cell in sorted(set(lab) - {ref}):
        idx = lab.index[lab.isin([cell, ref])]; y = (lab.loc[idx] == cell).to_numpy()
        if y.sum() < 1 or (~y).sum() < 1: continue
        V = scores.loc[idx].to_numpy(float); perms = np.stack([rng.permutation(y) for _ in range(n_perm)])
        for k, s in enumerate(scores.columns):
            v = V[:, k]; a, b = v[y], v[~y]
            if not (np.isfinite(a).all() and np.isfinite(b).all()): continue
            diff = a.mean() - b.mean(); null = (perms * v).sum(1) / y.sum() - ((~perms) * v).sum(1) / (~y).sum()
            rows.append(dict(comparison=f"{cell} vs {ref}", score=s, primary=s in primary, n_case=int(y.sum()), n_ref=int((~y).sum()), diff=diff,
                             auc=float(((a[:, None] > b[None]) + .5 * (a[:, None] == b[None])).mean()),
                             p_mannwhitney=float(sst.mannwhitneyu(a, b).pvalue), p_perm=float((np.abs(null) >= abs(diff) - 1e-12).mean())))
    T = pd.DataFrame(rows)
    if len(T):
        T["q_perm"] = T.groupby("comparison").p_perm.transform(bh)
        T.loc[T.primary, "q_perm_primary"] = T[T.primary].groupby("comparison").p_perm.transform(bh)
    return T
