"""Contrast tests (sample and markers modes; from cfdna_chromatin.enrichment stats.py) on a per-sample score matrix (samples x features).

contrast case:ref. If a patient column links >= 3 patients across both groups the test is PAIRED: d = case - ref per
patient, effect = mean d, exact sign-flip p (all 2^n flips up to n = 16, else 20,000 random flips), n_up = patients with d > 0.
Otherwise UNPAIRED: effect = AUC (case > ref), p from the rank-sum (equivalent to |AUC - 0.5|), exact over every relabelling
when there are <= 20,000 of them (small cohorts: the smallest possible p is reported), else 10,000 permutations. BH q within contrast.
"""
from __future__ import annotations

import itertools
from math import comb

import numpy as np
import pandas as pd

SEED = 17


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def parse_contrasts(spec: str | None, groups: list[str]) -> list[tuple[str, str]]:
    if spec:
        return [tuple(c.split(":")) for c in spec.split(",")]
    return [(groups[i + 1], groups[i]) for i in range(len(groups) - 1)]       # adjacent groups, later minus earlier


def pairs(meta: pd.DataFrame, case: str, ref: str, group_col: str, patient_col: str | None) -> list[tuple[int, int, str]]:
    if not patient_col or patient_col not in meta:
        return []
    out = []
    for p, h in meta.groupby(patient_col):
        if pd.isna(p) or str(p) == "":
            continue
        a, b = h.index[h[group_col] == case], h.index[h[group_col] == ref]
        if len(a) and len(b):
            out.append((meta.index.get_loc(a[0]), meta.index.get_loc(b[0]), str(p)))
    return out if len(out) >= 3 else []


class Contrast:
    """Pre-computes the null design (flips or relabellings) once; .test(v) for any feature vector over all samples."""

    def __init__(self, meta, case, ref, group_col="group", patient_col=None, seed=SEED):
        self.case, self.ref = case, ref; g = meta[group_col].to_numpy()
        self.pairs = pairs(meta, case, ref, group_col, patient_col); rng = np.random.default_rng(seed)
        if self.pairs:
            self.kind = "paired"; n = len(self.pairs)
            self.F = np.array(list(itertools.product([1, -1], repeat=n))) if n <= 16 else rng.choice([1, -1], size=(20000, n))
            self.ia = np.array([a for a, _, _ in self.pairs]); self.ib = np.array([b for _, b, _ in self.pairs])
            self.min_p = 2 / len(self.F) if n <= 16 else 1 / len(self.F)
        else:
            self.kind = "unpaired"; self.idx = np.where(np.isin(g, [case, ref]))[0]; y = (g[self.idx] == case).astype(int); self.y = y
            n, k = len(y), int(y.sum())
            if comb(n, k) <= 20000:
                L = np.zeros((comb(n, k), n), int)
                for r, c in enumerate(itertools.combinations(range(n), k)):
                    L[r, list(c)] = 1
            else:
                L = np.array([rng.permutation(y) for _ in range(10000)])
            self.L = L; self.min_p = 2 / len(L) if comb(n, k) <= 20000 else 1 / len(L)

    def label(self):
        return f"{self.case} − {self.ref}" if self.kind == "paired" else f"{self.case} vs {self.ref}"

    def n_text(self):
        return f"{len(self.pairs)} paired patients" if self.kind == "paired" else f"{int(self.y.sum())} {self.case} v {int((1 - self.y).sum())} {self.ref}"

    def test(self, v) -> dict:
        v = np.asarray(v, float)
        if self.kind == "paired":
            d = v[self.ia] - v[self.ib]; o = d.mean(); null = (self.F * d).mean(1)
            return dict(effect=o, p=float((np.abs(null) >= abs(o) - 1e-12).mean()), n_up=int((d > 0).sum()), n=len(d))
        from scipy.stats import rankdata
        x = v[self.idx]; y = self.y; r = rankdata(x); c = (len(x) + 1) / 2
        o = r[y == 1].mean() - c; null = (self.L * r).sum(1) / y.sum() - c
        a = x[y == 1][:, None]; b = x[y == 0][None, :]; auc = float(((a > b) + .5 * (a == b)).mean())
        return dict(effect=auc, diff=float(x[y == 1].mean() - x[y == 0].mean()), p=float((np.abs(null) >= abs(o) - 1e-12).mean()), n=len(x))

    def effect_name(self):
        return f"mean paired change ({self.case} − {self.ref})" if self.kind == "paired" else f"AUC, {self.case} vs {self.ref} (1 = every {self.case} above every {self.ref})"

    def null_center(self):
        return 0.0 if self.kind == "paired" else 0.5


def table(S: pd.DataFrame, con: Contrast, category: dict[str, str] | None = None) -> pd.DataFrame:
    rows = [dict(feature=f, category=(category or {}).get(f, ""), **con.test(S[f].to_numpy())) for f in S.columns]
    T = pd.DataFrame(rows); T["q"] = bh(T.p.to_numpy()); T["contrast"] = con.label(); T["test"] = con.kind; T["min_possible_p"] = con.min_p
    return T


def random_band(con: Contrast, rand_scores: dict[str, np.ndarray], obs: pd.Series) -> pd.DataFrame:
    """rand_scores[feature] = samples x n_random matrix of random-set scores -> 95% band + competitive p."""
    rows = []
    for f, R in rand_scores.items():
        e = np.array([con.test(R[:, j])["effect"] for j in range(R.shape[1])]); c = np.median(e)
        rows.append(dict(feature=f, rand_lo=np.quantile(e, .025), rand_hi=np.quantile(e, .975), rand_median=c,
                         random_p=(1 + (np.abs(e - c) >= abs(obs[f] - c)).sum()) / (1 + len(e))))
    return pd.DataFrame(rows)
