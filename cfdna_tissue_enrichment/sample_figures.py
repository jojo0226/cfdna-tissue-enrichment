"""Per-sample figure set (sample and markers modes; from cfdna_chromatin.enrichment figs.py): one figure set for any per-sample score matrix (tissue enrichment E or marker-set coverage).

Single-panel PNGs, same style as the cfdna-frag figures (fragmentomics.plots palette / ink / grid):
  01_pca.png                    every sample is a dot (shape = group; colour = patient when patients repeat, else group)
  02_pc1_loadings.png, 03_...   loadings coloured by feature category
  04_heatmap.png                samples x features, both clustered (Ward) with dendrograms; group bar on top
  05_effect_<contrast>.png      per-feature effect, sorted most positive -> most negative, p / q (and random-set band)
  06_strip_<contrast>_<feat>    the strongest features: every sample, all groups, paired lines when patients repeat
  07_pca_change_<contrast>.png  paired contrasts only: PCA of each patient's change (one dot per PATIENT, by design)
"""
from __future__ import annotations

import os
import re

import numpy as np
import pandas as pd

# plot style (cfdna_chromatin.fragmentomics.plots)
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def _setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                         "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titleweight": "bold", "axes.titlesize": 10,
                         "axes.titlelocation": "left", "legend.frameon": False})
    return plt


CAT_COLORS = {"blood / immune": PALETTE[0], "CNS": PALETTE[1], "glioma (TCGA)": PALETTE[5], "other tissue": "#9b9a95",
              "other tumour (TCGA)": "#8c6d31", "B": PALETTE[0], "CD4 T": PALETTE[1], "CD8 / γδ T": PALETTE[2], "NK": PALETTE[3],
              "myeloid": PALETTE[4], "tumour (TCGA)": PALETTE[5], "CNS cells": PALETTE[6], "microglia": PALETTE[7], "other": "#9b9a95"}
MARKERS = ["o", "s", "^", "D", "X", "v", "P", "*"]
PATIENT_PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948", "#8c564b", "#17becf"]
SINGLE = "#dcdbd6"


def cycle_colors(labels) -> dict:
    """Treatment-cycle labels (C1D1 ... C6D7, EOT, baseline ...) -> one sequential ramp in treatment order."""
    import matplotlib
    vals = sorted({str(x) for x in labels if str(x) not in ("", "nan")}, key=_cycle_key)
    cm = matplotlib.colormaps["YlGnBu"]
    return {v: cm(.25 + .7 * k / max(1, len(vals) - 1)) for k, v in enumerate(vals)}


def _cycle_key(v):
    m = re.match(r"C(\d+)D(\d+)", v)
    if m:
        return (1, int(m.group(1)), int(m.group(2)), v)
    return (0 if v.lower() in ("baseline", "screening", "pre") else 2, 0, 0, v)


def annotation_bars(fig, bars, order, x0, width, y0, bh, legend_x, legend_y, fontsize=7):
    """Stacked per-sample colour bars (top to bottom) + one legend per bar to the right.
    bars = [(title, colours per sample (list), legend [(label, colour)])]; order = sample order on the x axis."""
    from matplotlib.patches import Patch, Rectangle
    n = len(order); y = y0
    for title, cols, leg in bars:
        ax = fig.add_axes([x0, y, width, bh]); ax.set_xlim(-.5, n - .5); ax.set_ylim(0, 1)
        for x, i in enumerate(order):
            ax.add_patch(Rectangle((x - .5, 0), 1, 1, facecolor=cols[i], edgecolor="white", lw=.4))
        ax.set_xticks([]); ax.set_yticks([.5]); ax.set_yticklabels([title], fontsize=fontsize)
        ax.tick_params(left=False); [sp.set_visible(False) for sp in ax.spines.values()]
        y -= bh * 1.25
    ly = legend_y
    for title, cols, leg in bars:
        lg = fig.legend(handles=[Patch(facecolor=c, edgecolor=INK2, lw=.3, label=l) for l, c in leg], loc="upper left", bbox_to_anchor=(legend_x, ly),
                        fontsize=fontsize, title=title, title_fontsize=fontsize, ncol=1 if len(leg) <= 8 else 2)
        fig.add_artist(lg); ly -= (.035 + .022 * (len(leg) if len(leg) <= 8 else (len(leg) + 1) // 2))


def _safe(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")


class Style:
    def __init__(self, meta, group_col, groups, patient_col=None, label_col=None):
        import matplotlib
        self.meta, self.gcol, self.groups = meta.reset_index(drop=True), group_col, groups
        self.g = self.meta[group_col].astype(str).to_numpy()
        self.gcolor = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(groups)}
        self.gmark = {g: MARKERS[i % len(MARKERS)] for i, g in enumerate(groups)}
        pts = self.meta[patient_col].astype(str).where(self.meta[patient_col].notna(), "") if patient_col and patient_col in self.meta else pd.Series([""] * len(self.meta))
        self.patient = pts.to_numpy(); vc = pd.Series(self.patient[self.patient != ""]).value_counts()
        self.repeated = [p for p in vc.index if vc[p] > 1]
        tab = matplotlib.colormaps["tab20"]
        self.pcolor = {p: (PATIENT_PAL[i] if i < len(PATIENT_PAL) else tab(i % 20)) for i, p in enumerate(sorted(self.repeated))}
        if label_col and label_col in self.meta:
            lab = self.meta[label_col].astype(str).where(self.meta[label_col].notna() & (self.meta[label_col].astype(str) != ""), "")
        else:
            lab = pd.Series([""] * len(self.meta))
        sid = self.meta.sample_id.astype(str).str.extract(r"([^_]+)$")[0]
        self.short = np.array([self._pshort(p) if p in self.repeated else s for p, s in zip(self.patient, sid)])
        self.label = lab.to_numpy()
        self.ccolor = cycle_colors(self.label)

    def bars(self):
        """Annotation bars: group, patient (when patients repeat), cycle / label (when given)."""
        gc = self.gcolor
        if self.repeated:                     # patient colours take the categorical palette; groups get an ordinal grey ramp
            ramp = ["#d9d8d3", "#a3a29d", "#62615d", "#2b2a28", "#000000"]
            gc = {g: ramp[min(k, len(ramp) - 1)] if len(self.groups) > 1 else ramp[1] for k, g in enumerate(self.groups)}
        out = [("group", [gc[x] for x in self.g], [(x, gc[x]) for x in self.groups])]
        if self.repeated:
            out.append(("patient", [self.pcolor.get(p, SINGLE) for p in self.patient],
                        [(self._pshort(p), self.pcolor[p]) for p in sorted(self.repeated)] + [("single sample", SINGLE)]))
        if self.ccolor:
            out.append(("cycle", [self.ccolor.get(str(x), "#ffffff") for x in self.label], list(self.ccolor.items()) + [("none (Pre / Surgery)", "#ffffff")]))
        return out

    @staticmethod
    def _pshort(p):
        return re.sub(r"^\d+-\d+-", "", p)

    def color(self, i):
        return self.pcolor[self.patient[i]] if self.repeated and self.patient[i] in self.pcolor else (
            SINGLE if self.repeated else self.gcolor[self.g[i]])

    def tick(self, i):
        extra = self.label[i] or self.g[i]
        return f"{self.short[i]}\n{extra}"


def _save(fig, out, name):
    fig.savefig(os.path.join(out, name), dpi=200, bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)


def pca(S: pd.DataFrame, centre_rows=False):
    X = S.to_numpy(float)
    if centre_rows:
        X = X - X.mean(1, keepdims=True)
    sd = X.std(0); X = (X - X.mean(0)) / np.where(sd > 0, sd, 1)
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    return U * s, Vt, s ** 2 / (s ** 2).sum()


def orient(P, Vt, con):
    """Flip PC signs so the contrast's case group sits on the + side (sign is arbitrary in PCA)."""
    if con is None:
        return P, Vt
    for j in range(min(3, P.shape[1])):
        r = con.test(P[:, j])
        e = r["effect"] - con.null_center()
        if e < 0:
            P[:, j] *= -1; Vt[j] *= -1
    return P, Vt


def fig_pca(st: Style, P, ev, tests, out, title, connect=False, fname="01_pca.png"):
    plt = _setup(); fig, ax = plt.subplots(figsize=(7.6, 6.2))
    if connect and st.repeated:
        order = {g: k for k, g in enumerate(st.groups)}
        for p in st.repeated:
            ii = sorted(np.where(st.patient == p)[0], key=lambda i: order.get(st.g[i], 99))
            for a, b in zip(ii[:-1], ii[1:]):
                ax.annotate("", P[b, :2], P[a, :2], arrowprops=dict(arrowstyle="->", color=st.pcolor[p], lw=1.3, alpha=.75))
    for g in st.groups:
        for i in np.where(st.g == g)[0]:
            ax.scatter(P[i, 0], P[i, 1], marker=st.gmark[g], s=70, color=st.color(i), edgecolor=INK, lw=.5, zorder=3)
    for i in range(len(st.g)):
        t = st.short[i] + (f" {st.label[i]}" if st.label[i] else "")
        ax.annotate(t, (P[i, 0], P[i, 1]), xytext=(5, 3), textcoords="offset points", fontsize=6.5, color=INK2)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker=st.gmark[g], ls="", mfc="white" if st.repeated else st.gcolor[g], mec=INK, ms=7, label=f"{g} (n={(st.g == g).sum()})") for g in st.groups]
    if st.repeated:
        h += [Line2D([], [], marker="o", ls="", mfc=st.pcolor[p], mec=INK, ms=6, label=st._pshort(p)) for p in sorted(st.repeated)]
        h += [Line2D([], [], marker="o", ls="", mfc=SINGLE, mec=INK, ms=6, label="single sample")]
    ax.legend(handles=h, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1), title="shape = group" + ("\ncolour = patient" if st.repeated else ""), title_fontsize=7)
    ax.set_xlabel(f"PC1 ({ev[0] * 100:.0f}%)"); ax.set_ylabel(f"PC2 ({ev[1] * 100:.0f}%)"); ax.grid(color=GRID, lw=.5); ax.set_axisbelow(True)
    ax.set_title(title + ("\n" + tests if tests else ""), fontsize=9)
    _save(fig, out, fname)


def fig_loadings(Vt, ev, feats, cat, nice, out, tag="", suffix=""):
    plt = _setup(); from matplotlib.patches import Patch
    for j in range(2):
        l = pd.Series(Vt[j], index=feats).sort_values(); y = np.arange(len(l))
        fig, ax = plt.subplots(figsize=(5.8, max(4, .26 * len(l) + 1.2)))
        ax.barh(y, l.to_numpy(), color=[CAT_COLORS.get(cat.get(f, ""), "#9b9a95") for f in l.index])
        ax.set_yticks(y); ax.set_yticklabels([nice(f) for f in l.index], fontsize=8.5); ax.axvline(0, color=INK2, lw=.8)
        cats = list(dict.fromkeys(cat.get(f, "") for f in feats))
        ax.legend(handles=[Patch(color=CAT_COLORS.get(c, "#9b9a95"), label=c) for c in cats if c], fontsize=7, loc="lower right")
        ax.set_xlabel(f"PC{j + 1} loading"); ax.set_title(f"PC{j + 1} loadings ({ev[j] * 100:.0f}% of variance){tag}", fontsize=10)
        _save(fig, out, f"{2 + j:02d}_pc{j + 1}_loadings{suffix}.png")


def fig_heatmap(st: Style, S: pd.DataFrame, cat, nice, out, title, cbar="score (z across samples)", fname="04_heatmap.png", cluster_rows=True):
    from matplotlib.patches import Patch
    from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
    plt = _setup(); X = S.to_numpy(float); sd = X.std(0); Z = (X - X.mean(0)) / np.where(sd > 0, sd, 1)
    n, m = Z.shape; Lc = linkage(Z, "ward"); co = leaves_list(Lc)
    Lr = linkage(Z.T, "ward") if cluster_rows else None; ro = leaves_list(Lr) if cluster_rows else np.arange(m)
    bars = st.bars(); nb = len(bars)
    W = max(7, .32 * n + 5.2); H = max(5, .28 * m + 3) + .3 * (nb - 1); fig = plt.figure(figsize=(W, H))
    L, Wd = 2.6 / W, (W - 5.4) / W; bh = .16 / H; top = 1 - (1.3 + .2 * nb) / H; hh = top - 1.1 / H
    axt = fig.add_axes([L, top + bh * 1.25 * nb, Wd, .85 / H]); dendrogram(Lc, ax=axt, no_labels=True, color_threshold=0, above_threshold_color=INK2); axt.axis("off")
    annotation_bars(fig, bars, co, L, Wd, top + bh * 1.25 * (nb - 1) + .02 / H, bh, L + Wd + .2 / W, 1 - .5 / H)
    if cluster_rows:
        axl = fig.add_axes([.01, 1.1 / H, L - 1.75 / W, hh]); dendrogram(Lr, ax=axl, orientation="left", no_labels=True, color_threshold=0, above_threshold_color=INK2); axl.axis("off"); axl.invert_yaxis()
    axh = fig.add_axes([L, 1.1 / H, Wd, hh]); im = axh.imshow(np.clip(Z[co][:, ro].T, -3, 3), aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3, interpolation="nearest")
    axh.set_xticks(range(n)); axh.set_xticklabels([st.tick(i) for i in co], fontsize=6.5, rotation=90)
    axh.set_yticks(range(m)); axh.set_yticklabels([nice(S.columns[j]) for j in ro], fontsize=8)
    for t, j in zip(axh.get_yticklabels(), ro):
        t.set_color(CAT_COLORS.get(cat.get(S.columns[j], ""), INK))
    cax = fig.add_axes([L + Wd + 2.3 / W, 1.1 / H, .15 / W, min(hh, 2.2 / H)]); fig.colorbar(im, cax=cax).set_label(cbar, fontsize=8)
    fig.suptitle(title, x=L, ha="left", y=1 - .15 / H, fontsize=10, fontweight="bold")
    _save(fig, out, fname)


def fig_effect(T: pd.DataFrame, con, cat, nice, out, title, band: pd.DataFrame | None = None):
    plt = _setup(); from matplotlib.patches import Patch
    T = T.sort_values("effect").reset_index(drop=True); y = np.arange(len(T))
    fig, ax = plt.subplots(figsize=(7.6, max(4, .27 * len(T) + 1.6)))
    if band is not None:
        b = band.set_index("feature").reindex(T.feature)
        ax.hlines(y, b.rand_lo, b.rand_hi, color=GRID, lw=6, zorder=1)
    ax.scatter(T.effect, y, s=46, c=[CAT_COLORS.get(cat.get(f, ""), "#9b9a95") for f in T.feature], edgecolor=INK, lw=.4, zorder=3)
    ax.axvline(con.null_center(), color=INK2, lw=.8, ls=":")
    lo, hi = ax.get_xlim(); span = hi - lo
    if con.kind == "unpaired":
        lo, hi, span = 0, 1, 1
    for yi, r in zip(y, T.itertuples()):
        s = f"p {r.p:.2g} · q {r.q:.2f}" + (f" · {r.n_up}/{r.n} up" if con.kind == "paired" else "")
        if band is not None:
            s += f" · random-set p {band.set_index('feature').random_p[r.feature]:.2g}"
        ax.text(hi + .03 * span, yi, s, va="center", fontsize=6.8, color=INK2)
    ax.set_xlim(lo - .02 * span, hi + .62 * span)
    if con.kind == "unpaired":
        ax.set_xticks([0, .25, .5, .75, 1])
    else:
        ax.set_xticks([t for t in ax.get_xticks() if lo - .02 * span <= t <= hi])
    ax.set_yticks(y); ax.set_yticklabels([nice(f) for f in T.feature], fontsize=8.5); ax.set_xlabel(con.effect_name(), fontsize=8.5)
    cats = list(dict.fromkeys(cat.get(f, "") for f in T.feature))
    h = [Patch(color=CAT_COLORS.get(c, "#9b9a95"), label=c) for c in cats if c] + ([Patch(color=GRID, label="random GC-matched sets (95%)")] if band is not None else [])
    ax.legend(handles=h, fontsize=7, loc="upper center", bbox_to_anchor=(.4, -.08 if len(T) > 12 else -.16), ncol=4)
    ax.set_title(f"{title}\n{con.n_text()}; smallest possible p {con.min_p:.2g}", fontsize=9.5)
    _save(fig, out, f"05_effect_{_safe(con.label())}.png")


def fig_strip(st: Style, v, feat, nice, con, r, out, ylab):
    plt = _setup(); fig, ax = plt.subplots(figsize=(1.4 * len(st.groups) + 2.4, 4.6)); v = np.asarray(v, float)
    xs = {}
    rng = np.random.default_rng(0)
    for k, g in enumerate(st.groups):
        ii = np.where(st.g == g)[0]; jit = rng.uniform(-.12, .12, len(ii))
        for i, j in zip(ii, jit):
            xs[i] = k + j
        if len(ii):
            ax.hlines(np.median(v[ii]), k - .28, k + .28, color=INK, lw=1.8, zorder=4)
    if st.repeated:
        order = {g: k for k, g in enumerate(st.groups)}
        for p in st.repeated:
            ii = sorted(np.where(st.patient == p)[0], key=lambda i: order.get(st.g[i], 99))
            ax.plot([xs[i] for i in ii], v[ii], color=st.pcolor[p], lw=1, alpha=.7, zorder=2)
    for i, x in xs.items():
        ax.scatter(x, v[i], s=44, marker=st.gmark[st.g[i]], color=st.color(i), edgecolor=INK, lw=.4, zorder=3)
        if not st.repeated:
            ax.annotate(st.short[i], (x, v[i]), xytext=(5, 0), textcoords="offset points", fontsize=6.5, va="center", color=INK2)
    ax.set_xticks(range(len(st.groups))); ax.set_xticklabels([f"{g}\n(n={(st.g == g).sum()})" for g in st.groups])
    ax.set_ylabel(ylab); ax.grid(axis="y", color=GRID, lw=.5); ax.set_axisbelow(True)
    s = f"{con.label()}: p {r['p']:.2g}" + (f" ({r['n_up']}/{r['n']} up)" if con.kind == "paired" else f", AUC {r['effect']:.2f}")
    ax.set_title(f"{nice(feat)}\n{s}", fontsize=9.5)
    _save(fig, out, f"06_strip_{_safe(con.label())}_{_safe(feat)}.png")


def fig_pca_change(st: Style, S: pd.DataFrame, con, out, title):
    if con.kind != "paired":
        return None
    D = np.vstack([S.iloc[a].to_numpy(float) - S.iloc[b].to_numpy(float) for a, b, _ in con.pairs])
    P, Vt, ev = pca(pd.DataFrame(D, columns=S.columns))
    plt = _setup(); fig, ax = plt.subplots(figsize=(6.4, 5.2))
    for k, (a, b, p) in enumerate(con.pairs):
        ax.scatter(P[k, 0], P[k, 1], s=70, color=st.pcolor.get(p, PALETTE[0]), edgecolor=INK, lw=.5, zorder=3)
        ax.annotate(st._pshort(p) + (f" {st.label[a]}" if st.label[a] else ""), (P[k, 0], P[k, 1]), xytext=(5, 3), textcoords="offset points", fontsize=7, color=INK2)
    ax.axhline(0, color=GRID, lw=.8); ax.axvline(0, color=GRID, lw=.8)
    ax.set_xlabel(f"PC1 ({ev[0] * 100:.0f}%)"); ax.set_ylabel(f"PC2 ({ev[1] * 100:.0f}%)")
    ax.set_title(f"{title}\none dot per patient = that patient's change {con.case} − {con.ref}", fontsize=9)
    _save(fig, out, f"07_pca_change_{_safe(con.label())}.png")
    return P, Vt, ev


def fig_overlap(J: pd.DataFrame, nice, out):
    plt = _setup(); n = len(J); fig, ax = plt.subplots(figsize=(.28 * n + 2.5, .28 * n + 2))
    im = ax.imshow(J.to_numpy(), cmap="Greys", vmin=0, vmax=1); ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels([nice(c) for c in J.columns], rotation=90, fontsize=7); ax.set_yticklabels([nice(c) for c in J.index], fontsize=7)
    fig.colorbar(im, ax=ax, fraction=.04, pad=.02).set_label("Jaccard (shared marker bins)", fontsize=8)
    ax.set_title("Marker-set overlap: why related types move together", fontsize=9.5)
    _save(fig, out, "08_marker_overlap.png")


def lineage_scores(S: pd.DataFrame, cat: dict, split_tcga=True) -> tuple[pd.DataFrame, dict]:
    """Mean of each lineage's marker-set scores; glioma TCGA sets keep their own rows, other TCGA sets pool into one row."""
    groups = {}
    for f in S.columns:
        own = split_tcga and f.startswith("TCGA_") and cat.get(f) in ("tumour (TCGA)", "glioma (TCGA)")
        k = f if own else cat.get(f, "other")
        groups.setdefault(k, []).append(f)
    L = pd.DataFrame({k: S[v].mean(axis=1) for k, v in groups.items()})
    return L, {k: (cat.get(v[0], "other") if k.startswith("TCGA_") else k) for k, v in groups.items()}
