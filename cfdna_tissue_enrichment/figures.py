"""Figures (matplotlib, Agg backend)."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLS = ["#3a78c2", "#d9622b", "#7a5195", "#2a9d8f", "#888888"]


def _show(Z, primary, n_extra=12):
    p = [s for s in primary if s in Z.columns]
    rest = [s for s in Z.columns if s not in p]
    var = Z[rest].std().sort_values(ascending=False).index[:max(0, n_extra)] if rest else []
    return p + list(var)


def heatmap(path, Z, lab, primary, title, cbar):
    cols = _show(Z, primary); order = (list(lab.sort_values().index) + [s for s in Z.index if s not in lab.index]) if lab is not None else list(Z.index)
    M = Z.loc[order, cols].T.to_numpy(float); M = (M - np.nanmean(M, 1, keepdims=True)) / np.nanstd(M, 1, keepdims=True)
    fig, ax = plt.subplots(figsize=(max(6, .3 * len(order) + 4), .36 * len(cols) + 2.6), facecolor="white")
    im = ax.imshow(np.clip(M, -3, 3), cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto"); ax.set_yticks(range(len(cols))); ax.set_yticklabels(cols, fontsize=8)
    ax.set_xticks(range(len(order))); ax.set_xticklabels([f"{s}\n{lab.get(s, '')}" if lab is not None else s for s in order], rotation=90, fontsize=6.5)
    fig.colorbar(im, ax=ax, fraction=.03).set_label(cbar, fontsize=8); ax.set_title(title + "\n(primary tracks first, then the most variable)", fontsize=9.5)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def by_group(path, scores, lab, ref, primary, title):
    ps = [s for s in primary if s in scores.columns] or list(scores.std().sort_values(ascending=False).index[:8])
    lv = [ref] + sorted(set(lab) - {ref})
    fig, axs = plt.subplots(1, len(ps), figsize=(2.4 * len(ps), 3.9), facecolor="white", squeeze=False)
    for ax, s in zip(axs[0], ps):
        vals = [scores.loc[lab.index[lab == l], s].dropna() for l in lv]; ax.boxplot(vals, widths=.55, showfliers=False)
        for k, v in enumerate(vals): ax.scatter(np.full(len(v), k + 1) + np.random.default_rng(k).uniform(-.13, .13, len(v)), v, s=13, color=COLS[k % 5], zorder=3)
        ax.set_xticks(range(1, len(lv) + 1)); ax.set_xticklabels(lv, rotation=45, ha="right", fontsize=7); ax.set_title(s, fontsize=7.5)
    fig.suptitle(title, fontsize=9); fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def differential(path, sf, primary):
    eff = list(dict.fromkeys(sf.effect)); fig, axs = plt.subplots(1, len(eff), figsize=(5.4 * len(eff), 6.8), facecolor="white", squeeze=False)
    for ax, e in zip(axs[0], eff):
        d = sf[sf.effect == e]; p = d[d.track.isin(primary)]; o = d[~d.track.isin(primary)]
        d = pd.concat([p, o.reindex(o.rho.abs().sort_values(ascending=False).index).head(14)]).iloc[::-1]
        ax.barh(range(len(d)), d.rho, color=np.where(d.rho > 0, "#d9622b", "#3a78c2"))
        ax.set_yticks(range(len(d))); ax.set_yticklabels([f"{t}{' *' if p_ < .05 else ''}" for t, p_ in zip(d.track, d.p_perm)], fontsize=7)
        ax.axvline(0, color="#999", lw=.6); ax.set_xlabel("Spearman rho (bin t vs track specificity)", fontsize=8); ax.set_title(e, fontsize=9)
    fig.suptitle("Differential, all bins: + = bins specific to the track have higher coverage in cases (* permutation p < 0.05)", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def qc(path, Q, hists, lmax):
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.6), facecolor="white"); L = np.arange(lmax + 1)
    for h in hists.values(): axs[0].plot(L, h / max(h.sum(), 1), lw=.8)
    axs[0].set_xlim(50, lmax); axs[0].set_xlabel("fragment length (bp)"); axs[0].set_ylabel("fraction"); axs[0].set_title("fragment length", fontsize=9)
    axs[1].bar(range(len(Q)), Q.fragments_used / 1e6, color=np.where(Q.flag_low_depth, "#d9622b", "#888")); axs[1].set_xticks(range(len(Q)))
    axs[1].set_xticklabels(Q.index, rotation=90, fontsize=7); axs[1].set_ylabel("fragments used (M)"); axs[1].set_title("depth (orange = below minimum)", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


# ------------------------------------------------------------------ ranked bars (all tracks, two panels)
IMMUNE = ("B cell", "monocyte", "T cell", "natural killer", "neutrophil", "myeloid progenitor", "spleen", "thymus", "Peyer's patch", "germinal center")
KIND_COL = {"immune / blood": "#c0392b", "tissue / other cell": "#7f8c8d", "TCGA tumour": "#2471a3", "maternal-fetal": "#8e44ad", "derived (group / axis)": "#d4a017"}


def track_kind(t):
    if t.startswith("C_"): return "immune / blood" if any(k in t for k in IMMUNE) else "tissue / other cell"
    if t.startswith("T_"): return "TCGA tumour"
    if t.startswith("W_"): return "maternal-fetal"
    return "derived (group / axis)"


def track_label(t):
    return ("TCGA " + t[2:]) if t.startswith("T_") else t[2:] if t[:2] in ("C_", "W_") else t


def ranked_bars(path, d, val, p, title, xlabel, case="cases", ref="reference", top=None, q=None):
    """Two panels like the original cfdna-chromatin-xai plot: tracks lower in `case` (left) and higher (right), ranked by
    |val|, coloured by track type. * permutation p < 0.05, ** also q < 0.10. top=None shows every track."""
    d = d[np.isfinite(d[val])]
    lo = d[d[val] < 0].sort_values(val); hi = d[d[val] > 0].sort_values(val, ascending=False)
    if top: lo, hi = lo.head(top), hi.head(top)
    n = max(len(lo), len(hi), 1); xmax = float(np.nanmax(np.abs(d[val]))) * 1.15 if len(d) else 1
    fig, axs = plt.subplots(1, 2, figsize=(13, .19 * n + 2.0), facecolor="white")
    for ax, part, head in ((axs[0], lo, f"Lower coverage in {case}"), (axs[1], hi, f"Higher coverage in {case}")):
        part = part.iloc[::-1]; y = np.arange(len(part)) + (n - len(part))
        ax.barh(y, part[val].abs(), color=[KIND_COL[track_kind(t)] for t in part.track], height=.75)
        star = ["**" if (q is not None and qq < .10 and pp < .05) else "*" if pp < .05 else "" for pp, qq in
                zip(part[p], part[q] if q is not None else part[p])]
        ax.set_yticks(y); ax.set_yticklabels([f"{track_label(t)} {s}".strip() for t, s in zip(part.track, star)], fontsize=6.5)
        ax.set_ylim(-.7, n - .3); ax.set_xlim(0, xmax); ax.set_xlabel(f"|{xlabel}|", fontsize=8); ax.set_title(head, fontsize=10)
        ax.tick_params(axis="x", labelsize=7)
    hs = [plt.Rectangle((0, 0), 1, 1, color=c) for k, c in KIND_COL.items() if (d.track.map(track_kind) == k).any()]
    ks = [k for k in KIND_COL if (d.track.map(track_kind) == k).any()]
    fig.legend(hs, ks, loc="lower center", ncol=len(ks), fontsize=7.5, frameon=False)
    fig.suptitle(f"{title}\nreference: {ref}   * permutation p < 0.05" + ("   ** and q < 0.10" if q is not None else ""), fontsize=10)
    fig.tight_layout(rect=(0, .035, 1, 1)); fig.savefig(path, dpi=150); plt.close(fig)


def effect_parts(effect):
    """'col:case_vs_ref' -> (case, ref); falls back to ('cases', 'reference')."""
    e = effect.split(":", 1)[-1]
    return tuple(e.split("_vs_", 1)) if "_vs_" in e and " x " not in e else ("cases", "reference")


def ranked_outputs(out, top=None):
    """Draw ranked bar plots for every effect of a finished run folder; returns the files written."""
    import os, re
    done = []
    specs = [("sample", "effect_tests.csv", "score", "t", "Sample mode: group effect on E", "t (case vs reference)"),
             ("differential", "enrichment_all_bins.csv", "track", "rho", "Differential mode, all bins", "Spearman rho (bin t vs specificity)")]
    for mode, f, tcol, val, title, xl in specs:
        fp = os.path.join(out, mode, f)
        if not os.path.exists(fp): continue
        T = pd.read_csv(fp).rename(columns={tcol: "track"})
        for e, d in T.groupby("effect", sort=False):
            case, ref = effect_parts(e); fn = os.path.join(out, mode, f"ranked_{re.sub(r'[^A-Za-z0-9]+', '_', e).strip('_')}.png")
            ranked_bars(fn, d, val, "p_perm", f"{title} — {e}", xl, case, ref, top, "q_perm" if "q_perm" in d else None); done.append(fn)
    return done


def top_bins_activity(path, tp, loo=None, case="cases", ref="reference", top=20, primary=()):
    """Two panels from enrichment_top_bins.csv: tracks of the top bins with MORE cfDNA coverage in cases ('up') and in
    the reference ('down'), x = closed-form z (descriptive), small dots = leave-one-out runs when given. cfDNA logic:
    open / active chromatin loses coverage, so 'more coverage in cases' = cell types more active in the reference."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(15, 0.36 * min(top, tp.track.nunique()) + 2.2))
    for ax, (direction, g1, g2) in zip(axs, (("up", case, ref), ("down", ref, case))):
        d = tp[tp.direction == direction].copy()
        d = d[d.track.isin(primary) | d.z_closed_form.abs().rank(ascending=False).le(top)].sort_values("z_closed_form").tail(top)
        y = np.arange(len(d))
        if loo is not None:
            l = loo[loo.direction == direction]
            for k, t in enumerate(d.track):
                v = l.loc[l.track == t, "z_closed_form"]
                if len(v): ax.plot([v.min(), v.max()], [k, k], color="#bbbbbb", lw=3, zorder=1); ax.scatter(v, [k] * len(v), s=8, color="#777777", zorder=2)
        ax.scatter(d.z_closed_form, y, s=60, color=[{"T": "#6a3d9a", "C": "#d1453b", "W": "#1a7f2e"}.get(str(t)[0], "#333333") for t in d.track], zorder=3)
        ax.set_yticks(y); ax.set_yticklabels([track_label(t) for t in d.track], fontsize=9)
        ax.axvline(0, color="k", lw=1); ax.axvline(1.96, color="grey", ls="--", lw=1)
        ax.set_title(f"More cfDNA coverage in {g1}\n= cell types more active in {g2}", fontsize=11); ax.set_xlabel("specificity z, closed form (top bins; descriptive)")
        for sp_ in ("top", "right"): ax.spines[sp_].set_visible(False)
    fig.tight_layout(); fig.savefig(path, dpi=200); plt.close(fig)
