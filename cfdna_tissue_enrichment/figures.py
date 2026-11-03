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
