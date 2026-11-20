"""cfdna-tissue-enrich: tissue / cell-type enrichment of plasma cfDNA coverage against chromatin-accessibility atlases.

  cfdna-tissue-enrich run    --samples S.csv --out OUT --panel NAME [--modes sample,differential,markers] [--factor COL:REF:CASE ...]
  cfdna-tissue-enrich count  --samples S.csv --cache DIR [--index i]   per-sample bin counts (parallelisable; reused by run --cache)
  cfdna-tissue-enrich panels                       list panels and their layers
  cfdna-tissue-enrich tracks --panel NAME          list a panel's tracks (and which are primary)
  cfdna-tissue-enrich plot   --out RUN [--top N]   (re)draw ranked bar plots of all tracks for a finished run

modes  sample        per-sample enrichment: bin z vs a reference group (or all samples), Spearman with track specificity
       differential  per-bin group effect, then its enrichment in each track
       markers       per-sample mean coverage at cell-type marker bins (10 kb marker atlas) vs GC-matched random sets
sample and markers also give the per-sample report (contrasts: paired sign-flip or unpaired AUC; PCA; heatmaps) when the
samples file has the --group-col column. Samples: frag_path (fragment file) or coverage (per-sample CPM csv, key + *_cpm).
"""
import argparse, json, os, re, sys, time
import numpy as np, pandas as pd
from . import __version__
from . import reference as RF, counts as CN, stats as ST, modes as MD, figures as FG, report as RP, markers as MK


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def get_counts(sid, path, a, g, lut, bs):
    """Per-sample bin counts, read from / written to --cache when given (keyed by sample, build, bin and filters)."""
    tag = f"{a.build}.{a.bin}.q{a.min_mapq}.l{a.min_len}-{a.max_len}.m{a.mapq_col}"
    cf = os.path.join(a.cache, f"{sid}.{tag}.npz") if getattr(a, "cache", None) else None
    if cf and os.path.exists(cf):
        z = np.load(cf, allow_pickle=False)
        if len(z["counts"]) == len(g): return z["counts"], z["hist"], int(z["n_all"]), int(z["n_pass"]), bool(z["has_q"]), True
    cts, hist, n_all, n_pass, has_q = CN.count_fragments(path, lut, bs, len(g), a.min_mapq, a.min_len, a.max_len, a.mapq_col)
    if cf:
        os.makedirs(a.cache, exist_ok=True); tmp = cf + ".tmp.npz"
        np.savez_compressed(tmp, counts=cts.astype(np.int32), hist=hist, n_all=n_all, n_pass=n_pass, has_q=has_q, frag_path=path); os.replace(tmp, cf)
    return cts, hist, n_all, n_pass, has_q, False


def read_samples(path):
    meta = pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[""]); meta.columns = [c.strip() for c in meta.columns]
    if "sample_id" not in meta: sys.exit("samples file needs a 'sample_id' column")
    if "frag_path" not in meta and "coverage" not in meta: sys.exit("samples file needs a 'frag_path' (fragment file) or 'coverage' (CPM csv) column")
    if meta.sample_id.duplicated().any(): sys.exit("duplicate sample_id values")
    base = os.path.dirname(os.path.abspath(path))
    for c in ("frag_path", "coverage"):
        if c in meta: meta[c] = [p if not isinstance(p, str) or os.path.isabs(p) else os.path.join(base, p) for p in meta[c]]
    if "coverage" not in meta: meta["coverage"] = np.nan
    if "frag_path" not in meta: meta["frag_path"] = np.nan
    meta["input"] = [c if isinstance(c, str) else f for f, c in zip(meta.frag_path, meta.coverage)]
    if meta.input.isna().any(): sys.exit(f"samples without frag_path or coverage: {list(meta.sample_id[meta.input.isna()])[:5]}")
    return meta


def cmd_count(a):
    """Count fragments per bin for one sample (--index, 0-based row; for job arrays) or all samples, into --cache."""
    meta = read_samples(a.samples); rdir = RF.ref_dir(a.build, a.reference); bs = RF.BIN[a.bin]; g = RF.load_grid(rdir, a.bin); lut = RF.lookup_tables(g, bs)
    rows = [meta.iloc[a.index]] if a.index is not None else [r for _, r in meta.iterrows()]
    for r in rows:
        if not isinstance(r.frag_path, str): continue
        if not os.path.exists(r.frag_path): sys.exit(f"missing {r.frag_path}")
        cts, hist, n_all, n_pass, has_q, cached = get_counts(r.sample_id, r.frag_path, a, g, lut, bs)
        log(f"{r.sample_id}: {n_pass:,} fragments used of {n_all:,} (MAPQ filter {'on' if has_q else 'off'}){' [cached]' if cached else ''}")


def cmd_plot(a):
    for f in FG.ranked_outputs(a.out, a.top): log("wrote", f)


def cmd_panels(a):
    for p in RF.list_panels():
        P = RF.load_panel(p); print(f"{p:<10} layers {','.join(P['layers']):<8} {P.get('description', '')}")


def cmd_tracks(a):
    P = RF.load_panel(a.panel); rdir = RF.ref_dir(a.build, a.reference); g = RF.load_grid(rdir, a.bin)
    layers, meta = RF.load_layers(rdir, a.bin, P["layers"], g.key.to_numpy()); SP = RF.specificity(layers, P, a.specificity or P.get("specificity", "layer")); prim = RF.primary_tracks(P)
    for c in SP.columns: print(f"{'*' if c in prim else ' '} {c}")
    print(f"\n{SP.shape[1]} tracks; * = primary ({len(prim)})")


def cmd_run(a):
    os.makedirs(a.out, exist_ok=True); modes = [m.strip() for m in a.modes.split(",") if m.strip()]
    bad = set(modes) - {"sample", "differential", "markers"}
    if bad: sys.exit(f"unknown modes: {bad}")
    if "markers" in modes and not a.markers: sys.exit("markers mode needs --markers (marker atlas csv.gz: cell_type, key, rank)")
    meta = read_samples(a.samples)
    gcol = a.group_col or a.group or "group"
    if a.groups:
        if gcol not in meta: sys.exit(f"--groups: column '{gcol}' not in samples file")
        lv = a.groups.split(","); meta = meta[meta[gcol].isin(lv)].reset_index(drop=True)
    miss = [p for p in meta.input if not os.path.exists(p)]
    if miss: sys.exit(f"missing input files: {miss[:5]}")
    if a.group:
        lv = sorted(meta[a.group].dropna().astype(str).unique())
        if len(lv) < 2: sys.exit(f"--group {a.group} needs two levels")
        ref = a.reference_level or [l for l in lv if l != (a.case or lv[-1])][0]; case = a.case or [l for l in lv if l != ref][0]
        a.factor.append(f"{a.group}:{ref}:{case}")
    F = ST.parse_factors(a.factor, meta); covs = [c for c in a.covariates.split(",") if c]
    badc = [c.split(":")[0] for c in covs if c.split(":")[0] not in meta]
    if badc: sys.exit(f"covariates not in samples file: {badc}")
    pooled = a.diff_stat == "pooled" and "differential" in modes
    if pooled and (len(F) != 1 or covs or a.interaction):
        sys.exit("--diff-stat pooled needs exactly one two-level factor and no covariates / interaction (samples are pooled per group)")
    if not F and ("differential" in modes or ("sample" in modes and a.sample_reference == "group")):
        sys.exit("no groups given: differential mode needs --factor (or --group), and so does sample mode unless --sample-reference all")
    meta = meta.set_index("sample_id", drop=False)

    P = RF.load_panel(a.panel); rdir = RF.ref_dir(a.build, a.reference); bs = RF.BIN[a.bin]
    g = RF.load_grid(rdir, a.bin); lut = RF.lookup_tables(g, bs); spec_m = a.specificity or P.get("specificity", "layer"); SP = None; primary = []
    if "sample" in modes or "differential" in modes:
        layers, lmeta = RF.load_layers(rdir, a.bin, P["layers"], g.key.to_numpy()); SP = RF.specificity(layers, P, spec_m)
        primary = [t for t in RF.primary_tracks(P, a.primary) if t in SP.columns]
        log(f"panel '{a.panel}' ({', '.join(P['layers'])}), {a.build} {a.bin}: {len(g):,} bins, {SP.shape[1]} tracks ({len(primary)} primary); modes: {', '.join(modes)}")
    else:
        lmeta = {"layers": {}}; log(f"{a.build} {a.bin}: {len(g):,} bins; modes: {', '.join(modes)}")
    gcc = a.gc_correct == "yes"; win = a.detrend_bins or CN.DETREND[a.bin]
    log(f"preprocessing: GC correction {'on' if gcc else 'OFF'}, detrend {win} bins; specificity '{spec_m}'; differential statistic '{a.diff_stat}'")

    samples = list(meta.sample_id); Xp = np.zeros((len(samples), len(g))); Lev = np.zeros_like(Xp) if pooled else None; qc, hists = [], {}
    for i, sid in enumerate(samples):
        log(f"sample {i + 1}/{len(samples)} {sid}")
        if isinstance(meta.loc[sid, "frag_path"], str):
            cts, hist, n_all, n_pass, has_q, _ = get_counts(sid, meta.loc[sid, "frag_path"], a, g, lut, bs); cpm = None
            hists[sid] = hist; qc.append(CN.qc_row(sid, cts, hist, n_all, n_pass, has_q, a.min_fragments))
        else:
            cts = None; cpm = CN.read_coverage(meta.loc[sid, "coverage"], g.key)
        if pooled: Xp[i], Lev[i] = CN.preprocess(cts, g, a.bin, gcc, win, level=True, cpm=cpm)
        else: Xp[i] = CN.preprocess(cts, g, a.bin, gcc, win, cpm=cpm)
    if qc:
        Q = pd.DataFrame(qc).set_index("sample_id"); Q.to_csv(os.path.join(a.out, "qc.csv")); FG.qc(os.path.join(a.out, "qc.png"), Q, hists, a.max_len)
    rep_opts = None
    if gcol in meta:
        lv = a.groups.split(",") if a.groups else list(dict.fromkeys(meta[gcol].dropna().astype(str)))
        pcol = None if a.patient_col == "none" else a.patient_col
        rep_opts = argparse.Namespace(group_col=gcol, groups=lv, patient_col=pcol, label_col=a.label_col, contrasts=a.contrasts or (
            ",".join(f"{c}:{r}" for col, r, c in F if col == gcol) or None), connect=a.connect, top_features=a.top_features, pca_centre=None)
        rmeta = meta.reset_index(drop=True)

    lab = refcell = Xd = terms = None
    if F:
        Xd, terms, mm = ST.design_matrix(meta, F, covs, a.interaction); lab, refcell = ST.cells(mm, F)
        n_ref = int((lab == refcell).sum()); log(f"design: {len(Xd)} samples, effects {terms}, reference group '{refcell}' (n={n_ref}); groups {lab.value_counts().to_dict()}")
        if n_ref < 2: sys.exit(f"reference group '{refcell}' has {n_ref} samples; need >= 2")
        Xd.to_csv(os.path.join(a.out, "design_matrix.csv"))
    info = {}

    if "sample" in modes:
        od = os.path.join(a.out, "sample"); os.makedirs(od, exist_ok=True)
        if a.sample_reference == "all":
            ref_mask = np.ones(len(samples), bool); E, nb, loo = MD.mode_sample(Xp, samples, SP, ref_mask, loo_min=None); refname = "all samples"
        else:
            ref_mask = np.array([s in lab.index and lab[s] == refcell for s in samples]); refname = refcell
            E, nb, loo = MD.mode_sample(Xp, samples, SP, ref_mask)
            if not loo: log(f"note: {int(ref_mask.sum())} reference samples; they are not left out of their own reference (needs >= 4)")
        E.to_csv(os.path.join(od, "scores_E.csv"))
        if F:
            FG.heatmap(os.path.join(od, "scores_heatmap.png"), E, lab, primary, f"Sample: E vs reference '{refname}'", "row z of E")
            ST.effect_tests(E, Xd, terms, primary, a.n_perm, a.seed).to_csv(os.path.join(od, "effect_tests.csv"), index=False)
            ST.pairwise_tests(E.loc[lab.index], lab, refcell, primary, a.n_perm, a.seed).to_csv(os.path.join(od, "pairwise_vs_reference.csv"), index=False)
            FG.by_group(os.path.join(od, "primary_by_group.png"), E.loc[lab.index], lab, refcell, primary, f"Sample: E vs reference '{refname}'")
        if rep_opts:
            cat = {t: RP.track_category(t) for t in E.columns}; Sx = E[RP.order(list(E.columns), cat, RP.TISSUE_ORDER)].reset_index(drop=True)
            rep_opts.pca_centre = a.pca_centre or "no"
            RP.run(Sx, rmeta, rep_opts, "sample", cat, RP.track_nice, "per-sample enrichment (Spearman ρ)", os.path.join(od, "report"), log=log)
        info["sample"] = dict(bins_used=nb, reference=refname, n_reference=int(ref_mask.sum()), leave_one_out=loo); log("sample done")

    if "markers" in modes:
        od = os.path.join(a.out, "markers"); os.makedirs(od, exist_ok=True)
        Sm, rand, J, nbin = MD.mode_markers(Xp, g, a.markers, a.marker_top, a.n_random, a.include, a.exclude, a.seed)
        cat = RP.markers_categories(Sm.columns); Sm = Sm[RP.order(list(Sm.columns), cat, RP.MARKER_ORDER)]
        Sm.set_index(pd.Index(samples, name="sample_id")).to_csv(os.path.join(od, "scores.csv"))
        if rand:
            pd.DataFrame({c: rand[c].mean(1) for c in Sm.columns}, index=pd.Index(samples, name="sample_id")).to_csv(os.path.join(od, "random_sets_mean.csv"))
        if rep_opts:
            rep_opts.pca_centre = a.pca_centre or "both"
            RP.run(Sm, rmeta, rep_opts, "markers", cat, MK.nice, f"mean preprocessed coverage (top {a.marker_top} bins)", os.path.join(od, "report"), rand, J, log=log)
        info["markers"] = dict(atlas=a.markers, top=a.marker_top, n_sets=Sm.shape[1], bins_per_set_min=int(nbin.min()), n_random=a.n_random); log("markers done")

    if "differential" in modes:
        od = os.path.join(a.out, "differential"); os.makedirs(od, exist_ok=True)
        lo = None
        if pooled:
            sf, tp, bt, nb, lo = MD.mode_differential_pooled(Lev, samples, SP, Xd, terms[0], primary, g, a.n_perm_bins, a.top_n, a.n_random, a.seed, win, a.loo)
            info["differential"] = dict(bins_used=nb, statistic="pooled", n_case=int((Xd[terms[0]] == 1).sum()), n_reference=int((Xd[terms[0]] == 0).sum()))
            if lo is not None: lo.to_csv(os.path.join(od, "enrichment_top_bins_loo.csv"), index=False)
        else:
            sf, tp, bt, nb, df_ = MD.mode_differential(Xp, samples, SP, Xd, terms, primary, g, a.n_perm_bins, a.top_n, a.n_random, a.seed)
            info["differential"] = dict(bins_used=nb, statistic="moderated t", residual_df=float(df_))
            if df_ < 8: log(f"note: differential mode has {df_:.0f} residual df; bin-level t are unstable")
        sf.to_csv(os.path.join(od, "enrichment_all_bins.csv"), index=False); tp.to_csv(os.path.join(od, "enrichment_top_bins.csv"), index=False)
        bt.to_csv(os.path.join(od, "bin_effects.csv.gz"), index=False, float_format="%.5g"); FG.differential(os.path.join(od, "enrichment_all_bins.png"), sf, primary)
        for e_, t_ in tp.groupby("effect", sort=False):
            case_, ref_ = FG.effect_parts(e_); fn = os.path.join(od, f"top_bins_activity_{re.sub(r'[^A-Za-z0-9]+', '_', e_).strip('_')}.png")
            FG.top_bins_activity(fn, t_, None if lo is None else lo[lo.effect == e_], case_, ref_, a.bar_top or 20, primary)
        log("differential done")

    if F: FG.ranked_outputs(a.out, a.bar_top)
    json.dump(dict(version=__version__, panel=a.panel, panel_definition=P, layers={k: lmeta["layers"][k] for k in P["layers"] if k in lmeta["layers"]}, build=a.build, bin=a.bin,
                   modes=modes, factors=F, reference_group=refcell, effects=terms, covariates=covs, primary=primary, mode_info=info,
                   args=vars(a), n_samples=len(samples), time=time.strftime("%Y-%m-%d %H:%M:%S")),
              open(os.path.join(a.out, "run_manifest.json"), "w"), indent=1, default=str)
    log("done ->", a.out)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    ap = argparse.ArgumentParser(prog="cfdna-tissue-enrich", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "tracks", "count"):
        p = sub.add_parser(name)
        p.add_argument("--panel", default="pan", help="panel name (see 'panels') or a panel JSON path")
        p.add_argument("--build", default="hg38"); p.add_argument("--bin", default="50kb", choices=list(RF.BIN))
        p.add_argument("--reference", default=None, help="reference root (default: ./reference or $CFTE_REFERENCE)")
        p.add_argument("--specificity", default=None, choices=["layer", "combined", "combined-z"],
                       help="track specificity within each layer, across all panel layers (combined: cfdna-chromatin-xai; combined-z: C1D1 / Ahmed); default: the panel's, else layer")
        if name in ("run", "count"):
            p.add_argument("--samples", required=True, help="CSV with sample_id and frag_path or coverage (+ metadata columns)")
            p.add_argument("--cache", default=None, help="folder for per-sample bin counts; reused by later runs")
            p.add_argument("--min-mapq", type=int, default=30); p.add_argument("--mapq-col", type=int, default=None, help="1-based MAPQ column (default: header 'mapq', else column 4 or 5 if numeric; 0 = none)")
            p.add_argument("--min-len", type=int, default=50); p.add_argument("--max-len", type=int, default=500)
        if name == "count":
            p.add_argument("--index", type=int, default=None, help="0-based row of the samples file (for Slurm arrays)")
        if name == "run":
            p.add_argument("--out", required=True)
            p.add_argument("--modes", default="sample,differential", help="comma list of sample, differential, markers")
            p.add_argument("--factor", action="append", default=[], help="COLUMN:REFERENCE:CASE, repeatable")
            p.add_argument("--group", help="shortcut for one factor (with --case, --reference-level)"); p.add_argument("--case"); p.add_argument("--reference-level")
            p.add_argument("--interaction", action="store_true"); p.add_argument("--covariates", default="", help="comma-separated; :num / :cat to force type")
            p.add_argument("--primary", default=None, help="primary tracks, comma-separated (use ';' if a track name contains a comma); overrides the panel's list")
            p.add_argument("--min-fragments", type=float, default=5e6)
            p.add_argument("--n-random", type=int, default=200); p.add_argument("--n-perm", type=int, default=2000); p.add_argument("--n-perm-bins", type=int, default=200)
            p.add_argument("--top-n", type=int, default=2000); p.add_argument("--seed", type=int, default=17)
            p.add_argument("--bar-top", type=int, default=None, help="ranked bar plots: tracks per panel (default: all)")
            p.add_argument("--gc-correct", default="yes", choices=["yes", "no"], help="per-sample GC LOWESS before the detrend (default yes)")
            p.add_argument("--detrend-bins", type=int, default=None, help="rolling-median detrend window in bins (default ~2 Mb: 41 at 50kb, 205 at 10kb)")
            p.add_argument("--diff-stat", default="t", choices=["t", "pooled"], help="differential per-bin statistic: moderated t (default) or pooled-group difference (cfdna-chromatin-xai)")
            p.add_argument("--loo", action="store_true", help="pooled differential: leave-one-out re-pool + re-score of the top bins")
            p.add_argument("--sample-reference", default="group", choices=["group", "all"], help="sample mode: per-bin z against the reference group (default) or all samples (C1D1 / Ahmed)")
            p.add_argument("--markers", help="markers mode: marker atlas csv.gz (cell_type, key, rank; e.g. 10 kb immune / TCGA / CNS sets); use --bin 10kb")
            p.add_argument("--marker-top", type=int, default=1000, help="markers mode: top-N bins per set (default 1000)")
            p.add_argument("--include", help="markers mode: regex, keep matching sets"); p.add_argument("--exclude", help="markers mode: regex, drop matching sets")
            r = p.add_argument_group("per-sample report (sample and markers modes)")
            r.add_argument("--group-col", default=None, help="group column (default: --group, else 'group')")
            r.add_argument("--groups", help="comma list: keep only these groups, in this order")
            r.add_argument("--contrasts", help="case:ref[,case:ref] (default: the --factor/--group of the group column, else adjacent groups, later minus earlier)")
            r.add_argument("--patient-col", default="patient", help="links samples of one patient: paired sign-flip tests ('none' = unpaired)")
            r.add_argument("--label-col", help="extra per-sample label on plots (e.g. cycle)"); r.add_argument("--connect", action="store_true", help="arrows between a patient's samples in the PCA")
            r.add_argument("--top-features", type=int, default=6, help="strip plots for the N strongest features per contrast")
            r.add_argument("--pca-centre", default=None, choices=["no", "yes", "both"], help="PCA on raw and / or sample-centred scores (default: sample no, markers both)")
    sub.add_parser("panels")
    pp = sub.add_parser("plot"); pp.add_argument("--out", required=True, help="finished run folder"); pp.add_argument("--top", type=int, default=None, help="tracks per panel (default: all)")
    a = ap.parse_args(argv)
    {"run": cmd_run, "panels": cmd_panels, "tracks": cmd_tracks, "count": cmd_count, "plot": cmd_plot}[a.cmd](a)


if __name__ == "__main__":
    main()
