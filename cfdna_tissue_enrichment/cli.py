"""cfdna-tissue-enrich: tissue / cell-type enrichment of plasma cfDNA coverage against chromatin-accessibility atlases.

  cfdna-tissue-enrich run    --samples S.csv --out OUT --panel prenatal|tumor|pan [--factor COL:REF:CASE ...]
  cfdna-tissue-enrich count  --samples S.csv --cache DIR [--index i]   per-sample bin counts (parallelisable; reused by run --cache)
  cfdna-tissue-enrich panels                       list panels and their layers
  cfdna-tissue-enrich tracks --panel NAME          list a panel's tracks (and which are primary)
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from . import __version__
from . import reference as RF, counts as CN, stats as ST, modes as MD, figures as FG


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
    for c in ("sample_id", "frag_path"):
        if c not in meta: sys.exit(f"samples file needs a '{c}' column")
    if meta.sample_id.duplicated().any(): sys.exit("duplicate sample_id values")
    base = os.path.dirname(os.path.abspath(path)); meta["frag_path"] = [p if os.path.isabs(p) else os.path.join(base, p) for p in meta.frag_path]
    return meta


def cmd_count(a):
    """Count fragments per bin for one sample (--index, 0-based row; for job arrays) or all samples, into --cache."""
    meta = read_samples(a.samples); rdir = RF.ref_dir(a.build, a.reference); bs = RF.BIN[a.bin]; g = RF.load_grid(rdir, a.bin); lut = RF.lookup_tables(g, bs)
    rows = [meta.iloc[a.index]] if a.index is not None else [r for _, r in meta.iterrows()]
    for r in rows:
        if not os.path.exists(r.frag_path): sys.exit(f"missing {r.frag_path}")
        cts, hist, n_all, n_pass, has_q, cached = get_counts(r.sample_id, r.frag_path, a, g, lut, bs)
        log(f"{r.sample_id}: {n_pass:,} fragments used of {n_all:,} (MAPQ filter {'on' if has_q else 'off'}){' [cached]' if cached else ''}")


def cmd_panels(a):
    for p in RF.list_panels():
        P = RF.load_panel(p); print(f"{p:<10} layers {','.join(P['layers']):<8} {P.get('description', '')}")


def cmd_tracks(a):
    P = RF.load_panel(a.panel); rdir = RF.ref_dir(a.build, a.reference); g = RF.load_grid(rdir, a.bin)
    layers, meta = RF.load_layers(rdir, a.bin, P["layers"], g.key.to_numpy()); SP = RF.specificity(layers, P); prim = RF.primary_tracks(P)
    for c in SP.columns: print(f"{'*' if c in prim else ' '} {c}")
    print(f"\n{SP.shape[1]} tracks; * = primary ({len(prim)})")


def cmd_run(a):
    os.makedirs(a.out, exist_ok=True); modes = [m.strip() for m in a.modes.split(",") if m.strip()]
    bad = set(modes) - {"sets", "sample", "differential"}
    if bad: sys.exit(f"unknown modes: {bad}")
    meta = read_samples(a.samples)
    miss = [p for p in meta.frag_path if not os.path.exists(p)]
    if miss: sys.exit(f"missing fragment files: {miss[:5]}")
    if a.group:
        lv = sorted(meta[a.group].dropna().astype(str).unique())
        if len(lv) < 2: sys.exit(f"--group {a.group} needs two levels")
        ref = a.reference_level or [l for l in lv if l != (a.case or lv[-1])][0]; case = a.case or [l for l in lv if l != ref][0]
        a.factor.append(f"{a.group}:{ref}:{case}")
    F = ST.parse_factors(a.factor, meta); covs = [c for c in a.covariates.split(",") if c]
    badc = [c.split(":")[0] for c in covs if c.split(":")[0] not in meta]
    if badc: sys.exit(f"covariates not in samples file: {badc}")
    if not F and ("sample" in modes or "differential" in modes):
        log("no --factor given: running 'sets' only (sample mode needs a reference group, differential mode needs groups)"); modes = ["sets"]
    meta = meta.set_index("sample_id", drop=False)

    P = RF.load_panel(a.panel); rdir = RF.ref_dir(a.build, a.reference); bs = RF.BIN[a.bin]
    g = RF.load_grid(rdir, a.bin); lut = RF.lookup_tables(g, bs)
    layers, lmeta = RF.load_layers(rdir, a.bin, P["layers"], g.key.to_numpy()); SP = RF.specificity(layers, P)
    primary = [t for t in RF.primary_tracks(P, a.primary) if t in SP.columns]
    log(f"panel '{a.panel}' ({', '.join(P['layers'])}), {a.build} {a.bin}: {len(g):,} bins, {SP.shape[1]} tracks ({len(primary)} primary); modes: {', '.join(modes)}")

    samples = list(meta.sample_id); Xp = np.zeros((len(samples), len(g))); qc, hists = [], {}
    for i, sid in enumerate(samples):
        log(f"sample {i + 1}/{len(samples)} {sid}")
        cts, hist, n_all, n_pass, has_q, _ = get_counts(sid, meta.loc[sid, "frag_path"], a, g, lut, bs)
        Xp[i] = CN.preprocess(cts, g, a.bin); hists[sid] = hist; qc.append(CN.qc_row(sid, cts, hist, n_all, n_pass, has_q, a.min_fragments))
    Q = pd.DataFrame(qc).set_index("sample_id"); Q.to_csv(os.path.join(a.out, "qc.csv")); FG.qc(os.path.join(a.out, "qc.png"), Q, hists, a.max_len)

    lab = refcell = Xd = terms = None
    if F:
        Xd, terms, mm = ST.design_matrix(meta, F, covs, a.interaction); lab, refcell = ST.cells(mm, F)
        n_ref = int((lab == refcell).sum()); log(f"design: {len(Xd)} samples, effects {terms}, reference group '{refcell}' (n={n_ref}); groups {lab.value_counts().to_dict()}")
        if n_ref < 2: sys.exit(f"reference group '{refcell}' has {n_ref} samples; need >= 2")
        Xd.to_csv(os.path.join(a.out, "design_matrix.csv"))
    info = {}

    if "sets" in modes:
        od = os.path.join(a.out, "sets"); os.makedirs(od, exist_ok=True); top = a.set_size or RF.TOPN[a.bin]
        sets = RF.region_sets(SP, top); (Rw, Cw, Zw), sinfo = MD.mode_sets(Xp, samples, sets, g, a.n_random, a.seed)
        Rw.to_csv(os.path.join(od, "scores_raw.csv")); Cw.to_csv(os.path.join(od, "scores_corrected.csv")); Zw.to_csv(os.path.join(od, "scores_z.csv"))
        sinfo.assign(primary=sinfo.set.isin(primary)).to_csv(os.path.join(od, "region_sets_used.csv"), index=False)
        FG.heatmap(os.path.join(od, "scores_heatmap.png"), Zw, lab, primary, "Sets: score z vs matched random bins", "row z")
        if F:
            ST.effect_tests(Cw, Xd, terms, primary, a.n_perm, a.seed).to_csv(os.path.join(od, "effect_tests.csv"), index=False)
            ST.pairwise_tests(Cw.loc[lab.index], lab, refcell, primary, a.n_perm, a.seed).to_csv(os.path.join(od, "pairwise_vs_reference.csv"), index=False)
            FG.by_group(os.path.join(od, "primary_by_group.png"), Cw.loc[lab.index], lab, refcell, primary, "Sets: corrected scores")
        info["sets"] = dict(bins_per_set=top); log("sets done")

    if "sample" in modes:
        od = os.path.join(a.out, "sample"); os.makedirs(od, exist_ok=True)
        ref_mask = np.array([s in lab.index and lab[s] == refcell for s in samples])
        E, nb, loo = MD.mode_sample(Xp, samples, SP, ref_mask); E.to_csv(os.path.join(od, "scores_E.csv"))
        if not loo: log(f"note: {int(ref_mask.sum())} reference samples; they are not left out of their own reference (needs >= 4)")
        FG.heatmap(os.path.join(od, "scores_heatmap.png"), E, lab, primary, f"Sample: E vs reference '{refcell}'", "row z of E")
        ST.effect_tests(E, Xd, terms, primary, a.n_perm, a.seed).to_csv(os.path.join(od, "effect_tests.csv"), index=False)
        ST.pairwise_tests(E.loc[lab.index], lab, refcell, primary, a.n_perm, a.seed).to_csv(os.path.join(od, "pairwise_vs_reference.csv"), index=False)
        FG.by_group(os.path.join(od, "primary_by_group.png"), E.loc[lab.index], lab, refcell, primary, f"Sample: E vs reference '{refcell}'")
        info["sample"] = dict(bins_used=nb, reference_group=refcell, n_reference=int(ref_mask.sum()), leave_one_out=loo); log("sample done")

    if "differential" in modes:
        od = os.path.join(a.out, "differential"); os.makedirs(od, exist_ok=True)
        sf, tp, bt, nb, df_ = MD.mode_differential(Xp, samples, SP, Xd, terms, primary, g, a.n_perm_bins, a.top_n, a.n_random, a.seed)
        sf.to_csv(os.path.join(od, "enrichment_all_bins.csv"), index=False); tp.to_csv(os.path.join(od, "enrichment_top_bins.csv"), index=False)
        bt.to_csv(os.path.join(od, "bin_effects.csv.gz"), index=False, float_format="%.5g"); FG.differential(os.path.join(od, "enrichment_all_bins.png"), sf, primary)
        info["differential"] = dict(bins_used=nb, residual_df=float(df_))
        if df_ < 8: log(f"note: differential mode has {df_:.0f} residual df; bin-level t are unstable")
        log("differential done")

    json.dump(dict(version=__version__, panel=a.panel, panel_definition=P, layers={k: lmeta["layers"][k] for k in P["layers"]}, build=a.build, bin=a.bin,
                   modes=modes, factors=F, reference_group=refcell, effects=terms, covariates=covs, primary=primary, mode_info=info,
                   args=vars(a), n_samples=len(samples), time=time.strftime("%Y-%m-%d %H:%M:%S")),
              open(os.path.join(a.out, "run_manifest.json"), "w"), indent=1, default=str)
    log("done ->", a.out)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="cfdna-tissue-enrich", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "tracks", "count"):
        p = sub.add_parser(name)
        p.add_argument("--panel", default="pan", help="panel name (see 'panels') or a panel JSON path")
        p.add_argument("--build", default="hg38"); p.add_argument("--bin", default="50kb", choices=list(RF.BIN))
        p.add_argument("--reference", default=None, help="reference root (default: ./reference or $CFTE_REFERENCE)")
        if name in ("run", "count"):
            p.add_argument("--samples", required=True, help="CSV with sample_id, frag_path (+ metadata columns)")
            p.add_argument("--cache", default=None, help="folder for per-sample bin counts; reused by later runs")
            p.add_argument("--min-mapq", type=int, default=30); p.add_argument("--mapq-col", type=int, default=None, help="1-based MAPQ column (default: header 'mapq', else column 4 or 5 if numeric; 0 = none)")
            p.add_argument("--min-len", type=int, default=50); p.add_argument("--max-len", type=int, default=500)
        if name == "count":
            p.add_argument("--index", type=int, default=None, help="0-based row of the samples file (for Slurm arrays)")
        if name == "run":
            p.add_argument("--out", required=True)
            p.add_argument("--modes", default="sets,sample,differential")
            p.add_argument("--factor", action="append", default=[], help="COLUMN:REFERENCE:CASE, repeatable")
            p.add_argument("--group", help="shortcut for one factor (with --case, --reference-level)"); p.add_argument("--case"); p.add_argument("--reference-level")
            p.add_argument("--interaction", action="store_true"); p.add_argument("--covariates", default="", help="comma-separated; :num / :cat to force type")
            p.add_argument("--primary", default=None, help="comma-separated primary tracks (overrides the panel's list)")
            p.add_argument("--min-fragments", type=float, default=5e6); p.add_argument("--set-size", type=int, default=None, help="bins per set (default 1000 at 50kb, 5000 at 10kb)")
            p.add_argument("--n-random", type=int, default=200); p.add_argument("--n-perm", type=int, default=2000); p.add_argument("--n-perm-bins", type=int, default=200)
            p.add_argument("--top-n", type=int, default=2000); p.add_argument("--seed", type=int, default=17)
    sub.add_parser("panels")
    a = ap.parse_args(argv)
    {"run": cmd_run, "panels": cmd_panels, "tracks": cmd_tracks, "count": cmd_count}[a.cmd](a)


if __name__ == "__main__":
    main()
