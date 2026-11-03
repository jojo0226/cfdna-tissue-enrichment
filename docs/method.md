# Method

Every panel runs exactly these steps; only the atlas layers and the pre-registered tracks differ.

## 1. Coverage

1. **Grid.** hg38 autosomes in 50 kb (45,595 bins) or 10 kb (241,026 bins). Bins are removed if they overlap the
   ENCODE blacklist (ENCFF356LFX) by any amount, overlap an assembly gap, are > 10% N, or have 45-mer mappability < 0.80.
2. **Counting.** Fragment midpoints per bin; MAPQ ≥ 30 and length 50–500 bp. Midpoints in removed bins are dropped.
3. **Normalisation** (identical to `cfdna-enrich prep.py`): log2(CPM + 1) → LOWESS on bin GC (frac 0.10,
   3 iterations, level kept) → subtraction of a centred ~2 Mb rolling median per chromosome (removes copy number and
   other broad trends).

## 2. Track specificity

For each **layer** separately (one source, one assay), every bin's openness is standardised across that layer's
tracks, and each track's specificity is its value minus the mean of the layer's other tracks. Specificity says how
much more open a bin is in that track than in the other tracks *of the same layer*. Layers are never mixed in this
step, so DNase and ATAC values are not compared directly.

Panels can add derived tracks:

- **groups**: mean specificity of member tracks (e.g. `Wgroup_trophoblast_SCT_VCT`, `Tgroup_glioma`);
- **axes**: mean of one set of tracks minus another (e.g. `origin_fetal_specific`); a cross-layer axis
  (e.g. `glioma_vs_blood`, TCGA tracks minus ENCODE blood tracks) is allowed but marked exploratory.

!!! note "Within-layer means relative to the layer"
    `T_GBM` specificity is "open in GBM more than in the other TCGA cancer types", not "more than in blood".
    Use a cross-layer axis or the `C_` blood tracks alongside it when the question is tumour vs blood.

## 3. The three modes

| mode | per-sample score | method |
|---|---|---|
| `sets` | yes | Each track's top-N specific bins (1,000 at 50 kb, 5,000 at 10 kb). Score = mean normalised coverage over the set minus the mean of 200 random bin sets matched on GC decile × mappability quintile. This is the 10-kb marker method of `cfdna-enrich markers`, with sets defined from the atlas. |
| `sample` | yes | Each bin is turned into a z-score against the **reference group** (median, 1.4826 × MAD; reference samples left out of their own reference when ≥ 4). E = Spearman correlation of a sample's z over all bins with the track's specificity. Uses all bins, nothing selected. |
| `differential` | no | Each bin is fitted with the group model (moderated t). Per track: (a) **all bins**: Spearman correlation of bin t with specificity; (b) **top bins**: mean specificity of the 2,000 bins most up and most down. Both are tested against label permutations (bins refitted and top bins re-selected each time). A z against GC × mappability matched random bins is reported for description only. |

Positive values mean the same in every mode: more coverage, in cases (or in that sample), at bins specific to the track.

## 4. Group design and tests

- `--factor COLUMN:REFERENCE:CASE` defines a two-level factor; repeat for factorial designs. The **reference group**
  is every factor at its reference level.
- **Effect tests** (sets and sample modes): OLS of each score on the factors + covariates; p from OLS and from
  permutations of one factor's labels *within* levels of the other factors (covariates stay with their samples).
- **Pairwise tests**: every other group vs the reference group (difference, AUC, Mann–Whitney p, permutation p).
- **Multiple testing**: BH across all tracks (`q_perm`) and separately across the panel's primary tracks (`q_perm_primary`).

## 5. Why permutation p-values for top bins

Top differential bins are selected from noisy data, so they are not exchangeable with random bins even when GC and
mappability are matched. On null data, matched-random p-values flagged 6–14% of tracks at p < 0.05, against ~6% for
the permutation p-values. The permutation null repeats the selection on shuffled labels.
