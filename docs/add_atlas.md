# Add a new atlas

A new atlas is a new **layer**. The method does not change.

1. **Make the layer table** on the existing grid: `reference/hg38/layer_<ID>_50kb.csv.gz` (and `_10kb`), with a
   `key` column matching `grid_<bin>.csv.gz` (e.g. `chr1:900000-950000`) and one numeric column per track.
   Any openness scale works (peak coverage, signal, CPM); specificity re-standardises each bin across the layer's tracks.
   Rules:
    - one source and one assay per layer (do not mix DNase and ATAC in one layer);
    - build from hg38 data, or place windows on hg38 natively rather than lifting bin by bin where possible;
    - bins with no data stay NaN, never 0.
2. **Describe it** in `reference/hg38/layers.json` (name, assay, n_tracks, source, values, provenance, citation;
   `track_meta` for an optional per-track metadata CSV).
3. **Use it** from a panel (`"layers": ["C", "<ID>"]`) and add groups, axes and primary tracks as needed.
4. **Check it** with a known contrast before using it on outcomes, and record the result in the docs.

`build/build_layers.py` shows how the C, T and W layers were made and is the template for new ones.

Planned layers: Calderon 2019 immune cell ATAC (resting and stimulated), Corces 2016 haematopoiesis ATAC
(monocyte, erythroid, progenitors), peripheral-neutrophil ATAC (GSE153520), CATlas adult snATAC (endothelium,
hepatocyte), GSE174367 brain snATAC.
