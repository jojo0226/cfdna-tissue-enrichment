# Panels and layers

## Layers

hg38: C, T, W. hg19: C, T (the atlas's original build; W is not built for hg19). Choose with .

### hg38

| ID | Layer | Assay | Tracks | Source | Notes |
|---|---|---|---|---|---|
| `C` | ENCODE tissues and cell types | DNase-seq | 73 | ENCODE DNase narrowPeak (72 biosamples: tissues, blood cells, endothelium, trophoblast) + one neutrophil track; `cfdna-chromatin-xai` openness atlas | 10 kb lifted from hg19 by bin midpoint; 50 kb hg38 provenance not documented |
| `T` | TCGA tumours | ATAC-seq | 23 | TCGA-ATAC cancer-type peak calls (Corces et al., *Science* 2018) | split from the same atlas file as C |
| `W` | Maternal–fetal interface | snATAC-seq pseudobulk | 32 | Wang et al., *Nature* 653:167–179 (2026); 5 kb window counts prepared by Paz Polak | placed on hg38 natively; erythroid dropped (67 nuclei) |

Track names in outputs are `<layer>_<track>`, e.g. `C_placenta`, `T_GBM`, `W_SCT`.
Cell-type metadata for W (origin, purity, nuclei) is in `reference/hg38/tracks_W.csv`. Full provenance is in `reference/hg38/layers.json`.

## Panels

Panels live in `panels/*.json`. A panel lists layers, derived tracks and pre-registered primary tracks.

### `prenatal` (C + W)
- Groups: `Wgroup_trophoblast_SCT_VCT`, `Wgroup_EVT`, `Wgroup_decidual_stroma`, `Wgroup_maternal_immune`,
  `Wgroup_maternal_endothelial`, `Wgroup_fetal_endothelial_stroma`.
- Axis: `origin_fetal_specific` (11 clean fetal minus 7 clean maternal cell types; purity ≥ 0.90, ≥ 500 nuclei) and its mirror `origin_maternal_specific`.
- Primary: the four main W groups, the two origin tracks, `C_placenta`, `C_trophoblast cell`, `C_neutrophil`.

### `tumor` (C + T)
- Groups: `Tgroup_glioma`, `Tgroup_kidney`, `Tgroup_lung`, `Tgroup_gastrointestinal`, `Tgroup_squamous`.
- No primary tracks: set `--primary` for the tumour type under study.

### `glioma` (C + T)
- As `tumor`, plus the exploratory cross-layer axis `glioma_vs_blood` (TCGA GBM/LGG minus ENCODE blood cells).
- Primary: `Tgroup_glioma`, `T_GBM`, `T_LGG`, `C_brain`, `C_frontal cortex`, `C_spinal cord`, `C_neutrophil`, `C_CD14-positive monocyte`, `glioma_vs_blood`.

### `pan` (C + T + W)
All groups from the other panels, no primary tracks. For exploration.

## Make your own panel

Copy a panel JSON and edit it. Fields:

```json
{
 "name": "my_panel",
 "description": "one line",
 "layers": ["C", "T"],
 "groups": {"Tgroup_glioma": {"layer": "T", "members": ["GBM", "LGG"]}},
 "axes": {"glioma_vs_blood": {"plus_layer": "T", "plus": ["GBM", "LGG"],
                              "minus_layer": "C", "minus": ["neutrophil", "B cell"]}},
 "primary": ["Tgroup_glioma", "T_GBM", "C_neutrophil"]
}
```

Run with `--panel path/to/my_panel.json`. Fix the primary list **before** looking at outcome data.
