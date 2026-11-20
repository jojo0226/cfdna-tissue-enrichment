"""markers mode: cfDNA coverage at cell-type marker bins (10 kb; from cfdna_chromatin.enrichment markers.py).

marker table (csv.gz): cell_type,key,rank[,panel,...]; one row per (set, bin), rank 1 = most specific.
score[sample, set] = mean preprocessed coverage over the set's top-N bins.
Random null: sets of the same size drawn to match the markers' GC deciles; any genome-wide coverage shift moves them too.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

LINEAGE = {"B": ["Bulk_B", "Naive_B", "Mem_B", "Plasmablasts"],
           "CD4 T": ["Th1_precursors", "Th2_precursors", "Th17_precursors", "Regulatory_T", "Naive_Tregs", "Memory_Tregs",
                     "Naive_Teffs", "Memory_Teffs", "Follicular_T_Helper", "Effector_CD4pos_T"],
           "CD8 / γδ T": ["CD8pos_T", "Naive_CD8_T", "Effector_memory_CD8pos_T", "Central_memory_CD8pos_T", "Gamma_delta_T"],
           "NK": ["Mature_NK", "Immature_NK", "Memory_NK"], "myeloid": ["Monocytes", "Myeloid_DCs", "pDCs", "Neutrophil"],
           "CNS cells": ["Neurons", "Astrocyte", "OPC", "Oligodendrocyte"], "microglia": ["Microglia"]}
GLIOMA_TCGA = ("TCGA_GBM", "TCGA_LGG")
NICE = {"Bulk_B": "B (bulk)", "Naive_B": "B naive", "Mem_B": "B memory", "Plasmablasts": "Plasmablast", "Th1_precursors": "Th1 precursor",
        "Th2_precursors": "Th2 precursor", "Th17_precursors": "Th17 precursor", "Regulatory_T": "Treg", "Naive_Tregs": "Treg naive",
        "Memory_Tregs": "Treg memory", "Naive_Teffs": "CD4 Teff naive", "Memory_Teffs": "CD4 Teff memory", "Follicular_T_Helper": "Tfh",
        "Effector_CD4pos_T": "CD4 T effector", "CD8pos_T": "CD8 T", "Naive_CD8_T": "CD8 T naive", "Effector_memory_CD8pos_T": "CD8 T eff. memory",
        "Central_memory_CD8pos_T": "CD8 T central memory", "Gamma_delta_T": "γδ T", "Mature_NK": "NK mature", "Immature_NK": "NK immature",
        "Memory_NK": "NK memory", "Monocytes": "Monocyte", "Myeloid_DCs": "Myeloid DC", "pDCs": "pDC", "Neutrophil": "Neutrophil",
        "Neurons": "Neuron", "Astrocyte": "Astrocyte", "OPC": "OPC", "Oligodendrocyte": "Oligodendrocyte", "Microglia": "Microglia"}


def lineage(s: str) -> str:
    for l, cs in LINEAGE.items():
        if s in cs:
            return l
    if s in GLIOMA_TCGA:
        return "tumour (TCGA)"
    if s.startswith("TCGA_"):
        return "other tumour (TCGA)"
    return "other"


def nice(s: str) -> str:
    return NICE.get(s, s)


def load(path: str, keys: pd.Index, include: str | None = None, exclude: str | None = None) -> pd.DataFrame:
    mk = pd.read_csv(path)
    if include:
        mk = mk[mk.cell_type.str.contains(include, regex=True)]
    if exclude:
        mk = mk[~mk.cell_type.str.contains(exclude, regex=True)]
    pos = pd.Series(np.arange(len(keys)), index=keys)
    mk = mk[mk.key.isin(pos.index)].copy(); mk["i"] = mk.key.map(pos).astype(int)
    return mk


def scores(X: np.ndarray, mk: pd.DataFrame, top: int) -> pd.DataFrame:
    sets = list(dict.fromkeys(mk.cell_type))
    return pd.DataFrame({c: np.nanmean(X[:, mk[(mk.cell_type == c) & (mk["rank"] <= top)].i.to_numpy()], 1) for c in sets})


def random_sets(gc: np.ndarray, mk: pd.DataFrame, top: int, n: int, seed: int = 17) -> dict[str, list[np.ndarray]]:
    dec = pd.qcut(gc, 10, labels=False, duplicates="drop"); rng = np.random.default_rng(seed)
    pool = {d: np.where(dec == d)[0] for d in np.unique(dec)}
    out = {}
    for c in dict.fromkeys(mk.cell_type):
        idx = mk[(mk.cell_type == c) & (mk["rank"] <= top)].i.to_numpy(); cnt = pd.Series(dec[idx]).value_counts()
        out[c] = [np.concatenate([rng.choice(pool[d], k, replace=False) for d, k in cnt.items()]) for _ in range(n)]
    return out


def jaccard(mk: pd.DataFrame, top: int) -> pd.DataFrame:
    sets = {c: set(mk[(mk.cell_type == c) & (mk["rank"] <= top)].key) for c in dict.fromkeys(mk.cell_type)}
    ks = list(sets)
    return pd.DataFrame([[len(sets[a] & sets[b]) / max(1, len(sets[a] | sets[b])) for b in ks] for a in ks], index=ks, columns=ks)
