"""Module 4 plotting helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from scripts.common.plotting import save_figure
from scripts.nb04_de.load import gene_symbol_series, symbol_lookup_from_adata, resolve_gene_symbol


def plot_marker_dotplot(
    adata,
    marker_df: pd.DataFrame,
    groupby: str,
    path: Path | str,
    n_genes: int = 5,
    dpi: int = 150,
) -> Path | None:
    import scanpy as sc

    if marker_df.empty or groupby not in adata.obs.columns:
        return None
    # pick top n genes per group by padj
    df = marker_df.copy()
    if "pvals_adj" in df.columns:
        df = df.sort_values(["group", "pvals_adj"], ascending=[True, True])
    top = df.groupby("group", observed=True).head(n_genes)
    lookup = symbol_lookup_from_adata(adata)
    # Prefer genes with resolvable HUGO symbols for readable plots
    var_ids = []
    for gid in top["gene_id"].astype(str):
        if gid in adata.var_names and resolve_gene_symbol(gid, lookup):
            var_ids.append(gid)
    seen = set()
    genes = []
    for g in var_ids:
        if g not in seen:
            seen.add(g)
            genes.append(g)
    if len(genes) < 2:
        return None

    ad = adata[:, genes].copy()
    labels = [resolve_gene_symbol(g, lookup) or g for g in genes]
    ad.var_names = pd.Index(labels)
    if ad.var_names.duplicated().any():
        ad.var_names_make_unique()

    sc.pl.dotplot(ad, var_names=list(ad.var_names), groupby=groupby, show=False, dendrogram=False)
    return save_figure(path, dpi=dpi, close=True)


def plot_enrichment_bar(
    enrichment_df: pd.DataFrame,
    path: Path | str,
    n_terms: int = 8,
    dpi: int = 150,
) -> Path | None:
    import matplotlib.pyplot as plt

    if enrichment_df is None or enrichment_df.empty:
        return None
    df = enrichment_df.copy()
    if "Term" in df.columns:
        df = df[~df["Term"].astype(str).str.startswith("ENRICHMENT_FAILED")]
    padj_col = "Adjusted P-value" if "Adjusted P-value" in df.columns else None
    if padj_col is None or "group" not in df.columns:
        return None
    df = df[df[padj_col].astype(float) <= 0.05]
    if df.empty:
        return None

    groups = list(df["group"].astype(str).unique())
    n = len(groups)
    fig, axes = plt.subplots(n, 1, figsize=(10, max(2.8 * n, 3.5)), squeeze=False)
    for ax, group in zip(axes[:, 0], groups):
        sub = df[df["group"].astype(str) == group].sort_values(padj_col).head(n_terms)
        terms = [
            (t[:80] + "...") if len(t) > 80 else t
            for t in sub["Term"].astype(str).tolist()[::-1]
        ]
        vals = (-np.log10(sub[padj_col].astype(float).clip(lower=1e-300))).tolist()[::-1]
        ax.barh(terms, vals, color="#4C78A8")
        ax.set_xlabel("-log10(adjusted p)")
        ax.set_title(f"Enrichment -- {group}")
    fig.subplots_adjust(left=0.42, right=0.98, hspace=0.45)
    return save_figure(path, dpi=dpi, close=True)


def plot_known_marker_dotplot(
    adata,
    cfg: dict[str, Any],
    groupby: str,
    path: Path | str,
    dpi: int = 150,
) -> Path | None:
    """Dotplot of curated lung markers (schema Section 3 validation).

    Resolve markers by HUGO when present; otherwise by Ensembl ID in var_names
    (object symbols + curated/ortholog fallback). Never fills symbols from GTEx.
    """
    import matplotlib.pyplot as plt
    import scanpy as sc

    from scripts.nb04_de.markers import _symbol_to_ensembl_external

    lookup = symbol_lookup_from_adata(adata)
    inv = {v.upper(): k for k, v in lookup.items()}
    external = _symbol_to_ensembl_external(cfg)
    # map ensembl base -> actual var_name in object
    base_to_var = {ensembl_base_from_var(v): v for v in adata.var_names.astype(str)}
    genes_ids = []
    labels = []
    omitted = []
    for _cat, genes in (cfg["module2"].get("known_markers") or {}).items():
        for g in genes:
            base = inv.get(g.upper()) or ensembl_base_from_var(external.get(g.upper(), ""))
            if not base:
                omitted.append(g)
                continue
            vid = base_to_var.get(base)
            if vid and vid in adata.var_names:
                genes_ids.append(vid)
                labels.append(g)
            else:
                omitted.append(g)
    seen = set()
    ids_u, labs_u = [], []
    for gid, lab in zip(genes_ids, labels):
        if gid not in seen:
            seen.add(gid)
            ids_u.append(gid)
            labs_u.append(lab)
    if len(ids_u) < 2 or groupby not in adata.obs.columns:
        return None
    ad = adata[:, ids_u].copy()
    ad.var_names = pd.Index(labs_u)
    if ad.var_names.duplicated().any():
        ad.var_names_make_unique()
    sc.pl.dotplot(ad, var_names=list(ad.var_names), groupby=groupby, show=False)
    fig = plt.gcf()
    if omitted:
        fig.text(
            0.5,
            0.01,
            f"{len(omitted)} markers omitted (no HUGO/Ensembl match in object): "
            + ", ".join(omitted),
            ha="center",
            fontsize=8,
            color="#444444",
        )
        fig.subplots_adjust(bottom=0.12)
    return save_figure(path, dpi=dpi, close=True)


def ensembl_base_from_var(vid: str) -> str:
    from scripts.nb04_de.load import ensembl_base_id

    return ensembl_base_id(vid)
