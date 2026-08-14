"""Marker gene identification (scanpy rank_genes_groups)."""

from __future__ import annotations

from typing import Any
import re

import pandas as pd

from scripts.nb04_de.load import ensure_groupby, symbol_lookup_from_adata, resolve_gene_symbol


def _filter_groups(adata, groupby: str, min_cells: int) -> list[str]:
    counts = adata.obs[groupby].astype(str).value_counts()
    keep = counts[counts >= int(min_cells)].index.astype(str).tolist()
    if len(keep) < 2:
        raise ValueError(
            f"Need >=2 groups with >={min_cells} cells in '{groupby}'. Counts:\n{counts.head(20)}"
        )
    return keep


def run_rank_genes(
    adata,
    cfg: dict[str, Any],
    groupby: str | None = None,
    key_added: str | None = None,
    groups: list[str] | None = None,
) -> str:
    """
    Wilcoxon marker detection via scanpy.tl.rank_genes_groups.

    Uses log-normalized `.X` (Module 2/3 preprocessing). Raw counts remain in `.layers['counts']`.
    """
    import scanpy as sc

    m4 = cfg["module2"]
    de = m4["de"]
    groupby = ensure_groupby(adata, groupby or m4["groupby"])
    key = key_added or f"rank_genes_{groupby}"
    min_cells = int(de.get("min_cells_per_group", 50))
    keep = groups or _filter_groups(adata, groupby, min_cells)

    mask = adata.obs[groupby].astype(str).isin(keep)
    ad = adata[mask].copy()
    ad.obs[groupby] = ad.obs[groupby].astype(str).astype("category")

    sc.tl.rank_genes_groups(
        ad,
        groupby=groupby,
        groups=keep,
        reference="rest",
        method=str(de.get("method", "wilcoxon")),
        n_genes=int(de.get("n_genes", 100)),
        pts=True,
        key_added=key,
        use_raw=False,
    )
    adata.uns[key] = ad.uns[key]
    adata.uns[f"{key}_params"] = {
        "groupby": groupby,
        "groups": keep,
        "method": de.get("method", "wilcoxon"),
        "n_genes": int(de.get("n_genes", 100)),
        "min_cells_per_group": min_cells,
        "n_obs_used": int(ad.n_obs),
    }
    return key


def run_rank_genes_for_prerank(
    adata,
    cfg: dict[str, Any],
    groupby: str | None = None,
    *,
    key_added: str | None = None,
    groups: list[str] | None = None,
) -> str:
    """
    Full-universe Wilcoxon ranks for preranked GSEA only.

    Marker tables keep ``de.n_genes``; prerank needs both tails so n_genes = n_vars.
    """
    import scanpy as sc

    m4 = cfg["module2"]
    de = m4["de"]
    groupby = ensure_groupby(adata, groupby or m4["groupby"])
    key = key_added or f"rank_genes_{groupby}_prerank"
    min_cells = int(de.get("min_cells_per_group", 50))
    keep = groups or _filter_groups(adata, groupby, min_cells)

    mask = adata.obs[groupby].astype(str).isin(keep)
    ad = adata[mask].copy()
    ad.obs[groupby] = ad.obs[groupby].astype(str).astype("category")
    n_genes = int(ad.n_vars)

    sc.tl.rank_genes_groups(
        ad,
        groupby=groupby,
        groups=keep,
        reference="rest",
        method=str(de.get("method", "wilcoxon")),
        n_genes=n_genes,
        pts=False,
        key_added=key,
        use_raw=False,
    )
    adata.uns[key] = ad.uns[key]
    adata.uns[f"{key}_params"] = {
        "groupby": groupby,
        "groups": keep,
        "method": de.get("method", "wilcoxon"),
        "n_genes": n_genes,
        "min_cells_per_group": min_cells,
        "n_obs_used": int(ad.n_obs),
        "purpose": "prerank_full_universe",
    }
    return key


def rank_genes_to_frame(adata, key: str) -> pd.DataFrame:
    """Flatten scanpy rank_genes_groups into a tidy table with gene symbols."""
    import scanpy as sc

    df = sc.get.rank_genes_groups_df(adata, group=None, key=key)
    if "names" in df.columns and "gene_id" not in df.columns:
        df = df.rename(columns={"names": "gene_id"})
    lookup = symbol_lookup_from_adata(adata)
    df["gene_symbol"] = df["gene_id"].astype(str).map(lambda g: resolve_gene_symbol(g, lookup))
    # Keep rows even without symbol for DE tables; enrichment filter drops empties

    params = adata.uns.get(f"{key}_params", {})
    df["groupby"] = params.get("groupby", "")
    df["de_method"] = params.get("method", "")
    df["rank_key"] = key
    return df


# Ribosomal proteins with non-RP gene symbols (prefix filter cannot reach these).
_RIBOSOMAL_SYMBOL_EXCEPTIONS = frozenset({"UBA52", "FAU"})


def is_ribosomal_symbol(symbol: str) -> bool:
    """True for cytosolic/mito ribosomal protein gene symbols.

    Prefix rule covers RPS/RPL/MRPS/MRPL. Named exceptions: UBA52 (ubiquitin-L40
    fusion) and FAU (fubi-S30 fusion), which encode ribosomal proteins under
    non-RP names.
    """
    s = str(symbol).upper().strip()
    if not s or s in {"NAN", "NONE"}:
        return False
    if s in _RIBOSOMAL_SYMBOL_EXCEPTIONS:
        return True
    # Match RPS14, RPL34, RPLP0, MRPS6, MRPL12, etc.
    return bool(re.match(r"^(MRP[SL]|RP[SL])(\d+|[A-Z]\d*)", s))


def _want_exclude_ribosomal(cfg: dict[str, Any], for_enrichment: bool = False) -> bool:
    de_flag = bool(cfg["module2"]["de"].get("exclude_ribosomal_genes", True))
    if not for_enrichment:
        return de_flag
    enr_flag = cfg["module2"]["enrichment"].get("exclude_ribosomal_genes", None)
    if enr_flag is None:
        return de_flag
    return bool(enr_flag)


def drop_ribosomal_markers(marker_df: pd.DataFrame, cfg: dict[str, Any] | None = None) -> pd.DataFrame:
    """Remove ribosomal protein genes (and empty symbols) from a reported marker/DE table."""
    if marker_df is None or marker_df.empty:
        out = marker_df.copy() if marker_df is not None else pd.DataFrame()
        out.attrs["n_ribosomal_removed"] = 0
        out.attrs["n_unmapped_symbol_removed"] = 0
        return out

    df = marker_df.copy()
    n_unmapped = 0
    if "gene_symbol" in df.columns:
        missing = df["gene_symbol"].isna() | df["gene_symbol"].astype(str).isin(["", "nan", "None"])
        n_unmapped = int(missing.sum())
        df = df.loc[~missing].copy()

    if cfg is not None and not _want_exclude_ribosomal(cfg, for_enrichment=False):
        df = df.reset_index(drop=True)
        df.attrs["n_ribosomal_removed"] = 0
        df.attrs["n_unmapped_symbol_removed"] = n_unmapped
        df.attrs["exclude_ribosomal_genes"] = False
        return df

    if "gene_symbol" not in df.columns:
        df.attrs["n_ribosomal_removed"] = 0
        df.attrs["n_unmapped_symbol_removed"] = n_unmapped
        return df

    ribo_mask = df["gene_symbol"].map(is_ribosomal_symbol)
    n_removed = int(ribo_mask.sum())
    out = df.loc[~ribo_mask].reset_index(drop=True)
    out.attrs["n_ribosomal_removed"] = n_removed
    out.attrs["n_unmapped_symbol_removed"] = n_unmapped
    out.attrs["exclude_ribosomal_genes"] = True
    return out


def filter_markers_for_enrichment(marker_df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    """
    Build Enrichr input lists from DE markers.

    Applies significance/logFC cuts, requires HUGO symbols, and optionally removes
    ribosomal genes so pathway results are not dominated by generic translation terms.
    """
    de = cfg["module2"]["de"]
    enr = cfg["module2"]["enrichment"]
    max_padj = float(de.get("max_padj", enr.get("cutoff_padj", 0.05)))
    min_lfc = float(de.get("min_logfoldchange", 0.25))
    top_n = int(enr.get("top_genes_per_group", 100))
    drop_ribo = _want_exclude_ribosomal(cfg, for_enrichment=True)

    df = marker_df.copy()
    if "pvals_adj" in df.columns:
        df = df[df["pvals_adj"].astype(float) <= max_padj]
    if "logfoldchanges" in df.columns:
        df = df[df["logfoldchanges"].astype(float) >= min_lfc]
    df = df[df["gene_symbol"].notna()]
    df = df[df["gene_symbol"].astype(str).str.len() > 0]
    df = df[~df["gene_symbol"].astype(str).isin(["nan", "None", ""])]
    # Drop Ensembl-looking symbols for Enrichr (wants HUGO)
    df = df[~df["gene_symbol"].astype(str).str.startswith("ENSG")]

    n_before_ribo = int(df.shape[0])
    if drop_ribo:
        ribo_mask = df["gene_symbol"].map(is_ribosomal_symbol)
        df = df.loc[~ribo_mask].copy()
    n_ribo_removed = n_before_ribo - int(df.shape[0])

    df = (
        df.sort_values(["group", "pvals_adj", "logfoldchanges"], ascending=[True, True, False])
        .groupby("group", observed=True)
        .head(top_n)
        .reset_index(drop=True)
    )
    df.attrs["n_ribosomal_removed"] = n_ribo_removed
    df.attrs["exclude_ribosomal_genes"] = drop_ribo
    return df


def known_marker_presence(adata, cfg: dict[str, Any]) -> pd.DataFrame:
    """Check curated lung markers by HUGO resolution and by Ensembl ID in var_names."""
    from scripts.nb04_de.load import ensembl_base_id

    lookup = {v.upper(): k for k, v in symbol_lookup_from_adata(adata).items()}
    present_symbols = set(symbol_lookup_from_adata(adata).values())
    present_upper = {s.upper() for s in present_symbols}
    var_bases = {ensembl_base_id(str(v)) for v in adata.var_names.astype(str)}
    if "gene_id_no_version" in adata.var.columns:
        var_bases |= {
            ensembl_base_id(str(v)) for v in adata.var["gene_id_no_version"].astype(str)
        }

    # Extra symbol -> Ensembl map when hugo_symbol is empty/NaN in the DE object.
    # Prefer local orthologs reference; fall back to curated IDs for teaching markers.
    external = _symbol_to_ensembl_external(cfg)
    rows = []
    for category, genes in (cfg["module2"].get("known_markers") or {}).items():
        for g in genes:
            ens = lookup.get(g.upper(), "") or external.get(g.upper(), "")
            ens_base = ensembl_base_id(ens) if ens else ""
            ens_in_var = bool(ens_base) and ens_base in var_bases
            rows.append(
                {
                    "category": category,
                    "gene_symbol": g,
                    "symbol_resolved_in_de_object": g.upper() in present_upper,
                    "ensembl_id": ens_base or ens,
                    "ensembl_id_in_var_names": bool(ens_in_var),
                }
            )
    return pd.DataFrame(rows)


# Curated Ensembl IDs for known markers that lack mouse one2one orthologs / empty hugo.
_KNOWN_MARKER_ENSG_FALLBACK: dict[str, str] = {
    "SFTPA1": "ENSG00000122852",
    "SFTPA2": "ENSG00000185303",
    "SFTPC": "ENSG00000168484",
    "AGER": "ENSG00000204305",
    "SCGB1A1": "ENSG00000149021",
    "FOXJ1": "ENSG00000129654",
    "KRT5": "ENSG00000186081",
    "PECAM1": "ENSG00000261371",
    "VWF": "ENSG00000110799",
    "CLDN5": "ENSG00000184113",
    "COL1A1": "ENSG00000108821",
    "COL1A2": "ENSG00000164692",
    "DCN": "ENSG00000011465",
    "LUM": "ENSG00000139329",
    "PTPRC": "ENSG00000081237",
    "CD3D": "ENSG00000167286",
    "MS4A1": "ENSG00000156738",
    "LST1": "ENSG00000204482",
    "NKG7": "ENSG00000105374",
    "ACTA2": "ENSG00000107796",
    "RGS5": "ENSG00000232995",
    "PDGFRB": "ENSG00000113721",
}


def _symbol_to_ensembl_external(cfg: dict[str, Any]) -> dict[str, str]:
    """Map HUGO -> Ensembl using orthologs TSV when present, else curated fallback."""
    from scripts.common.paths import module_root

    out = dict(_KNOWN_MARKER_ENSG_FALLBACK)
    root = module_root(cfg)
    ortho = root / "data" / "reference" / "human_mouse_orthologs_one2one.tsv"
    if ortho.exists():
        try:
            df = pd.read_csv(ortho, sep="\t", comment="#")
            if {"human_symbol", "human_ensembl"}.issubset(df.columns):
                for sym, eid in zip(df["human_symbol"].astype(str), df["human_ensembl"].astype(str)):
                    if sym and not sym.upper().startswith("ENSG"):
                        out[sym.upper()] = eid.split(".")[0]
        except Exception:  # noqa: BLE001 - presence table must not fail the DE path
            pass
    return out


def run_azimuth_focus_de(adata, cfg: dict[str, Any]) -> tuple[str | None, pd.DataFrame]:
    focus = cfg["module2"].get("azimuth_focus") or []
    if not focus or "azimuth_label" not in adata.obs.columns:
        return None, pd.DataFrame()
    available = set(adata.obs["azimuth_label"].astype(str))
    groups = [g for g in focus if g in available]
    min_cells = int(cfg["module2"]["de"].get("min_cells_per_group", 50))
    counts = adata.obs["azimuth_label"].astype(str).value_counts()
    groups = [g for g in groups if int(counts.get(g, 0)) >= min_cells]
    if len(groups) < 2:
        return None, pd.DataFrame()
    key = run_rank_genes(
        adata,
        cfg,
        groupby="azimuth_label",
        key_added="rank_genes_azimuth_focus",
        groups=groups,
    )
    return key, rank_genes_to_frame(adata, key)


def groups_below_min_cells(adata, groupby: str, min_cells: int) -> pd.DataFrame:
    """Report groups dropped by min_cells_per_group (honesty requirement)."""
    counts = adata.obs[groupby].astype(str).value_counts()
    rows = [
        {"group": g, "n_cells": int(n), "min_cells_per_group": int(min_cells), "kept": bool(n >= min_cells)}
        for g, n in counts.items()
    ]
    return pd.DataFrame(rows).sort_values("n_cells", ascending=False)


def wilcoxon_score_ranks_for_prerank(
    adata,
    key: str,
    group: str,
    cfg: dict[str, Any] | None = None,
) -> pd.Series:
    """
    Build a gene-symbol -> Wilcoxon score Series for one group (for gseapy.prerank).

    Prefers HUGO symbols. When cfg requests ribosomal exclusion (same rule as ORA
    inputs), RPS/RPL/MRPS/MRPL and UBA52/FAU are dropped so prerank is not dominated by translation
    terms that the marker tables already filter.
    """
    import scanpy as sc

    df = sc.get.rank_genes_groups_df(adata, group=group, key=key)
    if df.empty:
        return pd.Series(dtype=float)
    if "names" in df.columns and "gene_id" not in df.columns:
        df = df.rename(columns={"names": "gene_id"})
    lookup = symbol_lookup_from_adata(adata)
    df["gene_symbol"] = df["gene_id"].astype(str).map(lambda g: resolve_gene_symbol(g, lookup))
    df = df[df["gene_symbol"].notna()]
    df = df[~df["gene_symbol"].astype(str).isin(["", "nan", "None"])]
    df = df[~df["gene_symbol"].astype(str).str.startswith("ENSG")]
    if cfg is not None and _want_exclude_ribosomal(cfg, for_enrichment=True):
        df = df[~df["gene_symbol"].map(is_ribosomal_symbol)]
    score_col = "scores" if "scores" in df.columns else None
    if score_col is None:
        return pd.Series(dtype=float)
    s = (
        df.groupby("gene_symbol", observed=True)[score_col]
        .mean()
        .sort_values(ascending=False)
    )
    s.name = f"wilcoxon_score_{group}"
    return s.astype(float)
