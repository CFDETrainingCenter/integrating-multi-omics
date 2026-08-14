"""Load Module 2 integrated AnnData for DE / pathways."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from scripts.common.io import read_h5ad
from scripts.common.paths import module_root


def module2_de_input_path(cfg: dict[str, Any]) -> Path:
    rel = cfg["module2"]["input_h5ad"]
    return module_root(cfg) / rel


def load_integrated(cfg: dict[str, Any]):
    path = module2_de_input_path(cfg)
    if not path.exists():
        raise FileNotFoundError(
            f"Module 2 integrated input missing: {path}. "
            "Run Module 2 integration first "
            "(scripts/run_modules_smoke.py --module 3 for the integration half)."
        )
    adata = read_h5ad(path)
    # Prefer Harmony UMAP for feature plots if present
    if "X_umap_harmony" in adata.obsm and "X_umap" in adata.obsm:
        adata.obsm["X_umap"] = adata.obsm["X_umap_harmony"].copy()
    feature_space_before = describe_feature_space(adata, label="integrated_pre_intron_filter")
    if cfg.get("module2", {}).get("exclude_intron_features", True):
        adata = exclude_intron_features(adata)
    feature_space_after = describe_feature_space(adata, label="de_object")
    adata.uns["module2_feature_space"] = {
        "integrated_object": feature_space_before,
        "de_object": feature_space_after,
        "exclude_intron_features": bool(cfg.get("module2", {}).get("exclude_intron_features", True)),
    }
    # Module 2 must not read GTEx (Phase 2). Optional GTEx Ensembl->HUGO fill is
    # only for Module 4 when gtex_context is enabled.
    if cfg.get("module2", {}).get("gtex_context", False):
        adata = annotate_hugo_from_gtex(adata, cfg)
    else:
        adata.uns["hugo_symbol_gtex_fill"] = {
            "ok": False,
            "skipped": True,
            "reason": "gtex_context=false (Module 2 uses object hugo_symbol only)",
        }
    return adata, path

def ensure_groupby(adata, groupby: str) -> str:
    if groupby not in adata.obs.columns:
        raise KeyError(
            f"groupby '{groupby}' not in adata.obs. "
            f"Available: {[c for c in adata.obs.columns if c in {'cell_class','azimuth_label','leiden','donor_label'}]}"
        )
    adata.obs[groupby] = adata.obs[groupby].astype("category")
    return groupby


def gene_symbol_series(adata):
    """Return per-gene symbols (hugo if present), aligned to var_names."""
    if "hugo_symbol" in adata.var.columns:
        sym = adata.var["hugo_symbol"].astype(str).replace({"nan": "", "None": ""})
    else:
        sym = pd.Series([""] * adata.n_vars, index=adata.var_names)
    empty = sym.isna() | (sym.str.len() == 0) | (sym == "nan")
    return sym.where(~empty, other="")


def ensembl_base_id(gene_id: str) -> str:
    """Strip version suffix and optional Salmon intron `-I` tag."""
    gid = str(gene_id)
    if gid.endswith("-I"):
        gid = gid[:-2]
    if "." in gid:
        gid = gid.split(".", 1)[0]
    return gid


def annotate_hugo_from_gtex(adata, cfg: dict[str, Any]):
    """
    Fill empty `adata.var['hugo_symbol']` using GTEx Ensembl->HUGO map.

    Does not overwrite non-empty existing symbols. Safe if GTEx file is missing
    (returns adata unchanged aside from metadata note).
    """
    from scripts.nb03_integration.gtex import load_gtex_lung_tpm, strip_ensembl_version

    try:
        gtex = load_gtex_lung_tpm(cfg)
    except Exception as exc:  # noqa: BLE001
        adata.uns["hugo_symbol_gtex_fill"] = {"ok": False, "error": str(exc)}
        return adata

    base_ids = strip_ensembl_version(pd.Index(adata.var_names.astype(str))).to_numpy()
    gtex_map = gtex["gene_symbol"].astype(str).to_dict()
    from_gtex = pd.Series(
        ["" if str(gtex_map.get(b, "")).lower() in {"", "nan", "none"} else str(gtex_map.get(b, "")) for b in base_ids],
        index=adata.var_names,
        dtype=object,
    )

    existing = gene_symbol_series(adata)
    n_empty_before = int((existing == "").sum())
    filled = existing.where(existing != "", other=from_gtex)
    n_empty_after = int((filled == "").sum())
    adata.var["hugo_symbol"] = filled.to_numpy()
    adata.uns["hugo_symbol_gtex_fill"] = {
        "ok": True,
        "n_empty_before": n_empty_before,
        "n_empty_after": n_empty_after,
        "n_filled": n_empty_before - n_empty_after,
    }
    return adata


def symbol_lookup_from_adata(adata) -> dict[str, str]:
    """Map versionless Ensembl ID -> HUGO using non-empty hugo_symbol rows."""
    sym = gene_symbol_series(adata)
    lookup: dict[str, str] = {}
    for vid, s in zip(adata.var_names.astype(str), sym.astype(str)):
        if not s or s.startswith("ENSG"):
            continue
        base = ensembl_base_id(vid)
        if base not in lookup or not str(vid).endswith("-I"):
            lookup[base] = s
    return lookup


def resolve_gene_symbol(gene_id: str, lookup: dict[str, str]) -> str:
    return lookup.get(ensembl_base_id(gene_id), "")


def exclude_intron_features(adata):
    """Drop Salmon intron quant features (`*-I`) for cleaner marker teaching."""
    from scripts.common.runtime import feature_space_counts

    before = feature_space_counts(adata.var_names)
    keep = ~pd.Index(adata.var_names.astype(str)).str.endswith("-I")
    out = adata[:, keep].copy()
    after = feature_space_counts(out.var_names)
    out.uns["feature_space"] = {
        "before_intron_filter": before,
        "after_intron_filter": after,
        "n_intron_removed": int(before["n_features_intron"]),
        "exclude_intron_features": True,
    }
    return out


def describe_feature_space(adata, *, label: str = "object") -> dict[str, Any]:
    """Report exonic vs intron feature counts for provenance / Module 4 bridge honesty."""
    from scripts.common.runtime import feature_space_counts

    counts = feature_space_counts(adata.var_names)
    stored = adata.uns.get("feature_space") if hasattr(adata, "uns") else None
    return {
        "label": label,
        **counts,
        "composition_weighted_features": (
            f"{counts['n_features_total']} "
            f"({counts['n_features_exonic']} exonic, {counts['n_features_intron']} intron)"
        ),
        "intron_filter": stored,
    }
