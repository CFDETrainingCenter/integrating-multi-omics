"""Summed HuBMAP counts-layer per donor for H_AGING_HUBMAP (Module 4 rebuild).

Deliberately separate from the composition-weighted / max-across level bridge.
The integrated ``layers['counts']`` values are not integer UMIs (fractional
pipeline outputs); sums are rounded to non-negative integers before pydeseq2,
same honesty rule as RSEM expected counts.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from scripts.common.paths import module_root


def _resolve_h5ad(cfg: dict[str, Any]) -> Path:
    data = (cfg.get("module4") or {}).get("data") or {}
    rel = data.get("hubmap_integrated_h5ad") or (
        "data/processed/integrated/module2_hubmap_lung_integrated.h5ad"
    )
    p = Path(rel)
    return p if p.is_absolute() else module_root(cfg) / p


def _resolve_out(cfg: dict[str, Any]) -> Path:
    data = (cfg.get("module4") or {}).get("data") or {}
    rel = data.get("hubmap_donor_counts") or (
        "outputs/tables/module4_hubmap_donor_summed_counts.tsv"
    )
    p = Path(rel)
    return p if p.is_absolute() else module_root(cfg) / p


def summed_counts_by_donor(
    cfg: dict[str, Any],
    *,
    force: bool = False,
) -> pd.DataFrame:
    """
    Return genes x donors integer count matrix (rounded sums of layers['counts']).

    Writes/reads ``module4_hubmap_donor_summed_counts.tsv`` unless ``force``.
    """
    out_path = _resolve_out(cfg)
    if out_path.exists() and not force:
        df = pd.read_csv(out_path, sep="\t", index_col=0)
        return df

    import anndata as ad
    import scipy.sparse as sp

    h5ad = _resolve_h5ad(cfg)
    if not h5ad.exists():
        raise FileNotFoundError(f"Missing HuBMAP integrated object: {h5ad}")

    adata = ad.read_h5ad(h5ad)
    if "counts" not in adata.layers:
        raise KeyError(f"No layers['counts'] in {h5ad}")
    if "donor_label" not in adata.obs.columns:
        raise KeyError(f"No obs['donor_label'] in {h5ad}")

    X = adata.layers["counts"]
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    else:
        X = X.tocsr()

    donors = adata.obs["donor_label"].astype(str)
    cols: dict[str, np.ndarray] = {}
    for donor in sorted(donors.unique()):
        idx = np.where(donors.values == donor)[0]
        summed = np.asarray(X[idx].sum(axis=0)).ravel().astype(np.float64)
        # Round like RSEM expected counts; layer is not integer UMIs.
        cols[donor] = np.maximum(np.rint(summed), 0).astype(int)

    genes = adata.var_names.astype(str)
    sym_col = "hugo_symbol" if "hugo_symbol" in adata.var.columns else (
        "gene_symbol" if "gene_symbol" in adata.var.columns else None
    )
    symbols = (
        adata.var[sym_col].astype(str).str.upper() if sym_col is not None else None
    )

    df = pd.DataFrame(cols, index=genes)
    df.index.name = "gene_id_no_version"
    if symbols is not None:
        df.insert(0, "gene_symbol", symbols.values)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t")
    return df


def donor_age_meta(cfg: dict[str, Any]) -> pd.DataFrame:
    """One row per donor_label with age_value from the integrated object (or teaching blocks)."""
    aging = (cfg.get("module4") or {}).get("hubmap_aging") or {}
    older = [str(x) for x in aging.get("older_donors") or ["Donor_3", "Donor_4"]]
    younger = [str(x) for x in aging.get("younger_donors") or ["Donor_1", "Donor_2"]]

    # Prefer ages from teaching blocks table (stable, small)
    blocks = module_root(cfg) / "outputs/tables/module1_snrna_teaching_blocks.tsv"
    if blocks.exists():
        tb = pd.read_csv(blocks, sep="\t")
        ages = (
            tb[["donor_label", "age_value"]]
            .drop_duplicates("donor_label")
            .set_index("donor_label")["age_value"]
            .astype(float)
        )
    else:
        ages = pd.Series(
            {"Donor_1": 37.0, "Donor_2": 25.0, "Donor_3": 56.8, "Donor_4": 52.96}
        )

    rows = []
    for d in younger:
        rows.append(
            {
                "donor_label": d,
                "condition": "younger",
                "age_value": float(ages.get(d, float("nan"))),
                "contrast_id": "H_AGING_HUBMAP",
                "arm_role": "younger",
            }
        )
    for d in older:
        rows.append(
            {
                "donor_label": d,
                "condition": "older",
                "age_value": float(ages.get(d, float("nan"))),
                "contrast_id": "H_AGING_HUBMAP",
                "arm_role": "older",
            }
        )
    return pd.DataFrame(rows)
