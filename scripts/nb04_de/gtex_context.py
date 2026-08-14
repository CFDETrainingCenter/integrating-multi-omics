"""Optional GTEx lung TPM context for Module 4 markers."""

from __future__ import annotations

from typing import Any

import pandas as pd

from scripts.nb03_integration.gtex import gtex_subject_means, load_gtex_lung_tpm, strip_ensembl_version


def marker_gtex_context(marker_df: pd.DataFrame, cfg: dict[str, Any], top_n: int = 20) -> pd.DataFrame:
    """Attach GTEx adult lung mean TPM to top markers (tissue-level context only)."""
    if not cfg["module2"].get("gtex_context", True) or marker_df.empty:
        return pd.DataFrame()

    gtex = load_gtex_lung_tpm(cfg)
    means = gtex_subject_means(gtex, cfg)

    # top markers per group by padj then lfc
    df = marker_df.copy()
    if "pvals_adj" in df.columns:
        df = df.sort_values(["group", "pvals_adj", "logfoldchanges"], ascending=[True, True, False])
    top = df.groupby("group", observed=True).head(top_n).copy()
    top["gene_id_no_version"] = strip_ensembl_version(top["gene_id"]).astype(str)
    top["gtex_lung_mean_tpm"] = top["gene_id_no_version"].map(means)
    top["gtex_note"] = (
        "Tissue-level GTEx adult lung mean TPM -- supports presence in lung bulk, "
        "not cell-type specificity"
    )
    return top
