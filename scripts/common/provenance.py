"""Provenance helpers shared across modules."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from scripts.common.io import donor_and_dataset_ids, metadata_kv, read_hubmap_metadata_tsv


ANALYSIS_FILE_ROLES = {
    "raw_expr.h5ad": "raw/reference counts for QC",
    "expr.h5ad": "processed expression matrix (often velocity-oriented)",
    "secondary_analysis.h5ad": "normalized object with Azimuth annotations",
    "scvelo_annotated.h5ad": "optional RNA velocity extension",
    "cell_by_bin.h5ad": "snATAC cell-by-bin matrix",
    "cell_by_gene.h5ad": "snATAC gene activity (often MAGIC-smoothed)",
}


def file_role(name: str) -> str:
    return ANALYSIS_FILE_ROLES.get(name, "other/supporting file")


def parse_block_metadata(tsv_path: Path) -> dict[str, Any]:
    df = read_hubmap_metadata_tsv(tsv_path)
    donor_ids, dataset_ids = donor_and_dataset_ids(df)
    donor_kv = metadata_kv(df, entity="Donor")
    assay_kv = metadata_kv(df, entity=None)
    # Prefer assay entity rows when available
    assay_entities = [e for e in df["Entity"].dropna().unique() if e != "Donor"]
    if assay_entities:
        assay_kv = metadata_kv(df, entity=assay_entities[0])
    return {
        "tsv_path": str(tsv_path),
        "donor_ids": donor_ids,
        "dataset_ids": dataset_ids,
        "primary_dataset_id": dataset_ids[0] if dataset_ids else "",
        "donor_id": donor_ids[0] if donor_ids else "",
        "age_value": donor_kv.get("age_value", ""),
        "sex": donor_kv.get("sex", ""),
        "race": donor_kv.get("race", ""),
        "ethnicity": donor_kv.get("ethnicity", ""),
        "assay_description": assay_kv.get("description", ""),
        "dataset_type": assay_kv.get("dataset_type", ""),
        "rnaseq_assay_method": assay_kv.get("rnaseq_assay_method", ""),
    }
