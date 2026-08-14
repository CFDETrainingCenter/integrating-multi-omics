"""Concatenate Module 2 QC'd AnnData objects across donors."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

from scripts.common.io import read_h5ad
from scripts.common.paths import resolve, portable_path


def module2_h5ad_for(cfg: dict[str, Any], block_id: str) -> Path:
    """Resolve Module 1 QC h5ad (module1_ prefix; legacy module2_ fallback)."""
    from scripts.nb02_qc.export import resolve_qc_h5ad

    return resolve_qc_h5ad(cfg, block_id)


def load_module2_inputs(cfg: dict[str, Any]) -> list[tuple[dict[str, str], Any]]:
    loaded = []
    for item in cfg["module2"]["inputs"]:
        path = module2_h5ad_for(cfg, item["block_id"])
        if not path.exists():
            raise FileNotFoundError(
                f"Missing Module 2 output for {item}: {path}. "
                "Run Module 2 for this donor/block first."
            )
        adata = read_h5ad(path)
        loaded.append((item, adata))
    return loaded


def concatenate_donors(cfg: dict[str, Any]):
    """Concatenate Module 2 objects on shared genes; preserve donor metadata."""
    loaded = load_module2_inputs(cfg)
    objects = []
    summary_rows = []

    for item, adata in loaded:
        # ensure required obs
        if "donor_id" not in adata.obs.columns:
            label = str(item["donor_label"])
            if label.lower().startswith("donor_"):
                dkey = f"donor_{label.split('_', 1)[1]}"
            else:
                dkey = "donor_1"
            adata.obs["donor_id"] = cfg.get(dkey, {}).get("donor_id", item["donor_label"])
        if "donor_label" not in adata.obs.columns:
            adata.obs["donor_label"] = item["donor_label"]
        if "block_id" not in adata.obs.columns:
            adata.obs["block_id"] = item["block_id"]
        if "sample_id" not in adata.obs.columns:
            adata.obs["sample_id"] = item["block_id"]

        # unique obs names across donors
        adata.obs_names = [f"{item['block_id']}_{x}" for x in adata.obs_names.astype(str)]
        objects.append(adata)
        summary_rows.append(
            {
                "donor_label": item["donor_label"],
                "donor_id": adata.obs["donor_id"].iloc[0],
                "block_id": item["block_id"],
                "n_nuclei": adata.n_obs,
                "n_genes": adata.n_vars,
                "source_h5ad": portable_path(cfg, module2_h5ad_for(cfg, item["block_id"])),
            }
        )

    # shared gene intersection (prefer gene_id_no_version if present)
    gene_sets = []
    for a in objects:
        if "gene_id_no_version" in a.var.columns:
            gene_sets.append(set(a.var["gene_id_no_version"].astype(str)))
        else:
            gene_sets.append(set(pd.Index(a.var_names.astype(str)).str.replace(r"\.\d+$", "", regex=True)))

    shared = set.intersection(*gene_sets) if gene_sets else set()
    filtered = []
    for a in objects:
        if "gene_id_no_version" in a.var.columns:
            keep = a.var["gene_id_no_version"].astype(str).isin(shared)
            b = a[:, keep].copy()
            # unify var index to versionless Ensembl for joins
            b.var_names = b.var["gene_id_no_version"].astype(str).values
        else:
            vn = pd.Index(a.var_names.astype(str)).str.replace(r"\.\d+$", "", regex=True)
            a.var["gene_id_no_version"] = vn
            keep = a.var["gene_id_no_version"].isin(shared)
            b = a[:, keep].copy()
            b.var_names = b.var["gene_id_no_version"].astype(str).values
        b.var_names_make_unique()
        filtered.append(b)

    combined = ad.concat(filtered, join="inner", merge="same", label="batch_concat", keys=None)
    combined.obs_names_make_unique()
    summary = pd.DataFrame(summary_rows)
    # Keep uns JSON-serializable / h5ad-writable (no list-of-dicts)
    combined.uns["module3_concat"] = {
        "n_shared_genes": int(len(shared)),
        "n_donors": int(summary.shape[0]),
        "block_ids": ",".join(summary["block_id"].astype(str).tolist()),
        "donor_ids": ",".join(summary["donor_id"].astype(str).tolist()),
    }
    return combined, summary, sorted(shared)
