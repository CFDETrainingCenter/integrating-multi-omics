"""QC metrics, filtering, and filter logging for Module 2."""

from __future__ import annotations

from typing import Any

import pandas as pd


def compute_qc_metrics(adata):
    import scanpy as sc

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    # standardize convenience aliases used in plots/filters
    if "pct_counts_mt" in adata.obs.columns:
        adata.obs["percent_mito"] = adata.obs["pct_counts_mt"]
    if "pct_counts_ribo" in adata.obs.columns:
        adata.obs["percent_ribo"] = adata.obs["pct_counts_ribo"]
    return adata


def qc_metrics_summary(adata) -> pd.DataFrame:
    cols = [
        c
        for c in [
            "total_counts",
            "n_genes_by_counts",
            "percent_mito",
            "percent_ribo",
        ]
        if c in adata.obs.columns
    ]
    return adata.obs[cols].describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95]).T.reset_index(
        names="metric"
    )


def apply_filters(adata, cfg: dict[str, Any]):
    """Apply configured QC filters; return filtered adata + step log."""
    import scanpy as sc

    qc_cfg = cfg["module1"]["qc"]
    log_rows = [
        {
            "step": "loaded",
            "n_nuclei": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "rationale": "raw_expr.h5ad as QC starting matrix",
        }
    ]

    # Gene filter first (cells expressing too few genes handled next)
    sc.pp.filter_genes(adata, min_cells=int(qc_cfg["min_cells_per_gene"]))
    log_rows.append(
        {
            "step": "min_cells_per_gene",
            "n_nuclei": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "rationale": f"Retain genes detected in >= {qc_cfg['min_cells_per_gene']} nuclei",
        }
    )

    sc.pp.filter_cells(adata, min_genes=int(qc_cfg["min_genes"]))
    log_rows.append(
        {
            "step": "min_genes",
            "n_nuclei": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "rationale": f"Remove low-complexity nuclei with < {qc_cfg['min_genes']} genes",
        }
    )

    if "total_counts" in adata.obs.columns:
        keep = adata.obs["total_counts"] >= float(qc_cfg["min_counts"])
        adata = adata[keep].copy()
        log_rows.append(
            {
                "step": "min_counts",
                "n_nuclei": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "rationale": f"Remove nuclei with < {qc_cfg['min_counts']} total counts",
            }
        )

    if "percent_mito" in adata.obs.columns:
        keep = adata.obs["percent_mito"] <= float(qc_cfg["max_mito_percent"])
        adata = adata[keep].copy()
        log_rows.append(
            {
                "step": "max_mito_percent",
                "n_nuclei": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "rationale": (
                    f"Remove nuclei with mitochondrial percent > {qc_cfg['max_mito_percent']} "
                    "(snRNA typically low mito; threshold is diagnostic)"
                ),
            }
        )

    if "n_genes_by_counts" in adata.obs.columns:
        keep = adata.obs["n_genes_by_counts"] <= float(qc_cfg["max_genes"])
        adata = adata[keep].copy()
        log_rows.append(
            {
                "step": "max_genes",
                "n_nuclei": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "rationale": f"Remove extreme high-gene outliers > {qc_cfg['max_genes']}",
            }
        )

    log_rows.append(
        {
            "step": "final_qc",
            "n_nuclei": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "rationale": "Input for normalization / HVG / PCA / UMAP",
        }
    )
    return adata, pd.DataFrame(log_rows)
