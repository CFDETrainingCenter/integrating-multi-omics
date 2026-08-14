"""Normalization, HVG, PCA, neighbors, UMAP for Module 2."""

from __future__ import annotations

from typing import Any

import numpy as np

from scripts.common.methods import module1_embedding_params


def normalize_log(adata, cfg: dict[str, Any]):
    import scanpy as sc

    # Preserve counts explicitly for downstream DE
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    pp = module1_embedding_params(cfg)
    target = float(pp["target_sum"])
    sc.pp.normalize_total(adata, target_sum=target)
    sc.pp.log1p(adata)
    adata.uns["preprocessing"] = {
        "normalized": True,
        "log1p": True,
        "target_sum": target,
        "counts_layer": "counts",
        "doublet_filtering": bool(cfg["module1"].get("doublet_filtering", False)),
    }
    return adata


def hvg_pca(adata, cfg: dict[str, Any]):
    import scanpy as sc

    pp = module1_embedding_params(cfg)
    # Use seurat flavor on log-normalized X to avoid seurat_v3/skmisc dependency
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=int(pp["n_top_genes"]),
        flavor="seurat",
        subset=False,
    )
    adata.var["highly_variable"] = adata.var["highly_variable"].astype(bool)

    # PCA on HVG subspace without permanently subsetting genes
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=float(pp.get("scale_max_value", 10)))
    sc.tl.pca(adata_hvg, n_comps=int(pp["n_pcs"]), svd_solver="arpack")
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.uns["pca"] = adata_hvg.uns.get("pca", {})
    adata.uns["embedding_params"] = {
        "n_top_genes": int(pp["n_top_genes"]),
        "n_pcs": int(pp["n_pcs"]),
        "scale_max_value": float(pp.get("scale_max_value", 10)),
        "hvg_flavor": "seurat",
        "svd_solver": "arpack",
    }
    return adata


def neighbors_umap(adata, cfg: dict[str, Any]):
    import scanpy as sc

    pp = module1_embedding_params(cfg)
    seed = int(cfg["module1"].get("random_seed", 0))
    n_pcs = int(pp["n_pcs"])
    n_pcs_use = int(pp.get("n_pcs_neighbors") or n_pcs)
    sc.pp.neighbors(
        adata,
        n_neighbors=int(pp["n_neighbors"]),
        n_pcs=min(n_pcs_use, adata.obsm["X_pca"].shape[1]),
        random_state=seed,
    )
    sc.tl.umap(
        adata,
        min_dist=float(pp.get("umap_min_dist", 0.3)),
        spread=float(pp.get("umap_spread", 1.0)),
        random_state=seed,
    )
    adata.uns["umap_params"] = {
        "n_neighbors": int(pp["n_neighbors"]),
        "n_pcs": min(n_pcs_use, adata.obsm["X_pca"].shape[1]),
        "min_dist": float(pp.get("umap_min_dist", 0.3)),
        "spread": float(pp.get("umap_spread", 1.0)),
        "random_state": seed,
    }
    return adata


def maybe_subsample(adata, cfg: dict[str, Any], random_seed: int | None = None):
    """Optional teaching-size cap after QC.

    Returns ``(adata, info)`` where ``info`` has enough fields for a filter-log row
    (``step``, ``n_nuclei``, ``rationale``, ``seed``, ``max_nuclei``).
    """
    max_n = cfg["module1"].get("max_nuclei")
    seed = int(cfg["module1"].get("random_seed", 0) if random_seed is None else random_seed)
    n_before = int(adata.n_obs)
    if max_n is None or n_before <= int(max_n):
        return adata, {
            "subsampled": False,
            "n_obs": n_before,
            "n_obs_before": n_before,
            "n_obs_after": n_before,
            "n_nuclei": n_before,
            "max_nuclei": None if max_n is None else int(max_n),
            "seed": seed,
            "step": "teaching_subsample",
            "rationale": (
                "no teaching subsample "
                f"(n_nuclei={n_before}"
                + (f" <= module1.max_nuclei={int(max_n)}" if max_n is not None else "")
                + ")"
            ),
        }

    rng = np.random.default_rng(seed)
    idx = rng.choice(adata.n_obs, size=int(max_n), replace=False)
    idx = np.sort(idx)
    out = adata[idx].copy()
    n_after = int(out.n_obs)
    return out, {
        "subsampled": True,
        "n_obs_before": n_before,
        "n_obs_after": n_after,
        "n_nuclei": n_after,
        "max_nuclei": int(max_n),
        "seed": seed,
        "step": "teaching_subsample",
        "rationale": f"capped at module1.max_nuclei={int(max_n)}, seed={seed}",
    }


def append_subsample_filter_row(filter_log, sub_info: dict[str, Any], adata):
    """Append a teaching_subsample row to a QC filter log DataFrame."""
    import pandas as pd

    row = {
        "step": str(sub_info.get("step", "teaching_subsample")),
        "n_nuclei": int(sub_info.get("n_nuclei", adata.n_obs)),
        "n_genes": int(adata.n_vars),
        "rationale": str(sub_info.get("rationale", "")),
    }
    return pd.concat([filter_log, pd.DataFrame([row])], ignore_index=True)
