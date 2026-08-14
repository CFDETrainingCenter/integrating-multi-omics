"""Donor-aware integration (Harmony) and embedding diagnostics."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from scripts.common.methods import module3_embedding_params


def recompute_pca_umap(adata, cfg: dict[str, Any], use_rep: str = "X_pca", key_added_umap: str = "X_umap"):
    import scanpy as sc

    emb = module3_embedding_params(cfg)
    seed = int(cfg["module2"].get("random_seed", 0))
    n_pcs = int(emb["n_pcs"])
    n_neighbors = int(emb["n_neighbors"])
    n_pcs_neighbors = int(emb.get("n_pcs_neighbors") or n_pcs)
    min_dist = float(emb["umap_min_dist"])
    spread = float(emb["umap_spread"])

    # HVG on combined object if needed
    if "highly_variable" not in adata.var.columns or int(adata.var["highly_variable"].sum()) < 100:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=int(emb["n_top_genes"]),
            flavor="seurat",
        )

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=float(emb.get("scale_max_value", 10)))
    sc.tl.pca(adata_hvg, n_comps=n_pcs, svd_solver="arpack")
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

    sc.pp.neighbors(
        adata,
        use_rep=use_rep,
        n_neighbors=n_neighbors,
        n_pcs=min(n_pcs_neighbors, adata.obsm[use_rep].shape[1]),
        random_state=seed,
    )
    sc.tl.umap(adata, min_dist=min_dist, spread=spread, random_state=seed)
    if key_added_umap != "X_umap":
        adata.obsm[key_added_umap] = adata.obsm["X_umap"].copy()
    adata.uns["module3_embedding_params"] = emb
    return adata


def run_harmony(adata, cfg: dict[str, Any]):
    """Run Harmony on PCA; store corrected embedding in obsm['X_pca_harmony']."""
    import scanpy as sc
    import harmonypy as hm

    emb = module3_embedding_params(cfg)
    key = cfg["module2"]["harmony"]["key"]
    max_iter = int(cfg["module2"]["harmony"].get("max_iter_harmony", 20))
    theta = cfg["module2"]["harmony"].get("theta")
    if key not in adata.obs.columns:
        raise KeyError(f"Harmony key '{key}' missing from adata.obs")

    if "X_pca" not in adata.obsm:
        adata = recompute_pca_umap(adata, cfg)

    # harmonypy expects cells x PCs matrix and a metadata dataframe
    meta = adata.obs[[key]].copy()
    seed = int(cfg["module2"].get("random_seed", 0))
    harmony_kwargs: dict[str, Any] = {
        "max_iter_harmony": max_iter,
        "random_state": seed,
    }
    if theta is not None:
        harmony_kwargs["theta"] = float(theta)

    ho = hm.run_harmony(
        np.asarray(adata.obsm["X_pca"]),
        meta,
        [key],
        **harmony_kwargs,
    )
    # Harmony returns corrected PCs; orientation can vary by version
    z = np.asarray(ho.Z_corr)
    if z.ndim != 2:
        raise ValueError(f"Unexpected Harmony Z_corr shape: {z.shape}")
    if z.shape[0] == adata.n_obs:
        adata.obsm["X_pca_harmony"] = z
    elif z.shape[1] == adata.n_obs:
        adata.obsm["X_pca_harmony"] = z.T
    else:
        raise ValueError(
            f"Harmony Z_corr shape {z.shape} incompatible with n_obs={adata.n_obs}"
        )

    # UMAP on harmony space - use Module 3 embedding knobs for cluster separation
    seed = int(cfg["module2"].get("random_seed", 0))
    n_neighbors = int(emb["n_neighbors"])
    n_pcs = adata.obsm["X_pca_harmony"].shape[1]
    n_pcs_neighbors = min(int(emb.get("n_pcs_neighbors") or n_pcs), n_pcs)
    sc.pp.neighbors(
        adata,
        use_rep="X_pca_harmony",
        n_neighbors=n_neighbors,
        n_pcs=n_pcs_neighbors,
        random_state=seed,
        key_added="harmony",
    )
    sc.tl.umap(
        adata,
        neighbors_key="harmony",
        min_dist=float(emb["umap_min_dist"]),
        spread=float(emb["umap_spread"]),
        random_state=seed,
    )
    adata.obsm["X_umap_harmony"] = adata.obsm["X_umap"].copy()
    adata.uns["harmony_params"] = {
        "key": key,
        "max_iter_harmony": max_iter,
        "theta": theta,
        "random_state": seed,
        "n_neighbors": n_neighbors,
        "n_pcs_neighbors": n_pcs_neighbors,
        "umap_min_dist": float(emb["umap_min_dist"]),
        "umap_spread": float(emb["umap_spread"]),
    }
    return adata


# Coarse compartments for teaching figures (Azimuth fine labels -> broad class).
# Four buckets only. Completeness of this dict is the teaching point: unmapped
# labels are named loudly and never silently bucketed to a catch-all.
_AZIMUTH_TO_COARSE: dict[str, str] = {
    "AT1": "epithelial",
    "AT2": "epithelial",
    "AT2 proliferating": "epithelial",
    "Basal resting": "epithelial",
    "Suprabasal": "epithelial",
    "Multiciliated (non-nasal)": "epithelial",
    "Ionocyte": "epithelial",
    "Tuft": "epithelial",
    # Former airway_secretory members folded into epithelium (P0-6).
    "Transitional Club-AT2": "epithelial",
    "Club (non-nasal)": "epithelial",
    "Club (nasal)": "epithelial",
    "Goblet (bronchial)": "epithelial",
    "Goblet (nasal)": "epithelial",
    "SMG duct": "epithelial",
    "SMG serous (bronchial)": "epithelial",
    "SMG mucous": "epithelial",
    "Deuterosomal": "epithelial",
    "EC general capillary": "endothelial",
    "EC aerocyte capillary": "endothelial",
    "EC arterial": "endothelial",
    "EC venous pulmonary": "endothelial",
    "EC venous systemic": "endothelial",
    "Lymphatic EC mature": "endothelial",
    "Lymphatic EC differentiating": "endothelial",
    "Alveolar fibroblasts": "stromal",
    "Adventitial fibroblasts": "stromal",
    "Peribronchial fibroblasts": "stromal",
    "Myofibroblasts": "stromal",
    "Pericytes": "stromal",
    "Smooth muscle": "stromal",
    "SM activated stress response": "stromal",
    "Alveolar macrophages": "immune",
    "Interstitial Mφ perivascular": "immune",
    "Interstitial Mphi perivascular": "immune",
    "Monocyte-derived Mφ": "immune",
    "Monocyte-derived Mphi": "immune",
    "Classical monocytes": "immune",
    "Non-classical monocytes": "immune",
    "CD4 T cells": "immune",
    "CD8 T cells": "immune",
    "T cells proliferating": "immune",
    "NK cells": "immune",
    "B cells": "immune",
    "Plasma cells": "immune",
    "Mast cells": "immune",
    "DC2": "immune",
    "Plasmacytoid DCs": "immune",
}

_VALID_CELL_CLASSES = frozenset({"epithelial", "endothelial", "stromal", "immune"})


def assert_azimuth_mapping_coverage(labels: list[str] | set[str] | pd.Series) -> list[str]:
    """Return unmapped labels (empty list = full coverage). Never silently buckets."""
    present = sorted({str(x) for x in labels if str(x) not in {"", "nan", "None", "unlabeled"}})
    return [lab for lab in present if lab not in _AZIMUTH_TO_COARSE]


def cell_class_membership_frame(
    labels: pd.Series,
    *,
    source_col: str = "azimuth_label",
) -> pd.DataFrame:
    """Per-compartment membership table (fine label counts inside each cell_class)."""
    raw = labels.astype("string").fillna("unlabeled")
    rows: list[dict[str, Any]] = []
    for lab, n in raw.value_counts().items():
        s = str(lab)
        if s == "unlabeled" or s.lower() in {"nan", "none", ""}:
            compartment = "unlabeled"
        elif s in _AZIMUTH_TO_COARSE:
            compartment = _AZIMUTH_TO_COARSE[s]
        else:
            compartment = f"unmapped:{s}"
        rows.append(
            {
                source_col: s,
                "cell_class": compartment,
                "n_nuclei": int(n),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["fraction_of_compartment"] = out.groupby("cell_class")["n_nuclei"].transform(
        lambda s: s / s.sum() if float(s.sum()) else 0.0
    )
    return out.sort_values(["cell_class", "n_nuclei"], ascending=[True, False]).reset_index(drop=True)


def assign_coarse_cell_class(
    adata,
    source_col: str = "azimuth_label",
    key_added: str = "cell_class",
    *,
    on_unmapped: str = "warn",
) -> str:
    """
    Map fine Azimuth labels to broad compartments.

    ``unlabeled`` stays its own category (Azimuth had nothing to say).
    Unmapped labels never silently become a biological class: they warn (default)
    or raise (``on_unmapped='raise'``) and are stored as ``unmapped:<label>``.
    """
    if source_col not in adata.obs.columns:
        adata.obs[key_added] = "unlabeled"
        return key_added

    raw = adata.obs[source_col].astype("string").fillna("unlabeled")
    present = sorted({str(x) for x in raw.unique()})
    unmapped = assert_azimuth_mapping_coverage(present)

    if unmapped:
        msg = (
            "Azimuth labels with no coarse mapping (not silently bucketed to 'other'): "
            + ", ".join(unmapped)
        )
        if on_unmapped == "raise":
            raise KeyError(msg)
        print(f"WARNING: {msg}")

    def _map_one(x: str) -> str:
        s = str(x)
        if s == "unlabeled" or s.lower() in {"nan", "none", ""}:
            return "unlabeled"
        if s in _AZIMUTH_TO_COARSE:
            return _AZIMUTH_TO_COARSE[s]
        return f"unmapped:{s}"

    mapped = raw.map(_map_one)
    adata.obs[key_added] = pd.Categorical(mapped)
    adata.uns["cell_class_unmapped_labels"] = unmapped
    adata.uns["cell_class_mapping_n"] = {
        "n_labels_present": len(present),
        "n_mapped": len(present) - len(unmapped) - (1 if "unlabeled" in present else 0),
        "n_unmapped": len(unmapped),
    }
    adata.uns["cell_class_membership"] = cell_class_membership_frame(
        raw, source_col=source_col
    )
    return key_added


def compute_leiden(adata, cfg: dict[str, Any], key_added: str = "leiden") -> str:
    """Always compute Leiden on the current neighbor graph (Harmony if available)."""
    import scanpy as sc

    emb = module3_embedding_params(cfg)
    seed = int(cfg["module2"].get("random_seed", 0))
    resolution = float(emb.get("leiden_resolution", 0.8))
    neighbors_key = "harmony" if "harmony" in adata.uns else None
    if neighbors_key:
        sc.tl.leiden(
            adata,
            neighbors_key=neighbors_key,
            key_added=key_added,
            resolution=resolution,
            random_state=seed,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
    else:
        sc.tl.leiden(
            adata,
            key_added=key_added,
            resolution=resolution,
            random_state=seed,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
    return key_added


def ensure_cluster_labels(adata, cfg: dict[str, Any], label_col: str | None = None) -> str:
    """Return a usable categorical label column, computing Leiden if needed."""
    preferred = label_col or cfg["module2"].get("label_column", "azimuth_label")
    if preferred in adata.obs.columns and adata.obs[preferred].notna().any():
        # fill missing with unknown for cross-donor consistency
        adata.obs[preferred] = adata.obs[preferred].astype("string").fillna("unlabeled")
        return preferred

    compute_leiden(adata, cfg, key_added="leiden")
    adata.obs["cell_label"] = adata.obs["leiden"].astype(str)
    return "cell_label"


def post_harmony_label_views(adata, cfg: dict[str, Any]) -> list[str]:
    """
    Prepare comparison obs columns and return preferred UMAP color order:
    donor | azimuth | leiden | coarse cell class.
    """
    label_col = cfg["module2"].get("label_column", "azimuth_label")
    if label_col in adata.obs.columns:
        adata.obs[label_col] = adata.obs[label_col].astype("string").fillna("unlabeled")
    compute_leiden(adata, cfg, key_added="leiden")
    assign_coarse_cell_class(adata, source_col=label_col if label_col in adata.obs.columns else "leiden")
    colors = [
        c
        for c in ["donor_label", label_col, "leiden", "cell_class"]
        if c in adata.obs.columns
    ]
    return colors


def donor_composition(adata, label_col: str = "azimuth_label") -> pd.DataFrame:
    """Per-donor nucleus counts by label (counts, not only fractions)."""
    if "donor_label" not in adata.obs.columns and "donor_id" not in adata.obs.columns:
        return pd.DataFrame()
    donor_col = "donor_label" if "donor_label" in adata.obs.columns else "donor_id"
    if label_col not in adata.obs.columns:
        label_col = "leiden" if "leiden" in adata.obs.columns else None
    if label_col is None or label_col not in adata.obs.columns:
        return (
            adata.obs[donor_col]
            .value_counts()
            .rename_axis(donor_col)
            .reset_index(name="n_nuclei")
        )
    ct = pd.crosstab(adata.obs[donor_col], adata.obs[label_col])
    return ct.reset_index()


def donor_composition_fractions(adata, label_col: str = "azimuth_label") -> pd.DataFrame:
    """Row-normalized donor x label fractions."""
    if "donor_label" not in adata.obs.columns and "donor_id" not in adata.obs.columns:
        return pd.DataFrame()
    donor_col = "donor_label" if "donor_label" in adata.obs.columns else "donor_id"
    if label_col not in adata.obs.columns:
        return pd.DataFrame()
    ct = pd.crosstab(adata.obs[donor_col], adata.obs[label_col], normalize="index")
    return ct.reset_index()


def harmony_mixing_metrics(
    adata,
    *,
    donor_key: str = "donor_id",
    k: int = 30,
    max_cells: int = 8000,
    random_seed: int = 0,
) -> pd.DataFrame:
    """
    Quantify Harmony effect: mean fraction of kNN from a *different* donor
    on X_pca vs X_pca_harmony (higher = more mixed).
    """
    from sklearn.neighbors import NearestNeighbors

    if donor_key not in adata.obs.columns:
        donor_key = "donor_label" if "donor_label" in adata.obs.columns else None
    if donor_key is None:
        raise KeyError("No donor key for Harmony mixing metrics")
    if "X_pca" not in adata.obsm or "X_pca_harmony" not in adata.obsm:
        raise KeyError("Need X_pca and X_pca_harmony")

    rng = np.random.default_rng(random_seed)
    n = adata.n_obs
    if n > max_cells:
        idx = np.sort(rng.choice(n, size=max_cells, replace=False))
    else:
        idx = np.arange(n)
    donors = adata.obs[donor_key].astype(str).to_numpy()[idx]

    rows = []
    for rep in ("X_pca", "X_pca_harmony"):
        X = np.asarray(adata.obsm[rep][idx], dtype=np.float64)
        nn = NearestNeighbors(n_neighbors=min(k + 1, len(idx)), algorithm="auto")
        nn.fit(X)
        ind = nn.kneighbors(return_distance=False)[:, 1:]  # drop self
        frac = []
        for i, neigh in enumerate(ind):
            frac.append(float(np.mean(donors[neigh] != donors[i])))
        rows.append(
            {
                "embedding": rep,
                "k": int(ind.shape[1]),
                "n_cells_scored": int(len(idx)),
                "mean_foreign_donor_knn_fraction": float(np.mean(frac)),
                "median_foreign_donor_knn_fraction": float(np.median(frac)),
            }
        )
    return pd.DataFrame(rows)


def harmony_silhouette_by_celltype(
    adata,
    *,
    label_col: str = "cell_class",
    donor_key: str = "donor_id",
    max_cells: int = 8000,
    random_seed: int = 0,
) -> pd.DataFrame:
    """
    Per-label silhouette of donor identity on X_pca vs X_pca_harmony.
    Lower silhouette after Harmony => donors less separable within that cell class.
    """
    from sklearn.metrics import silhouette_score

    if label_col not in adata.obs.columns:
        return pd.DataFrame()
    if donor_key not in adata.obs.columns:
        donor_key = "donor_label" if "donor_label" in adata.obs.columns else None
    if donor_key is None:
        return pd.DataFrame()

    rng = np.random.default_rng(random_seed)
    rows = []
    for lab in sorted(adata.obs[label_col].astype(str).unique()):
        mask = (adata.obs[label_col].astype(str) == lab).to_numpy()
        idxs = np.where(mask)[0]
        if len(idxs) < 30:
            continue
        if len(idxs) > max_cells:
            idxs = np.sort(rng.choice(idxs, size=max_cells, replace=False))
        donors = adata.obs[donor_key].astype(str).to_numpy()[idxs]
        if len(set(donors)) < 2:
            continue
        for rep in ("X_pca", "X_pca_harmony"):
            if rep not in adata.obsm:
                continue
            X = np.asarray(adata.obsm[rep][idxs], dtype=np.float64)
            try:
                sil = float(silhouette_score(X, donors, metric="euclidean"))
            except Exception:  # noqa: BLE001
                sil = float("nan")
            rows.append(
                {
                    "label": lab,
                    "embedding": rep,
                    "n_cells": int(len(idxs)),
                    "n_donors": int(len(set(donors))),
                    "silhouette_donor": sil,
                }
            )
    return pd.DataFrame(rows)


def leiden_resolution_sweep(
    adata,
    cfg: dict[str, Any],
    resolutions: list[float] | None = None,
) -> pd.DataFrame:
    """Compute Leiden at a few resolutions; store columns leiden_r{res}."""
    import scanpy as sc

    emb = module3_embedding_params(cfg)
    seed = int(cfg["module2"].get("random_seed", 0))
    resolutions = resolutions or list(
        cfg.get("module2", {}).get("embedding", {}).get("leiden_resolutions")
        or [0.4, 0.8, 1.2]
    )
    neighbors_key = "harmony" if "harmony" in adata.uns else None
    rows = []
    for res in resolutions:
        key = f"leiden_r{res}"
        kwargs = dict(
            key_added=key,
            resolution=float(res),
            random_state=seed,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
        if neighbors_key:
            sc.tl.leiden(adata, neighbors_key=neighbors_key, **kwargs)
        else:
            sc.tl.leiden(adata, **kwargs)
        n = int(adata.obs[key].nunique())
        rows.append({"resolution": float(res), "key": key, "n_clusters": n})
    # keep default leiden at configured resolution for downstream
    default_res = float(emb.get("leiden_resolution", 0.8))
    default_key = f"leiden_r{default_res}"
    if default_key in adata.obs.columns:
        adata.obs["leiden"] = adata.obs[default_key].astype(str)
    return pd.DataFrame(rows)
