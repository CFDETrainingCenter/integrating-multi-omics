"""Load configured Muon/MOFA multiome assets without materializing the full MuData."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd

from scripts.common.paths import module_root
from scripts.nb03_integration.gtex import strip_ensembl_version


def multiome_dir(cfg: dict[str, Any]) -> Path:
    return module_root(cfg) / cfg["module3"]["multiome_dir"]


def multiome_paths(cfg: dict[str, Any]) -> dict[str, Path]:
    root = multiome_dir(cfg)
    files = cfg["module3"]["files"]
    return {
        "dir": root,
        "secondary_h5mu": root / files["secondary_h5mu"],
        "mudata_raw_h5mu": root / files.get("mudata_raw_h5mu", "mudata_raw.h5mu"),
        "mofa_hdf5": root / files["mofa_hdf5"],
        "pdf_combined": root / files["pdf_combined"],
        "pdf_rna": root / files["pdf_rna"],
        "pdf_atac": root / files["pdf_atac"],
        "label_cache": module_root(cfg) / files["label_cache"]
        if files.get("label_cache")
        else None,
    }


def _decode_arr(values) -> np.ndarray:
    out = []
    for v in values:
        if isinstance(v, bytes):
            out.append(v.decode("utf-8", errors="replace"))
        else:
            out.append(str(v))
    return np.asarray(out, dtype=object)


def _read_obs_column(group: h5py.Group, key: str) -> pd.Series | None:
    if key not in group:
        return None
    obj = group[key]
    if isinstance(obj, h5py.Dataset):
        vals = obj[:]
        if vals.dtype.kind in ("S", "O"):
            vals = _decode_arr(vals)
        return pd.Series(vals)
    # categorical encoding used by anndata
    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        cats = _decode_arr(obj["categories"][:])
        codes = np.asarray(obj["codes"][:])
        # anndata uses -1 for NA
        labels = np.array([cats[c] if c >= 0 else None for c in codes], dtype=object)
        return pd.Series(labels)
    return None


_METADATA_ALIASES = {
    "rna:leiden": "rna_leiden",
    "atac_cbg:leiden": "atac_leiden",
    "atac_cbg:Clusters": "atac_clusters",
    "atac_cbg:nFrags": "nFrags",
    "atac_cbg:ReadsInTSS": "ReadsInTSS",
    "atac_cbg:ReadsInPromoter": "ReadsInPromoter",
    "atac_cbg:PromoterRatio": "PromoterRatio",
    "atac_cbg:DoubletEnrichment": "DoubletEnrichment",
    "atac_cbg:DoubletScore": "DoubletScore",
    "atac_cbg:TSSEnrichment": "TSSEnrichment",
    "atac_cbg:BlacklistRatio": "BlacklistRatio",
    "atac_cbg:NucleosomeRatio": "NucleosomeRatio",
    "atac_cbg:PassQC": "PassQC",
    "atac_cbg:Sample": "Sample",
    "atac_cbg:nDiFrags": "nDiFrags",
    "atac_cbg:nMonoFrags": "nMonoFrags",
    "atac_cbg:nMultiFrags": "nMultiFrags",
    "atac_cbg:ReadsInBlacklist": "ReadsInBlacklist",
}


def mofa_sample_metadata(cfg: dict[str, Any]) -> pd.DataFrame:
    """Read every ``samples_metadata/group1/*`` key, indexed by MOFA sample order."""
    paths = multiome_paths(cfg)
    mofa = paths["mofa_hdf5"]
    if not mofa.exists():
        raise FileNotFoundError(f"Missing MOFA file: {mofa}")
    with h5py.File(mofa, "r") as f:
        samples = _decode_arr(f["samples/group1"][:]).astype(str)
        meta_grp = f["samples_metadata/group1"]
        data: dict[str, np.ndarray] = {}
        for key in meta_grp.keys():
            obj = meta_grp[key]
            if not isinstance(obj, h5py.Dataset) or obj.ndim != 1:
                continue
            vals = obj[:]
            if vals.dtype.kind in ("S", "O"):
                data[key] = _decode_arr(vals)
            else:
                data[key] = np.asarray(vals)
    frame = pd.DataFrame(data, index=pd.Index(samples, name="barcode"))
    for src, dst in _METADATA_ALIASES.items():
        if src in frame.columns and dst not in frame.columns:
            frame[dst] = frame[src]
    return frame


def load_label_cache(cfg: dict[str, Any]) -> pd.DataFrame:
    """Load the committed Azimuth / WNN cache, aligned to MOFA sample order."""
    paths = multiome_paths(cfg)
    cache_path = paths.get("label_cache")
    if cache_path is None:
        raise FileNotFoundError(
            "module3.files.label_cache is not set. Rebuild the cache with "
            "scripts.nb05_multimodal.cache.build_label_cache."
        )
    if not cache_path.exists():
        raise FileNotFoundError(
            f"Missing Module 3 label cache: {cache_path}. "
            "Generate it with scripts.nb05_multimodal.cache.build_label_cache "
            "(authoring machine only; learners receive the committed TSV)."
        )
    mofa = paths["mofa_hdf5"]
    if not mofa.exists():
        raise FileNotFoundError(f"Missing MOFA file: {mofa}")
    with h5py.File(mofa, "r") as f:
        samples = _decode_arr(f["samples/group1"][:]).astype(str)
    cache = pd.read_csv(cache_path, sep="\t", comment="#")
    if "barcode" not in cache.columns:
        raise ValueError(f"{cache_path} has no barcode column")
    cache["barcode"] = cache["barcode"].astype(str)
    cache_set = set(cache["barcode"].tolist())
    mofa_set = set(samples.tolist())
    if cache_set != mofa_set:
        raise ValueError(
            "Label cache barcodes do not match the MOFA sample list. "
            f"only_mofa={sorted(mofa_set - cache_set)[:5]} "
            f"only_cache={sorted(cache_set - mofa_set)[:5]}"
        )
    cache = cache.set_index("barcode").loc[samples]
    cache.index.name = "barcode"
    return cache


def assemble_nucleus_frame(cfg: dict[str, Any]) -> pd.DataFrame:
    """Join MOFA sample metadata and the committed label cache on barcode."""
    meta = mofa_sample_metadata(cfg)
    cache = load_label_cache(cfg)
    if not meta.index.equals(cache.index):
        cache = cache.reindex(meta.index)
    return meta.join(cache, how="left")


def modality_shapes(cfg: dict[str, Any]) -> pd.DataFrame:
    """Lightweight shape/provenance table from secondary_analysis.h5mu."""
    paths = multiome_paths(cfg)
    h5mu = paths["secondary_h5mu"]
    if not h5mu.exists():
        raise FileNotFoundError(f"Missing MuData: {h5mu}")
    rows = []
    with h5py.File(h5mu, "r") as f:
        for mod in f["mod"].keys():
            g = f["mod"][mod]
            n_obs = int(g["obs"]["_index"].shape[0]) if "_index" in g["obs"] else None
            if "X" in g and isinstance(g["X"], h5py.Dataset):
                n_vars = int(g["X"].shape[1])
                x_type = "dense"
            elif "X" in g and isinstance(g["X"], h5py.Group) and "shape" in g["X"].attrs:
                shape = list(g["X"].attrs["shape"])
                n_vars = int(shape[1])
                x_type = "sparse"
            else:
                n_vars = int(g["var"]["_index"].shape[0]) if "_index" in g.get("var", {}) else None
                x_type = "unknown"
            rows.append(
                {
                    "modality": mod,
                    "n_nuclei": n_obs,
                    "n_features": n_vars,
                    "matrix_type": x_type,
                    "has_umap": "X_umap" in g.get("obsm", {}),
                    "has_pca": "X_pca" in g.get("obsm", {}),
                }
            )
    return pd.DataFrame(rows)


def barcode_sets(cfg: dict[str, Any]) -> dict[str, set[str]]:
    paths = multiome_paths(cfg)
    out: dict[str, set[str]] = {}
    with h5py.File(paths["secondary_h5mu"], "r") as f:
        for mod in f["mod"].keys():
            idx = _decode_arr(f[f"mod/{mod}/obs/_index"][:])
            out[mod] = set(map(str, idx))
    return out


def load_joint_embedding(cfg: dict[str, Any]) -> ad.AnnData:
    """
    Build a lightweight AnnData with joint MOFA/UMAP embeddings + RNA labels.

    Avoids loading modality count matrices from the 10 GB MuData.
    """
    paths = multiome_paths(cfg)
    h5mu = paths["secondary_h5mu"]
    if not h5mu.exists():
        raise FileNotFoundError(f"Missing MuData: {h5mu}")

    with h5py.File(h5mu, "r") as f:
        barcodes = _decode_arr(f["obs"]["_index"][:]).astype(str)
        obsm = {}
        if "X_mofa" in f["obsm"]:
            obsm["X_mofa"] = np.asarray(f["obsm"]["X_mofa"][:], dtype=np.float64)
        if "X_umap" in f["obsm"]:
            obsm["X_umap"] = np.asarray(f["obsm"]["X_umap"][:], dtype=np.float32)

        # RNA modality labels / UMAP (may be in different barcode order)
        rna_bc = _decode_arr(f["mod/rna/obs/_index"][:]).astype(str)
        rna_pos = {b: i for i, b in enumerate(rna_bc)}
        order = np.array([rna_pos[b] for b in barcodes])

        obs = pd.DataFrame(index=pd.Index(barcodes, name="barcode"))
        for col in ("azimuth_label", "predicted_label", "leiden", "prediction_score"):
            s = _read_obs_column(f["mod/rna/obs"], col)
            if s is not None:
                obs[col] = s.to_numpy()[order]

        # ATAC QC / cluster aligned to joint barcodes
        atac_bc = _decode_arr(f["mod/atac_cbg/obs/_index"][:]).astype(str)
        atac_pos = {b: i for i, b in enumerate(atac_bc)}
        atac_order = np.array([atac_pos[b] for b in barcodes])
        for col, rename in (
            ("Clusters", "atac_cluster"),
            ("leiden", "atac_leiden"),
            ("TSSEnrichment", "tss_enrichment"),
            ("nFrags", "n_frags"),
        ):
            s = _read_obs_column(f["mod/atac_cbg/obs"], col)
            if s is not None:
                obs[rename] = s.to_numpy()[atac_order]

        if "X_umap" in f["mod/rna/obsm"]:
            obsm["X_umap_rna"] = np.asarray(f["mod/rna/obsm"]["X_umap"][:], dtype=np.float32)[order]
        if "X_umap" in f["mod/atac_cbg/obsm"]:
            obsm["X_umap_atac"] = np.asarray(f["mod/atac_cbg/obsm"]["X_umap"][:], dtype=np.float32)[
                atac_order
            ]

    adata = ad.AnnData(obs=obs)
    for k, v in obsm.items():
        adata.obsm[k] = v
    adata.uns["hubmap_donor_id"] = cfg["module3"]["donor_id"]
    adata.uns["hubmap_dataset_id"] = cfg["module3"]["dataset_id"]
    adata.uns["module3_source"] = str(h5mu)
    return adata


def _column_means(dset: h5py.Dataset, chunk: int = 2000) -> np.ndarray:
    n, p = dset.shape
    acc = np.zeros(p, dtype=np.float64)
    for i in range(0, n, chunk):
        acc += np.asarray(dset[i : i + chunk], dtype=np.float64).sum(axis=0)
    return acc / float(n)


def load_mofa(cfg: dict[str, Any]) -> dict[str, Any]:
    """Load MOFA factors, variance explained, and feature annotations."""
    paths = multiome_paths(cfg)
    mofa = paths["mofa_hdf5"]
    if not mofa.exists():
        raise FileNotFoundError(f"Missing MOFA file: {mofa}")

    with h5py.File(mofa, "r") as f:
        samples = _decode_arr(f["samples/group1"][:]).astype(str)
        z = np.asarray(f["expectations/Z/group1"][:], dtype=np.float64).T  # nuclei x factors
        views = _decode_arr(f["views/views"][:]).astype(str).tolist()
        r2_total = np.asarray(f["variance_explained/r2_total/group1"][:], dtype=np.float64)
        r2_per_factor = np.asarray(
            f["variance_explained/r2_per_factor/group1"][:], dtype=np.float64
        )  # views x factors

        features = {
            "rna": _decode_arr(f["features/rna"][:]).astype(str),
            "atac_cbg": _decode_arr(f["features/atac_cbg"][:]).astype(str),
        }
        hugo = _decode_arr(f["features_metadata/rna/hugo_symbol"][:]).astype(str)

        # streaming means + std for gene bridge (MOFA inputs are centered; use std)
        rna_means = _column_means(f["data/rna/group1"])
        atac_means = _column_means(f["data/atac_cbg/group1"])
        rna_std = _column_std(f["data/rna/group1"], rna_means)
        atac_std = _column_std(f["data/atac_cbg/group1"], atac_means)

    factors = pd.DataFrame(
        z,
        index=pd.Index(samples, name="barcode"),
        columns=[f"Factor{i+1}" for i in range(z.shape[1])],
    )
    # MOFA stores r2 on a percent scale here, so 1.03 means 1.03 percent, not 103
    # percent. The units column ships with the table so the value cannot be misread
    # as a fraction (an R2 above 1 would otherwise look impossible).
    var_df = pd.DataFrame(
        {
            "view": views,
            "r2_total": r2_total,
            "units": "percent",
        }
    )
    r2_factor_df = pd.DataFrame(
        r2_per_factor.T,
        columns=views,
        index=[f"Factor{i+1}" for i in range(r2_per_factor.shape[1])],
    )
    return {
        "path": mofa,
        "factors": factors,
        "views": views,
        "variance_total": var_df,
        "variance_per_factor": r2_factor_df,
        "features": features,
        "rna_hugo": hugo,
        "rna_means": rna_means,
        "atac_means": atac_means,
        "rna_std": rna_std,
        "atac_std": atac_std,
    }


def recompute_mofa_r2_total(cfg: dict[str, Any], view: str = "rna") -> dict[str, Any]:
    """
    Recompute total R^2 for one view from Z, W, and Y in multiome_mofa.hdf5.

    Removes unsupported 'scaling artifact' claims unless recomputation matches stored value.
    """
    paths = multiome_paths(cfg)
    mofa = paths["mofa_hdf5"]
    with h5py.File(mofa, "r") as f:
        views = _decode_arr(f["views/views"][:]).astype(str).tolist()
        if view not in views:
            raise KeyError(f"view {view} not in {views}")
        v_idx = views.index(view)
        Z = np.asarray(f["expectations/Z/group1"][:], dtype=np.float64).T  # n x k
        W = np.asarray(f[f"expectations/W/{view}"][:], dtype=np.float64)
        # W is typically features x factors or factors x features
        if W.shape[0] == Z.shape[1]:
            W_use = W.T
        elif W.shape[1] == Z.shape[1]:
            W_use = W
        else:
            return {
                "ok": False,
                "view": view,
                "error": f"W shape {W.shape} incompatible with Z factors={Z.shape[1]}",
            }
        Y = np.asarray(f[f"data/{view}/group1"][:], dtype=np.float64)
        stored = float(np.asarray(f["variance_explained/r2_total/group1"][:])[v_idx])
        # Match nuclei x features
        if Y.shape[0] != Z.shape[0] and Y.shape[1] == Z.shape[0]:
            Y = Y.T
        pred = Z @ W_use.T if W_use.shape[0] == Y.shape[1] else Z @ W_use
        if pred.shape != Y.shape:
            # try alternate
            pred = Z @ W_use.T
        if pred.shape != Y.shape:
            return {
                "ok": False,
                "view": view,
                "error": f"pred shape {pred.shape} != Y shape {Y.shape}",
                "W_shape": list(W.shape),
                "Z_shape": list(Z.shape),
                "Y_shape": list(Y.shape),
            }
        ss_res = float(np.nansum((Y - pred) ** 2))
        ss_tot = float(np.nansum(Y**2))  # MOFA uses centered Y; SS around 0
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        n_nan = int(np.isnan(Y).sum())
        matches_direct = abs(r2 - stored) < 1e-3
        matches_percent = abs(r2 * 100 - stored) < 1e-3
        if matches_direct and matches_percent:
            scale_branch = "both"
        elif matches_direct:
            scale_branch = "direct"
        elif matches_percent:
            scale_branch = "percent"
        else:
            scale_branch = "neither"
        matches = matches_direct or matches_percent
        if stored > 1.0 and matches_percent and not matches_direct:
            label = (
                "Recomputation matches stored/100 (percent scale); "
                "Gate 2 confirms convention before swapping published %."
            )
        elif stored > 1.0 and matches_direct:
            label = (
                "Recomputation matches stored value (>1); treat as model/scaling convention, not % variance."
            )
        elif stored > 1.0 and not matches:
            label = (
                "Stored R2 > 1 is NOT reproduced by Z@W reconstruction "
                "- do not call it a scaling artifact."
            )
        elif matches:
            label = "Recomputation matches stored R2."
        else:
            label = "Recomputation differs from stored R2; report both."
        return {
            "ok": True,
            "view": view,
            "rna_r2_recomputed" if view == "rna" else "r2_recomputed": r2,
            "rna_r2_stored" if view == "rna" else "r2_stored": stored,
            "ss_res": ss_res,
            "ss_tot": ss_tot,
            "n_nan_in_Y": n_nan,
            "matches_stored": matches,
            "matches_direct": bool(matches_direct),
            "matches_percent": bool(matches_percent),
            "scale_branch": scale_branch,
            "label": label,
            "W_shape": list(W.shape),
            "Z_shape": list(Z.shape),
            "Y_shape": list(Y.shape),
        }


def mofa_nondominant_r2_summary(
    variance_per_factor: pd.DataFrame,
    *,
    scale_branch: str | None = None,
    r2_recompute: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Max non-dominant R^2 and which factors are numerically dead.

    ``scale_branch`` must come from ``recompute_mofa_r2_total`` (``percent`` or
    ``direct``). Stored per-factor values share units with stored totals. When the
    branch is ``percent``, ``max_nondominant_percent`` equals the stored table
    value (do **not** multiply by 100).
    """
    if variance_per_factor is None or variance_per_factor.empty:
        return {}
    arr = variance_per_factor.to_numpy(dtype=float)
    nondom = []
    for i, row in enumerate(arr):
        order = np.argsort(row)
        nondom.append(float(row[order[-2]]) if len(row) > 1 else 0.0)
    max_nd = float(max(nondom)) if nondom else 0.0
    # dead: both views tiny
    dead = [
        str(variance_per_factor.index[i])
        for i, row in enumerate(arr)
        if float(np.max(row)) < 1e-3
    ]

    branch = scale_branch or (r2_recompute or {}).get("scale_branch")
    if branch == "both":
        stored = (r2_recompute or {}).get("rna_r2_stored")
        branch = "percent" if stored is not None and float(stored) > 1.0 else "direct"
    if branch == "percent":
        max_pct = max_nd
        max_frac = max_nd / 100.0
        units = "percent"
    elif branch == "direct":
        max_pct = 100.0 * max_nd
        max_frac = max_nd
        units = "fraction"
    else:
        raise ValueError(
            "mofa_nondominant_r2_summary requires scale_branch from "
            f"recompute_mofa_r2_total (got {branch!r}). Pass r2_recompute=..."
        )

    return {
        "max_nondominant_r2": max_nd,
        "max_nondominant_percent": max_pct,
        "max_nondominant_fraction": max_frac,
        "r2_units": units,
        "scale_branch": branch,
        "n_factors": int(len(variance_per_factor)),
        "dead_factors": dead,
        "active_factor_cutoff_note": (
            "Factors with max view R^2 < 1e-3 in stored units are numerically dead "
            f"(cutoff interpreted on the {units} scale)"
        ),
    }


def adata_from_mofa_factors(mofa: dict[str, Any]) -> ad.AnnData:
    """Lightweight AnnData from MOFA factors alone (no MuData required)."""
    factors = mofa["factors"]
    adata = ad.AnnData(obs=pd.DataFrame(index=factors.index))
    adata.obsm["X_mofa"] = factors.to_numpy(dtype=np.float64)
    return adata


def resolve_rna_layer(cfg: dict[str, Any], available: list[str]) -> str:
    analysis = cfg["module3"].get("analysis") or {}
    preferred = str(analysis.get("rna_layer", "spliced_unspliced_sum"))
    fallback = str(analysis.get("rna_layer_fallback", "spliced"))
    if preferred in available:
        return preferred
    if fallback in available:
        return fallback
    if available:
        return available[0]
    raise KeyError("No RNA layers found in MuData")


def focus_gene_matrix(
    cfg: dict[str, Any],
    mofa: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """
    Paired per-nucleus RNA vs ATAC values for focus genes from secondary_analysis.h5mu.

    Uses configured ``analysis.rna_layer`` (default ``spliced_unspliced_sum``).
    """
    from scipy import sparse

    focus = [g.upper() for g in (cfg["module3"]["analysis"].get("focus_genes") or [])]
    if not focus:
        return pd.DataFrame()

    paths = multiome_paths(cfg)
    h5mu = paths["secondary_h5mu"]
    if not h5mu.exists():
        raise FileNotFoundError(
            f"MuData required for focus_gene_matrix: {h5mu} "
            "(~10.4 GiB optional extension; set load_h5mu or use focus_gene_matrix_from_mofa)"
        )

    with h5py.File(h5mu, "r") as f:
        joint_bc = _decode_arr(f["obs"]["_index"][:]).astype(str)
        rna_bc = _decode_arr(f["mod/rna/obs/_index"][:]).astype(str)
        atac_bc = _decode_arr(f["mod/atac_cbg/obs/_index"][:]).astype(str)
        rna_pos = {b: i for i, b in enumerate(rna_bc)}
        atac_pos = {b: i for i, b in enumerate(atac_bc)}
        rna_order = np.array([rna_pos[b] for b in joint_bc])
        atac_order = np.array([atac_pos[b] for b in joint_bc])

        rna_hugo = _read_categorical_or_string(f["mod/rna/var"], "hugo_symbol")
        rna_ids = strip_ensembl_version(
            pd.Index(_decode_arr(f["mod/rna/var/_index"][:]).astype(str))
        ).astype(str)
        atac_syms = _decode_arr(f["mod/atac_cbg/var/_index"][:]).astype(str)

        atac_grp = f["mod/atac_cbg/X"]
        shape = tuple(int(x) for x in atac_grp.attrs["shape"])
        atac = sparse.csr_matrix(
            (
                np.asarray(atac_grp["data"][:]),
                np.asarray(atac_grp["indices"][:]),
                np.asarray(atac_grp["indptr"][:]),
            ),
            shape=shape,
        )
        available_layers = list(f["mod/rna/layers"].keys())
        layer = resolve_rna_layer(cfg, available_layers)
        rna_grp = f[f"mod/rna/layers/{layer}"]
        rna_shape = tuple(int(x) for x in rna_grp.attrs["shape"])
        enc = rna_grp.attrs.get("encoding-type", b"")
        if isinstance(enc, bytes):
            enc = enc.decode()
        rna = sparse.csr_matrix(
            (
                np.asarray(rna_grp["data"][:]),
                np.asarray(rna_grp["indices"][:]),
                np.asarray(rna_grp["indptr"][:]),
            ),
            shape=rna_shape,
        )
        layout = "CSR" if rna.indptr.shape[0] == rna_shape[0] + 1 else "unknown"

        azimuth = _read_obs_column(f["mod/rna/obs"], "azimuth_label")
        leiden = _read_obs_column(f["mod/rna/obs"], "leiden")
        azimuth_vals = (
            azimuth.to_numpy()[rna_order] if azimuth is not None else np.array([None] * len(joint_bc))
        )
        leiden_vals = (
            leiden.to_numpy()[rna_order] if leiden is not None else np.array([None] * len(joint_bc))
        )

        frames = []
        rna_hugo_up = np.char.upper(rna_hugo.astype(str))
        atac_syms_up = np.char.upper(atac_syms.astype(str))
        for g in focus:
            rna_hits = np.where(rna_hugo_up == g)[0]
            atac_hits = np.where(atac_syms_up == g)[0]
            if len(rna_hits) == 0 or len(atac_hits) == 0:
                continue
            rc = int(rna_hits[0])
            ac = int(atac_hits[0])
            rna_col = np.asarray(rna[:, rc].todense(), dtype=np.float32).ravel()[rna_order]
            atac_col = np.asarray(atac[:, ac].todense(), dtype=np.float32).ravel()[atac_order]
            frames.append(
                pd.DataFrame(
                    {
                        "barcode": joint_bc,
                        "gene_symbol": str(rna_hugo[rc]),
                        "gene_id_no_version": str(rna_ids[rc]),
                        "rna_value": rna_col,
                        "atac_value": atac_col,
                        "azimuth_label": azimuth_vals,
                        "leiden": leiden_vals,
                    }
                )
            )

    if not frames:
        out = pd.DataFrame()
    else:
        out = pd.concat(frames, ignore_index=True)
    out.attrs["rna_layer"] = layer
    out.attrs["rna_layer_encoding"] = enc
    out.attrs["rna_layer_layout"] = layout
    out.attrs["available_rna_layers"] = available_layers
    out.attrs["source"] = "secondary_analysis.h5mu"
    return out


def focus_gene_matrix_from_mofa(cfg: dict[str, Any], mofa: dict[str, Any]) -> pd.DataFrame:
    """
    Focus-gene paired values from MOFA HDF5 inputs (centered; interpretation-only path).

    Cell-type labels are joined from the committed cache in the notebook / smoke
    runner. Useful when ``load_h5mu: false``.
    """
    focus = [g.upper() for g in (cfg["module3"]["analysis"].get("focus_genes") or [])]
    if not focus:
        return pd.DataFrame()
    paths = multiome_paths(cfg)
    with h5py.File(paths["mofa_hdf5"], "r") as f:
        samples = _decode_arr(f["samples/group1"][:]).astype(str)
        rna_ids = strip_ensembl_version(
            pd.Index(_decode_arr(f["features/rna"][:]).astype(str))
        ).astype(str)
        hugo = _decode_arr(f["features_metadata/rna/hugo_symbol"][:]).astype(str)
        atac_syms = _decode_arr(f["features/atac_cbg"][:]).astype(str)
        hugo_up = np.char.upper(hugo.astype(str))
        atac_up = np.char.upper(atac_syms.astype(str))
        frames = []
        for g in focus:
            rna_hits = np.where(hugo_up == g)[0]
            atac_hits = np.where(atac_up == g)[0]
            if len(rna_hits) == 0 or len(atac_hits) == 0:
                continue
            rc, ac = int(rna_hits[0]), int(atac_hits[0])
            rna_col = np.asarray(f["data/rna/group1"][:, rc], dtype=np.float32)
            atac_col = np.asarray(f["data/atac_cbg/group1"][:, ac], dtype=np.float32)
            frames.append(
                pd.DataFrame(
                    {
                        "barcode": samples,
                        "gene_symbol": str(hugo[rc]),
                        "gene_id_no_version": str(rna_ids[rc]),
                        "rna_value": rna_col,
                        "atac_value": atac_col,
                    }
                )
            )
    if not frames:
        out = pd.DataFrame()
    else:
        out = pd.concat(frames, ignore_index=True)
    out.attrs["rna_layer"] = "mofa_input_rna"
    out.attrs["source"] = "multiome_mofa.hdf5"
    out.attrs["means_are_centered"] = True
    out.attrs["note"] = (
        "Centered MOFA inputs; not spliced counts. "
        "Do not plot as mean expression."
    )
    return out


def _read_categorical_or_string(group: h5py.Group, key: str) -> np.ndarray:
    s = _read_obs_column(group, key)
    if s is None:
        raise KeyError(key)
    return s.astype(str).to_numpy()


def _column_std(dset: h5py.Dataset, means: np.ndarray, chunk: int = 2000) -> np.ndarray:
    n, p = dset.shape
    acc = np.zeros(p, dtype=np.float64)
    for i in range(0, n, chunk):
        block = np.asarray(dset[i : i + chunk], dtype=np.float64)
        acc += ((block - means) ** 2).sum(axis=0)
    return np.sqrt(acc / float(max(n - 1, 1)))

