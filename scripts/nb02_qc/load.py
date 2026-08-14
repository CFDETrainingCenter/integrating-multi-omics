"""Load Donor_1 snRNA-seq AnnData and attach metadata / optional labels."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from scripts.common.io import read_h5ad
from scripts.common.paths import snrna_block_dir
from scripts.common.provenance import parse_block_metadata


def block_paths(
    cfg: dict[str, Any],
    block_id: str | None = None,
    donor_label: str | None = None,
) -> dict[str, Path]:
    block = block_id or cfg["module1"]["block_id"]
    label = donor_label or cfg["module1"].get("donor_label", "Donor_1")
    d = snrna_block_dir(cfg, block, donor_label=label)
    tsvs = list(d.glob("*.tsv"))
    return {
        "block_dir": d,
        "block_id": block,
        "donor_label": label,
        "raw_expr": d / "raw_expr.h5ad",
        "expr": d / "expr.h5ad",
        "secondary_analysis": d / "secondary_analysis.h5ad",
        "tsv": tsvs[0] if tsvs else d / f"{block}.tsv",
    }


def load_raw_expr(
    cfg: dict[str, Any],
    block_id: str | None = None,
    donor_label: str | None = None,
):
    paths = block_paths(cfg, block_id, donor_label=donor_label)
    if not paths["raw_expr"].exists():
        raise FileNotFoundError(f"Missing raw_expr.h5ad at {paths['raw_expr']}")
    adata = read_h5ad(paths["raw_expr"])
    adata.uns["hubmap_block_id"] = paths["block_id"]
    adata.uns["hubmap_donor_label"] = paths["donor_label"]
    adata.uns["hubmap_source_file"] = str(paths["raw_expr"])
    return adata, paths


def _donor_cfg(cfg: dict[str, Any], donor_label: str) -> dict[str, Any]:
    label = str(donor_label)
    if label.lower().startswith("donor_"):
        key = f"donor_{label.split('_', 1)[1]}"
    else:
        key = "donor_1"
    return cfg.get(key, cfg["donor_1"])


def attach_donor_metadata(
    adata,
    cfg: dict[str, Any],
    tsv_path: Path,
    donor_label: str | None = None,
):
    meta = parse_block_metadata(tsv_path)
    label = donor_label or adata.uns.get("hubmap_donor_label", "Donor_1")
    donor = _donor_cfg(cfg, label)
    adata.obs["donor_id"] = meta.get("donor_id") or donor["donor_id"]
    adata.obs["donor_label"] = donor.get("label", label)
    adata.obs["sample_id"] = meta.get("primary_dataset_id") or adata.uns.get("hubmap_block_id")
    adata.obs["block_id"] = adata.uns.get("hubmap_block_id")
    adata.obs["sex"] = meta.get("sex") or donor.get("sex")
    adata.obs["age_value"] = meta.get("age_value") or str(donor.get("age", ""))
    adata.obs["race"] = meta.get("race") or donor.get("race")
    adata.obs["organ"] = donor.get("organ", "Lung")
    adata.obs["assay"] = meta.get("rnaseq_assay_method") or "snRNAseq (10x Genomics v3)"
    adata.uns["hubmap_donor_metadata"] = meta
    return adata


def _strip_ensembl_version(ids: pd.Index) -> pd.Index:
    return ids.str.replace(r"\.\d+$", "", regex=True)


def annotate_mito_ribo(adata):
    """Mark MT / ribo genes using var index or hugo_symbol if present."""
    symbols = None
    if "hugo_symbol" in adata.var.columns:
        symbols = adata.var["hugo_symbol"].astype(str)
    else:
        # fall back to gene id / index string match
        symbols = pd.Series(adata.var_names.astype(str), index=adata.var_names)

    adata.var["mt"] = symbols.str.upper().str.startswith("MT-")
    adata.var["ribo"] = symbols.str.upper().str.match(r"^(RPS|RPL)\d")
    return adata


def join_azimuth_labels_from_cache(adata, cfg: dict[str, Any], donor_label: str):
    """Join Azimuth labels from the committed cache TSV.

    The minimum learner download deliberately excludes secondary_analysis.h5ad
    (5.67 GB across the four teaching donors), so the cache is the only label
    source on that path. Long format written by
    scripts.common.remote_obs.write_azimuth_label_cache.
    """
    import pandas as pd

    from scripts.common.paths import module_root

    rel = ((cfg.get("module2") or {}).get("labels") or {}).get("cache")
    if not rel:
        return adata, {"joined": False, "reason": "module2.labels.cache not configured"}
    path = module_root(cfg) / rel
    if not path.exists():
        return adata, {"joined": False, "reason": f"label cache missing: {rel}"}

    df = pd.read_csv(path, sep="\t")
    sub = df[df["donor_label"].astype(str) == str(donor_label)]
    if sub.empty:
        return adata, {
            "joined": False,
            "reason": f"label cache has no rows for {donor_label}",
            "path": rel,
        }
    mapping = dict(zip(sub["barcode"].astype(str), sub["azimuth_label"].astype(str)))
    labels = adata.obs_names.astype(str).map(mapping)
    n_overlap = int(labels.notna().sum())
    adata.obs["azimuth_label"] = labels.fillna("unlabeled").astype(str).values
    return adata, {
        "joined": True,
        "source": "cache",
        "path": rel,
        "n_raw": int(adata.n_obs),
        "n_overlap": n_overlap,
        "n_raw_without_label": int(adata.n_obs) - n_overlap,
        "fraction_raw_with_label": round(n_overlap / adata.n_obs, 6) if adata.n_obs else 0.0,
    }


def join_azimuth_labels(adata, secondary_path: Path, label_cols: list[str] | None = None):
    """Join annotation columns from secondary_analysis onto raw barcodes."""
    if not secondary_path.exists():
        return adata, {"joined": False, "reason": "secondary_analysis missing"}

    label_cols = label_cols or [
        "azimuth_label",
        "predicted_label",
        "leiden",
    ]
    try:
        sec = read_h5ad(secondary_path)
    except OSError as exc:
        return adata, {
            "joined": False,
            "reason": f"secondary_analysis unreadable/truncated: {exc}",
            "path": str(secondary_path),
        }

    available = [c for c in label_cols if c in sec.obs.columns]
    if not available:
        return adata, {"joined": False, "reason": "no requested label columns present"}

    overlap = adata.obs_names.intersection(sec.obs_names)
    n_raw = int(adata.n_obs)
    n_overlap = int(len(overlap))
    n_raw_without_label = n_raw - n_overlap
    info = {
        "joined": True,
        "n_raw": n_raw,
        "n_secondary": int(sec.n_obs),
        "n_overlap": n_overlap,
        "n_raw_without_label": n_raw_without_label,
        "fraction_raw_with_label": round(n_overlap / n_raw, 6) if n_raw else None,
        "columns": available,
        "note": (
            "Labels live on secondary_analysis barcodes; raw_expr may carry additional "
            "barcodes that remain unlabeled until QC drops them."
        ),
    }
    for col in available:
        joined = pd.Series(pd.NA, index=adata.obs_names, dtype="string")
        joined.loc[overlap] = sec.obs.loc[overlap, col].astype(str).values
        adata.obs[col] = joined
        if col == "azimuth_label":
            n_labeled = int(joined.notna().sum())
            n_missing = int(joined.isna().sum())
            info["n_azimuth_labeled"] = n_labeled
            info["n_azimuth_missing"] = n_missing
    return adata, info


def load_module2_adata(
    cfg: dict[str, Any],
    block_id: str | None = None,
    donor_label: str | None = None,
):
    """Full Module 2 load: raw_expr + donor metadata + optional Azimuth labels."""
    label = donor_label or cfg["module1"].get("donor_label", "Donor_1")
    adata, paths = load_raw_expr(cfg, block_id, donor_label=label)
    adata = attach_donor_metadata(adata, cfg, paths["tsv"], donor_label=label)
    adata = annotate_mito_ribo(adata)

    join_info = {"joined": False, "reason": "disabled in config"}
    if cfg["module1"].get("join_azimuth_labels", True):
        adata, join_info = join_azimuth_labels(adata, paths["secondary_analysis"])
        if not join_info.get("joined"):
            # Minimum-download path: secondary_analysis.h5ad is not fetched, so fall
            # back to the committed cache before giving up.
            adata, cache_info = join_azimuth_labels_from_cache(adata, cfg, label)
            if cache_info.get("joined"):
                cache_info["local_join_reason"] = join_info.get("reason")
                join_info = cache_info
            else:
                # Never continue unlabelled. cell_class would silently become
                # "unlabeled" for every nucleus and Modules 2 and 4 would report
                # plausible but wrong composition.
                raise FileNotFoundError(
                    "Azimuth labels unavailable for "
                    f"{label}: local secondary_analysis ({join_info.get('reason')}) "
                    f"and cache ({cache_info.get('reason')}). Build the cache with "
                    "python -c \"from scripts.common.paths import load_config; "
                    "from scripts.common.remote_obs import write_azimuth_label_cache; "
                    "print(write_azimuth_label_cache(load_config()))\" "
                    "or set module1.join_azimuth_labels: false to run without labels."
                )

    adata.var["gene_id_no_version"] = _strip_ensembl_version(adata.var_names.to_series()).values
    return adata, paths, join_info
