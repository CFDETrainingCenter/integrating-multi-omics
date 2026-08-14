"""Export Module 1 QC / preprocess artifacts (learner-facing module1_ prefix)."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from scripts.common.io import environment_versions, write_h5ad
from scripts.common.methods import methods_narrative_module2, methods_table_module2
from scripts.common.paths import ensure_output_dirs, resolve, portable_path
from scripts.common.runtime import sha256_optional


def module1_h5ad_path(cfg: dict[str, Any], block_id: str) -> Path:
    return resolve(cfg, "processed_hubmap") / f"module1_{block_id}_qc_preprocessed.h5ad"


# Backward-compatible alias used by older call sites / smoke runner
def module2_h5ad_path(cfg: dict[str, Any], block_id: str) -> Path:
    return module1_h5ad_path(cfg, block_id)


def resolve_qc_h5ad(cfg: dict[str, Any], block_id: str) -> Path:
    """Prefer module1_ outputs; fall back to legacy module2_ filenames if present."""
    new = module1_h5ad_path(cfg, block_id)
    if new.exists():
        return new
    legacy = resolve(cfg, "processed_hubmap") / f"module2_{block_id}_qc_preprocessed.h5ad"
    return legacy if legacy.exists() else new


def save_module1_outputs(
    cfg: dict[str, Any],
    adata,
    filter_log: pd.DataFrame,
    metrics_summary: pd.DataFrame,
    figure_paths: dict[str, Path],
    run_extras: dict[str, Any] | None = None,
) -> dict[str, Path]:
    """Write per-block QC tables/figures so multi-donor runs do not overwrite (M7)."""
    ensure_output_dirs(cfg)
    block_id = str(adata.uns.get("hubmap_block_id", cfg["module1"]["block_id"]))
    donor_label = str(
        adata.uns.get("hubmap_donor_label")
        or (adata.obs["donor_label"].iloc[0] if "donor_label" in adata.obs.columns else "")
        or cfg["module1"].get("donor_label", "")
    )

    h5ad_path = module1_h5ad_path(cfg, block_id)
    write_h5ad(adata, h5ad_path)

    tables = resolve(cfg, "outputs_tables")
    # Namespace by block_id (M7). Also keep a copy of filter log columns with donor_label.
    filter_log = filter_log.copy()
    if "donor_label" not in filter_log.columns:
        filter_log.insert(0, "donor_label", donor_label)
    if "block_id" not in filter_log.columns:
        filter_log.insert(1, "block_id", block_id)

    metrics_summary = metrics_summary.copy()
    if "donor_label" not in metrics_summary.columns:
        metrics_summary.insert(0, "donor_label", donor_label)
    if "block_id" not in metrics_summary.columns:
        metrics_summary.insert(1, "block_id", block_id)

    filter_path = tables / f"module1_qc_filter_log_{block_id}.tsv"
    metrics_path = tables / f"module1_qc_metrics_summary_{block_id}.tsv"
    filter_log.to_csv(filter_path, sep="\t", index=False)
    metrics_summary.to_csv(metrics_path, sep="\t", index=False)

    params = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "module": 1,
        "block_id": block_id,
        "donor_label": donor_label,
        "donor_id": adata.obs["donor_id"].iloc[0] if "donor_id" in adata.obs else None,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "qc": cfg["module1"]["qc"],
        "preprocess": cfg["module1"]["preprocess"],
        "doublet_filtering": bool(cfg["module1"].get("doublet_filtering", False)),
        "max_nuclei": cfg["module1"].get("max_nuclei"),
        "random_seed": cfg["module1"].get("random_seed"),
        "source_file": portable_path(cfg, adata.uns.get("hubmap_source_file") or ""),
        "inputs": {
            "raw_expr": {
                "path": portable_path(cfg, adata.uns.get("hubmap_source_file") or ""),
                "sha256": sha256_optional(adata.uns.get("hubmap_source_file")),
            }
        },
        "methods_narrative": methods_narrative_module2(cfg),
        "methods_table": methods_table_module2(cfg).to_dict(orient="records"),
        "environment": environment_versions(),
        "figures": {k: portable_path(cfg, v) for k, v in figure_paths.items()},
        "extras": run_extras or {},
        "_compute_seconds": (run_extras or {}).get("_compute_seconds"),
        "peak_rss_mb": (run_extras or {}).get("peak_rss_mb"),
    }
    params_path = resolve(cfg, "outputs_reports") / f"module1_run_params_{block_id}.json"
    params_path.write_text(json.dumps(params, indent=2))

    methods_path = tables / "module1_methods_parameters.tsv"
    methods_table_module2(cfg).to_csv(methods_path, sep="\t", index=False)

    # Refresh combined multi-donor QC log after each save
    combined = combine_qc_filter_logs(cfg)
    coverage_path = upsert_azimuth_join_coverage(
        cfg,
        donor_label=donor_label,
        block_id=block_id,
        join_info=(run_extras or {}).get("label_join") or {},
        n_obs_saved=int(adata.n_obs),
    )

    return {
        "h5ad": h5ad_path,
        "filter_log": filter_path,
        "metrics_summary": metrics_path,
        "methods_parameters": methods_path,
        "run_params": params_path,
        "filter_log_all_donors": combined,
        "azimuth_join_coverage": coverage_path,
        **{f"figure_{k}": Path(v) for k, v in figure_paths.items()},
    }


def upsert_azimuth_join_coverage(
    cfg: dict[str, Any],
    *,
    donor_label: str,
    block_id: str,
    join_info: dict[str, Any],
    n_obs_saved: int,
) -> Path:
    """Append/replace one donor row in module1_azimuth_join_coverage.tsv."""
    tables = resolve(cfg, "outputs_tables")
    path = tables / "module1_azimuth_join_coverage.tsv"
    cols = [
        "donor_label",
        "block_id",
        "n_raw",
        "n_secondary",
        "n_overlap",
        "n_raw_without_label",
        "fraction_raw_with_label",
        "n_obs_saved",
        "joined",
        "note",
    ]
    row = {
        "donor_label": donor_label,
        "block_id": block_id,
        "n_raw": join_info.get("n_raw"),
        "n_secondary": join_info.get("n_secondary"),
        "n_overlap": join_info.get("n_overlap"),
        "n_raw_without_label": join_info.get("n_raw_without_label"),
        "fraction_raw_with_label": join_info.get("fraction_raw_with_label"),
        "n_obs_saved": int(n_obs_saved),
        "joined": join_info.get("joined"),
        "note": "live from save_module1_outputs / label_join",
    }
    # Back-fill coverage fields when older run_params only stored overlap counts.
    try:
        n_raw = int(row["n_raw"]) if row["n_raw"] is not None and str(row["n_raw"]) != "nan" else None
        n_overlap = (
            int(row["n_overlap"])
            if row["n_overlap"] is not None and str(row["n_overlap"]) != "nan"
            else None
        )
    except (TypeError, ValueError):
        n_raw, n_overlap = None, None
    if n_raw is not None and n_overlap is not None:
        if row["n_raw_without_label"] is None or str(row["n_raw_without_label"]) in {"", "nan", "<NA>"}:
            row["n_raw_without_label"] = int(n_raw - n_overlap)
        if row["fraction_raw_with_label"] is None or str(row["fraction_raw_with_label"]) in {
            "",
            "nan",
            "<NA>",
        }:
            row["fraction_raw_with_label"] = float(n_overlap / n_raw) if n_raw else float("nan")

    if path.exists():
        df = pd.read_csv(path, sep="\t")
        for c in cols:
            if c not in df.columns:
                df[c] = pd.NA
        df = df[df["donor_label"].astype(str) != str(donor_label)]
        df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
    else:
        df = pd.DataFrame([row], columns=cols)
    # Stable donor order if recognizable
    order = {"Donor_1": 0, "Donor_2": 1, "Donor_3": 2, "Donor_4": 3}
    df["_ord"] = df["donor_label"].map(lambda x: order.get(str(x), 99))
    df = df.sort_values(["_ord", "donor_label"]).drop(columns=["_ord"])
    df.to_csv(path, sep="\t", index=False)
    return path


# Alias for smoke / notebooks still calling the old name
def save_module2_outputs(*args, **kwargs) -> dict[str, Path]:
    return save_module1_outputs(*args, **kwargs)


def combine_qc_filter_logs(cfg: dict[str, Any]) -> Path:
    """Concatenate per-block module1_qc_filter_log_*.tsv -> all_donors table."""
    tables = resolve(cfg, "outputs_tables")
    paths = sorted(tables.glob("module1_qc_filter_log_*.tsv"))
    paths = [p for p in paths if p.name != "module1_qc_filter_log_all_donors.tsv"]
    out = tables / "module1_qc_filter_log_all_donors.tsv"
    if not paths:
        pd.DataFrame(columns=["donor_label", "block_id", "step", "n_nuclei", "n_genes", "rationale"]).to_csv(
            out, sep="\t", index=False
        )
        return out
    frames = [pd.read_csv(p, sep="\t") for p in paths]
    pd.concat(frames, ignore_index=True).to_csv(out, sep="\t", index=False)
    return out
