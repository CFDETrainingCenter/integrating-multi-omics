"""Export Module 2 integration artifacts (HuBMAP-only)."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from scripts.common.io import environment_versions, write_h5ad
from scripts.common.methods import methods_narrative_module3, methods_table_module3, module3_embedding_params
from scripts.common.paths import ensure_output_dirs, resolve, portable_path
from scripts.common.runtime import feature_space_counts, sha256_optional


def save_module2_integration_outputs(
    cfg: dict[str, Any],
    adata,
    donor_summary: pd.DataFrame,
    donor_label_crosstab: pd.DataFrame,
    composition_weighted: pd.Series,
    composition_fractions: pd.Series,
    max_across_labels: pd.Series,
    marker_matrix: pd.DataFrame,
    harmony_mixing: pd.DataFrame | None,
    harmony_silhouette: pd.DataFrame | None,
    leiden_sweep: pd.DataFrame,
    figure_paths: dict[str, Path],
    extras: dict[str, Any] | None = None,
) -> dict[str, Path]:
    ensure_output_dirs(cfg)
    tables = resolve(cfg, "outputs_tables")
    integrated_dir = resolve(cfg, "processed_integrated")
    reports = resolve(cfg, "outputs_reports")

    h5ad_path = integrated_dir / "module2_hubmap_lung_integrated.h5ad"
    write_h5ad(adata, h5ad_path)

    paths: dict[str, Path] = {
        "h5ad": h5ad_path,
        "donor_summary": tables / "module2_donor_composition_summary.tsv",
        "donor_label_crosstab": tables / "module2_donor_by_label_crosstab.tsv",
        "composition_weighted": tables / "module2_composition_weighted_pseudobulk.tsv",
        "composition_fractions": tables / "module2_celltype_fractions.tsv",
        "max_across_labels": tables / "module2_max_across_labels_pseudobulk.tsv",
        "marker_matrix": tables / "module2_hubmap_marker_means.tsv",
        "harmony_mixing": tables / "module2_harmony_mixing_metrics.tsv",
        "harmony_silhouette": tables / "module2_harmony_silhouette_by_celltype.tsv",
        "leiden_sweep": tables / "module2_leiden_resolution_sweep.tsv",
        "cell_class_membership": tables / "module2_cell_class_membership.tsv",
        "run_params": reports / "module2_integration_run_params.json",
        "methods_parameters": tables / "module2_integration_methods_parameters.tsv",
    }

    donor_summary.to_csv(paths["donor_summary"], sep="\t", index=False)
    donor_label_crosstab.to_csv(paths["donor_label_crosstab"], sep="\t", index=False)

    cw = composition_weighted.rename("composition_weighted_linear").to_frame()
    fs = feature_space_counts(composition_weighted.index)
    cw.attrs = {}  # pandas Series/Frame attrs not always written; record in params
    cw.to_csv(paths["composition_weighted"], sep="\t")
    composition_fractions.rename("fraction").to_frame().to_csv(
        paths["composition_fractions"], sep="\t"
    )
    max_across_labels.rename("max_across_labels_linear").to_frame().to_csv(
        paths["max_across_labels"], sep="\t"
    )
    marker_matrix.to_csv(paths["marker_matrix"], sep="\t", index=False)
    if harmony_mixing is None:
        pd.DataFrame(
            [{"status": "skipped", "reason": "module2.diagnostics.enabled=false"}]
        ).to_csv(paths["harmony_mixing"], sep="\t", index=False)
    else:
        harmony_mixing.to_csv(paths["harmony_mixing"], sep="\t", index=False)
    if harmony_silhouette is None:
        pd.DataFrame(
            [{"status": "skipped", "reason": "module2.diagnostics.enabled=false"}]
        ).to_csv(paths["harmony_silhouette"], sep="\t", index=False)
    else:
        harmony_silhouette.to_csv(paths["harmony_silhouette"], sep="\t", index=False)
    leiden_sweep.to_csv(paths["leiden_sweep"], sep="\t", index=False)
    membership = adata.uns.get("cell_class_membership")
    if isinstance(membership, pd.DataFrame):
        membership.to_csv(paths["cell_class_membership"], sep="\t", index=False)
    elif "azimuth_label" in adata.obs.columns:
        from scripts.nb03_integration.integrate import cell_class_membership_frame

        cell_class_membership_frame(adata.obs["azimuth_label"]).to_csv(
            paths["cell_class_membership"], sep="\t", index=False
        )
    else:
        pd.DataFrame().to_csv(paths["cell_class_membership"], sep="\t", index=False)
    methods_table_module3(cfg).to_csv(paths["methods_parameters"], sep="\t", index=False)

    input_hashes = []
    for item in cfg["module2"]["inputs"]:
        from scripts.nb02_qc.export import resolve_qc_h5ad

        bid = item.get("block_id") or item.get("primary_id")
        hpath = resolve_qc_h5ad(cfg, str(bid))
        input_hashes.append(
            {
                **item,
                "qc_h5ad": portable_path(cfg, hpath),
                "sha256": sha256_optional(hpath),
            }
        )

    extras_out = dict(extras or {})
    extras_out["feature_space"] = {
        "composition_weighted_features": (
            f"{fs['n_features_total']} "
            f"({fs['n_features_exonic']} exonic, {fs['n_features_intron']} intron)"
        ),
        **fs,
        "note": (
            "Module 2 pseudobulk is on the full integrated gene space; "
            "Module 2 DE may drop *-I intron features when exclude_intron_features=true."
        ),
    }

    params = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": input_hashes,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "random_seed": cfg["module2"].get("random_seed", 0),
        "exclude_intron_features": cfg["module2"].get("exclude_intron_features", True),
        "embedding": module3_embedding_params(cfg),
        "harmony": {
            **(cfg["module2"].get("harmony") or {}),
            "random_state": cfg["module2"].get("random_seed", 0),
        },
        "diagnostics": cfg["module2"].get("diagnostics") or {},
        "methods_narrative": methods_narrative_module3(cfg),
        "methods_table": methods_table_module3(cfg).to_dict(orient="records"),
        "environment": environment_versions(),
        "figures": {k: portable_path(cfg, v) for k, v in figure_paths.items()},
        "extras": extras_out,
        "_compute_seconds": extras_out.get("_compute_seconds"),
        "peak_rss_mb": extras_out.get("peak_rss_mb"),
        "note": (
            "Pseudobulk vectors are on the linear (expm1) scale. "
            "max_across_labels is a contrast to composition-weighted, not a bulk average. "
            "GTEx is not used in Module 2."
        ),
    }
    paths["run_params"].write_text(json.dumps(params, indent=2, default=str))
    for k, v in figure_paths.items():
        paths[f"figure_{k}"] = Path(v)
    return paths


# Backward-compatible alias (old smoke / notebook names)
def save_module3_outputs(*args, **kwargs):
    raise RuntimeError(
        "save_module3_outputs was renamed to save_module2_integration_outputs "
        "(Module 2; no GTEx artifacts)."
    )
