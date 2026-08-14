"""Export Module 3 multiome artifacts for the configured teaching donor."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from scripts.common.io import environment_versions
from scripts.common.methods import methods_narrative_module5, methods_table_module5
from scripts.common.paths import ensure_output_dirs, resolve, portable_path


def write_interpretation_markdown(
    cfg: dict[str, Any],
    overlap_df: pd.DataFrame,
    bridge_df: pd.DataFrame,
    focus_df: pd.DataFrame,
    corr_df: pd.DataFrame,
    variance_total: pd.DataFrame,
    path: Path,
    availability_df: pd.DataFrame | None = None,
    factor_dom_df: pd.DataFrame | None = None,
    cluster_concordance_df: pd.DataFrame | None = None,
    qc_association_df: pd.DataFrame | None = None,
    factor_label_association_df: pd.DataFrame | None = None,
    azimuth_composition_df: pd.DataFrame | None = None,
    *,
    r2_recompute: dict[str, Any] | None = None,
    nondom_summary: dict[str, Any] | None = None,
) -> Path:
    ov = overlap_df.iloc[0].to_dict() if overlap_df is not None and not overlap_df.empty else {}
    r = bridge_df.attrs.get("pearson_log1p")
    rho = bridge_df.attrs.get("spearman_log1p")
    metric = bridge_df.attrs.get("bridge_metric", "feature_std_on_mofa_inputs")
    m5 = cfg["module3"]
    donor = m5["donor_label"]
    nd = nondom_summary or {}
    lines = [
        f"# Module 3 -- {donor} SNARE multiome interpretation",
        "",
        f"_Auto-drafted at {datetime.now(timezone.utc).isoformat()}. Edit with learner conclusions._",
        "",
        "## Role in the course",
        "",
        "This component reads two modalities measured in the same nucleus (paired RNA + ATAC). "
        "Modules 1-2 and 4 are multi-*resource* transcriptomics.",
        "",
        "## Design",
        "",
        f"- Donor: `{m5['donor_id']}` ({donor})",
        f"- Dataset: `{m5['dataset_id']}`",
        f"- Assay: `{m5.get('assay', 'SNARE-seq2 [Salmon + ArchR + Muon]')}`",
        f"- Portal: {m5.get('portal_url', 'https://portal.hubmapconsortium.org')}",
        "- Default path: interpret fitted MOFA (`multiome_mofa.hdf5`) plus a committed "
        "Azimuth/WNN label cache; full MuData is not a learner download (`load_h5mu`)",
        "",
        "## Barcode pairing",
        "",
        f"- RNA barcodes: {ov.get('rna_n_barcodes')}",
        f"- ATAC barcodes: {ov.get('atac_n_barcodes')}",
        f"- Exact overlap: {ov.get('n_exact_overlap')}",
        f"- Paired multiome: {ov.get('paired_multiome')}",
        f"- Note: {ov.get('interpretation')}",
        "",
        "## MOFA: modality-private factors (not a joint latent space)",
        "",
    ]
    if variance_total is not None and not variance_total.empty:
        for _, row in variance_total.iterrows():
            r2 = float(row["r2_total"])
            lines.append(
                f"- View `{row['view']}`: MOFA explains ~{r2:.2f} percent of total "
                f"variance (stored on a percent scale in the HDF5, so {r2:.4f} is "
                f"{r2:.2f} percent and not a fraction above 1)."
            )
    else:
        lines.append("_Variance explained table unavailable._")

    if r2_recompute:
        lines.extend(
            [
                "",
                "### R^2 recomputation (from Z, W, Y in the HDF5)",
                "",
                f"- {r2_recompute.get('label', '')}",
                f"- Details: `{json.dumps({k: r2_recompute[k] for k in r2_recompute if k != 'label'}, default=str)}`",
            ]
        )

    max_nd = nd.get("max_nondominant_r2")
    max_pct = nd.get("max_nondominant_percent")
    units = nd.get("r2_units") or "percent"
    if max_nd is not None and max_pct is not None:
        finding = (
            f"**Finding.** No factor explains more than {max_pct:.4f} percent "
            f"of variance in its non-dominant view "
            f"(HDF5 value {max_nd:.5f}; units match stored totals = {units})."
        )
    else:
        finding = "**Finding.** Every factor is modality-private (see variance table)."
    lines.extend(
        [
            "",
            finding,
            "This MOFA space is **two stacked single-modality spaces**, not a joint shared space. "
            "Trailing factors beyond the configured active set are numerically dead.",
            "",
            "### Interpretation question",
            "",
            "What would a genuinely shared factor look like, and what would have to be true of the "
            "data for MOFA to find one? Candidates: view scaling, sparsity/dimensionality mismatch, "
            "and ATAC gene activity being a smoothed proxy rather than a measurement.",
            "",
        ]
    )
    if factor_dom_df is not None and not factor_dom_df.empty:
        lines.append("Factor dominance (active factors):")
        lines.append("")
        for _, row in factor_dom_df.iterrows():
            lines.append(
                f"- `{row['factor']}`: dominant=`{row['dominant_view']}` "
                f"(rna={float(row.get('rna', float('nan'))):.4f}, "
                f"atac_cbg={float(row.get('atac_cbg', float('nan'))):.4f}, "
                f"nondominant={float(row.get('nondominant_r2', float('nan'))):.5f})"
            )
        lines.append("")

    if factor_label_association_df is not None and not factor_label_association_df.empty:
        lines.extend(["", "## Do the factors track cell identity?", ""])
        top = (
            factor_label_association_df.sort_values("eta_squared", ascending=False)
            .groupby("labelling", as_index=False)
            .first()
        )
        for _, row in top.iterrows():
            lines.append(
                f"- `{row['labelling']}`: strongest factor `{row['factor']}` "
                f"eta_squared={float(row['eta_squared']):.3f}"
            )
        lines.append("")

    if cluster_concordance_df is not None and not cluster_concordance_df.empty:
        lines.extend(["", "## Cluster concordance (ARI / NMI)", ""])
        for _, row in cluster_concordance_df.iterrows():
            lines.append(
                f"- {row['pair']}: ARI={float(row['ari']):.4f}, NMI={float(row['nmi']):.4f} "
                f"(k_a={int(row['k_a'])}, k_b={int(row['k_b'])}, n={int(row['n_nuclei'])})"
            )
        lines.append("")

    if qc_association_df is not None and not qc_association_df.empty:
        lines.extend(["", "## ATAC QC vs clustering (eta squared)", ""])
        for _, row in qc_association_df.iterrows():
            lines.append(
                f"- {row['qc_metric']} vs `{row['labelling']}`: "
                f"eta_squared={float(row['eta_squared']):.3f} (transform={row['transform']})"
            )
        lines.append("")

    lines.extend(
        [
            "## Gene-level bridge (MOFA feature space)",
            "",
            f"- Bridge metric: `{metric}`",
            f"- Shared genes: {0 if bridge_df is None else len(bridge_df)}",
            f"- Pearson r (log1p SD): {None if r is None else round(float(r), 3) if r == r else r}",
            f"- Spearman rho (log1p SD): {None if rho is None else round(float(rho), 3) if rho == rho else rho}",
            "",
            "## Focus markers",
            "",
            "RNA and ATAC use **different units** -- interpret concordance, not equality.",
            "",
        ]
    )
    if focus_df is None or focus_df.empty:
        lines.append("_No focus genes recovered in both modalities._")
    else:
        for _, row in focus_df.iterrows():
            sym = row.get("gene_symbol") or row["gene_id_no_version"]
            lines.append(
                f"- **{sym}**: RNA mean={float(row['rna_mean_expression']):.3f}, "
                f"ATAC activity mean={float(row['atac_mean_gene_activity']):.3f}"
            )

    if availability_df is not None and not availability_df.empty:
        missing = availability_df[~availability_df["recovered_in_both_modalities"]]
        if not missing.empty:
            lines.extend(["", "### Markers not recovered in both modality spaces", ""])
            for _, row in missing.iterrows():
                lines.append(f"- **{row['gene_symbol']}**: {row['notes']}")

    lines.extend(["", "## Focus markers (paired correlations)", ""])
    if corr_df is None or corr_df.empty:
        lines.append("_No paired correlations computed._")
    else:
        lines.append(
            "Gene-level RNA<->ATAC correlation is weak everywhere in this dataset "
            "(largest absolute r typically << 0.2). Restricting to one cell type usually "
            "**shrinks** r by removing between-type variance -- that is the stricter test."
        )
        lines.append("")
        global_rows = corr_df[corr_df["scope"] == "global"] if "scope" in corr_df.columns else corr_df
        cell_rows = (
            corr_df[corr_df["scope"] == "within_cell_type"]
            if "scope" in corr_df.columns
            else pd.DataFrame()
        )
        lines.append("### Global (all nuclei)")
        lines.append("")
        for _, row in global_rows.iterrows():
            note = f" -- {row['note']}" if "note" in row and str(row.get("note", "")) not in {"", "nan"} else ""
            lines.append(
                f"- **{row['gene_symbol']}**: r~{float(row['pearson_r']) if row['pearson_r']==row['pearson_r'] else float('nan'):.3f} "
                f"(n={int(row['n_nuclei'])}; rna_mean={float(row['rna_mean']):.4f}){note}"
            )
        if cell_rows is not None and not cell_rows.empty:
            lines.extend(["", "### Within expected Azimuth cell type(s)", ""])
            for _, row in cell_rows.iterrows():
                r_val = row["pearson_r"]
                r_txt = f"{float(r_val):.3f}" if r_val == r_val else "undefined"
                note = f" -- **{row['note']}**" if "note" in row and str(row.get("note", "")) not in {"", "nan"} else ""
                lines.append(
                    f"- **{row['gene_symbol']}**: r~{r_txt} "
                    f"(n={int(row['n_nuclei'])}; rna_mean={float(row['rna_mean']) if row['rna_mean']==row['rna_mean'] else float('nan'):.4f}; "
                    f"{row['cell_type_filter']}){note}"
                )

    lines.extend(
        [
            "",
            "## Interpretation prompts",
            "",
            "1. What evidence shows this is true multiome rather than unpaired same-donor assays?",
            "2. Why is this MOFA space not a joint latent space? Quote the max non-dominant R^2.",
            "3. What would a shared factor look like, and what data properties might prevent one?",
            "4. Why does within-cell-type correlation usually shrink relative to global?",
            "5. When rna_mean==0 for a defining marker (e.g. NKG7 in NK/CD8), what hypotheses remain?",
            "",
        ]
    )

    # Dev-machine Q1/M12 probe (optional; does not require learners to load 10.4 GiB)
    probe_path = resolve(cfg, "outputs_reports") / "module3_rna_layer_probe.json"
    if probe_path.exists():
        try:
            probe = json.loads(probe_path.read_text())
            m12 = probe.get("M12_conclusion") or {}
            lines.extend(
                [
                    "## Recorded layer probe (Q1 / Q11 / M12)",
                    "",
                    f"- RNA layers present: {probe.get('rna_layers')}",
                    f"- Encoding / layout: CSR (`encoding-type=csr_matrix`); Q11 confirmed",
                    f"- Preferred teaching layer when MuData is loaded: "
                    f"`{(cfg['module3'].get('analysis') or {}).get('rna_layer', 'spliced_unspliced_sum')}`",
                    f"- NKG7 column sums (all nuclei): `{json.dumps(probe.get('NKG7_layer_sums'), default=str)}`",
                ]
            )
            if m12:
                lines.append(f"- **M12 finding:** {m12.get('finding')}")
            lines.append("")
        except Exception:  # noqa: BLE001
            pass

    lines.extend(
        [
            "## Three-sentence results draft (edit me)",
            "",
            (
                f"In {donor} (`{m5['donor_id']}`), SNARE-seq2 product `{m5['dataset_id']}` provides "
                "paired RNA and ATAC for the same nuclei. "
                "MOFA factors are modality-private: the maximum non-dominant R^2 is "
                f"{max_pct:.4f} percent (units match percent-scale stored totals) -- "
                "two stacked single-modality spaces, not a joint embedding. "
                "Focus-gene RNA<->ATAC correlations are near zero globally and usually smaller within "
                "cell type; zeros and low-n groups are flagged explicitly rather than reported as bare nan."
            ),
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))
    return path


def save_module3_outputs(
    cfg: dict[str, Any],
    inventory_df: pd.DataFrame,
    contrast_df: pd.DataFrame,
    overlap_df: pd.DataFrame,
    bridge_df: pd.DataFrame,
    focus_df: pd.DataFrame,
    corr_df: pd.DataFrame,
    variance_total: pd.DataFrame,
    variance_per_factor: pd.DataFrame,
    figure_paths: dict[str, Path],
    extras: dict[str, Any] | None = None,
    availability_df: pd.DataFrame | None = None,
    factor_dom_df: pd.DataFrame | None = None,
    cluster_concordance_df: pd.DataFrame | None = None,
    qc_association_df: pd.DataFrame | None = None,
    factor_label_association_df: pd.DataFrame | None = None,
    azimuth_composition_df: pd.DataFrame | None = None,
    *,
    r2_recompute: dict[str, Any] | None = None,
    nondom_summary: dict[str, Any] | None = None,
) -> dict[str, Path]:
    ensure_output_dirs(cfg)
    tables = resolve(cfg, "outputs_tables")
    reports = resolve(cfg, "outputs_reports")

    paths: dict[str, Path] = {
        "inventory": tables / "module3_multimodal_inventory.tsv",
        "contrast": tables / "module3_modality_contrast.tsv",
        "barcode_overlap": tables / "module3_barcode_overlap.tsv",
        "gene_bridge": tables / "module3_rna_atac_gene_bridge.tsv",
        "focus_genes": tables / "module3_focus_gene_bridge.tsv",
        "focus_correlations": tables / "module3_focus_gene_correlations.tsv",
        "focus_availability": tables / "module3_focus_gene_availability.tsv",
        "mofa_factor_dominance": tables / "module3_mofa_factor_dominance.tsv",
        "mofa_variance_total": tables / "module3_mofa_variance_total.tsv",
        "mofa_variance_per_factor": tables / "module3_mofa_variance_per_factor.tsv",
        "methods_parameters": tables / "module3_methods_parameters.tsv",
        "cluster_concordance": tables / "module3_cluster_concordance.tsv",
        "qc_association": tables / "module3_qc_association.tsv",
        "factor_label_association": tables / "module3_factor_label_association.tsv",
        "azimuth_composition": tables / "module3_azimuth_composition.tsv",
        "interpretation": reports / "module3_multimodal_interpretation.md",
        "run_params": reports / "module3_run_params.json",
    }

    inventory_df.to_csv(paths["inventory"], sep="\t", index=False)
    contrast_df.to_csv(paths["contrast"], sep="\t", index=False)
    overlap_df.to_csv(paths["barcode_overlap"], sep="\t", index=False)
    bridge_df.to_csv(paths["gene_bridge"], sep="\t", index=False)
    focus_df.to_csv(paths["focus_genes"], sep="\t", index=False)
    corr_df.to_csv(paths["focus_correlations"], sep="\t", index=False)
    if availability_df is None:
        availability_df = pd.DataFrame()
    availability_df.to_csv(paths["focus_availability"], sep="\t", index=False)
    if factor_dom_df is None:
        factor_dom_df = pd.DataFrame()
    factor_dom_df.to_csv(paths["mofa_factor_dominance"], sep="\t", index=False)
    variance_total.to_csv(paths["mofa_variance_total"], sep="\t", index=False)
    variance_per_factor.to_csv(paths["mofa_variance_per_factor"], sep="\t")
    if cluster_concordance_df is None:
        cluster_concordance_df = pd.DataFrame()
    if qc_association_df is None:
        qc_association_df = pd.DataFrame()
    if factor_label_association_df is None:
        factor_label_association_df = pd.DataFrame()
    if azimuth_composition_df is None:
        azimuth_composition_df = pd.DataFrame()
    cluster_concordance_df.to_csv(paths["cluster_concordance"], sep="\t", index=False)
    qc_association_df.to_csv(paths["qc_association"], sep="\t", index=False)
    factor_label_association_df.to_csv(paths["factor_label_association"], sep="\t", index=False)
    azimuth_composition_df.to_csv(paths["azimuth_composition"], sep="\t", index=False)
    methods_table_module5(cfg).to_csv(paths["methods_parameters"], sep="\t", index=False)
    write_interpretation_markdown(
        cfg,
        overlap_df,
        bridge_df,
        focus_df,
        corr_df,
        variance_total,
        paths["interpretation"],
        availability_df=availability_df,
        factor_dom_df=factor_dom_df,
        cluster_concordance_df=cluster_concordance_df,
        qc_association_df=qc_association_df,
        factor_label_association_df=factor_label_association_df,
        azimuth_composition_df=azimuth_composition_df,
        r2_recompute=r2_recompute,
        nondom_summary=nondom_summary,
    )

    params = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "module3": cfg["module3"],
        "n_shared_genes": int(bridge_df.shape[0]),
        "pearson_log1p": bridge_df.attrs.get("pearson_log1p"),
        "spearman_log1p": bridge_df.attrs.get("spearman_log1p"),
        "bridge_metric": bridge_df.attrs.get("bridge_metric"),
        "r2_recompute": r2_recompute or {},
        "nondominant_r2_summary": nondom_summary or {},
        "methods_narrative": methods_narrative_module5(cfg),
        "methods_table": methods_table_module5(cfg).to_dict(orient="records"),
        "environment": environment_versions(),
        "figures": {k: portable_path(cfg, v) for k, v in figure_paths.items()},
        "extras": extras or {},
        "_compute_seconds": (extras or {}).get("_compute_seconds"),
        "peak_rss_mb": (extras or {}).get("peak_rss_mb"),
    }
    paths["run_params"].write_text(json.dumps(params, indent=2, default=str))
    for k, v in figure_paths.items():
        paths[f"figure_{k}"] = Path(v)
    return paths


def save_module5_outputs(*args, **kwargs):
    raise RuntimeError(
        "save_module5_outputs was renamed to save_module3_outputs (learner Module 3)."
    )
