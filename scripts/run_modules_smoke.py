#!/usr/bin/env python
"""Non-interactive smoke runs for Modules 1-5."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

MODULE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_ROOT))


def run_module1() -> None:
    from scripts.common.paths import load_config
    from scripts.nb01_discovery.access_table import write_access_provenance_table
    from scripts.nb01_discovery.composition_gate import write_composition_gate
    from scripts.nb01_discovery.inventory import write_inventory_outputs
    from scripts.nb01_discovery.verify import write_verification_outputs

    cfg = load_config()
    inv = write_inventory_outputs(cfg)
    access = write_access_provenance_table(cfg)
    ver = write_verification_outputs(cfg)
    gate = write_composition_gate(cfg)
    print("Module 1 discovery/verify/composition-gate OK")
    for k, p in {**inv, **ver, **gate}.items():
        print(f"  {k}: {p} exists={p.exists()}")
    print(f"  access: {access} exists={access.exists()}")


def run_module2(block_id: str | None = None, donor_label: str | None = None) -> Path:
    """QC + preprocess for one teaching block (learner Module 1 analysis half)."""
    from scripts.common.paths import load_config, resolve
    from scripts.common.plotting import apply_figure_style, qc_scatter, qc_violin_panel, umap_panel
    from scripts.common.runtime import finalize_timing, measure_run
    from scripts.nb02_qc.export import save_module1_outputs
    from scripts.nb02_qc.load import load_module2_adata
    from scripts.nb02_qc.preprocess import (
        append_subsample_filter_row,
        hvg_pca,
        maybe_subsample,
        neighbors_umap,
        normalize_log,
    )
    from scripts.nb02_qc.qc import apply_filters, compute_qc_metrics, qc_metrics_summary

    with measure_run() as timing:
        cfg = load_config()
        apply_figure_style()
        block = block_id or cfg["module1"]["block_id"]
        donor = donor_label or cfg["module1"].get("donor_label", "Donor_1")
        adata, paths, join_info = load_module2_adata(cfg, block_id=block, donor_label=donor)
        print("Loaded", paths["raw_expr"], adata.shape, "join", join_info)
        adata.uns["hubmap_donor_label"] = donor

        adata = compute_qc_metrics(adata)
        fig_dir = resolve(cfg, "outputs_figures")
        tag = block
        fig_paths = {
            "qc_violins": qc_violin_panel(
                adata,
                keys=["total_counts", "n_genes_by_counts", "percent_mito"],
                path=fig_dir / f"module1_{tag}_qc_violins.png",
            ),
            "counts_vs_genes": qc_scatter(
                adata,
                x="total_counts",
                y="n_genes_by_counts",
                path=fig_dir / f"module1_{tag}_qc_counts_vs_genes.png",
            ),
            "counts_vs_mito": qc_scatter(
                adata,
                x="total_counts",
                y="percent_mito",
                path=fig_dir / f"module1_{tag}_qc_counts_vs_mito.png",
            ),
        }

        adata, filter_log = apply_filters(adata, cfg)
        print(filter_log)
        adata = normalize_log(adata, cfg)
        adata = hvg_pca(adata, cfg)
        adata, sub_info = maybe_subsample(adata, cfg)
        filter_log = append_subsample_filter_row(filter_log, sub_info, adata)
        adata = neighbors_umap(adata, cfg)

        colors = [
            c
            for c in ["total_counts", "percent_mito", "azimuth_label", "sample_id"]
            if c in adata.obs.columns
        ]
        fig_paths["umap"] = umap_panel(adata, colors=colors, path=fig_dir / f"module1_{tag}_umap.png")

        # Snapshot timing before run_params is written (finally alone is too late).
        finalize_timing(timing)
        out = save_module1_outputs(
            cfg,
            adata,
            filter_log=filter_log,
            metrics_summary=qc_metrics_summary(adata),
            figure_paths=fig_paths,
            run_extras={
                "label_join": join_info,
                "subsample": sub_info,
                "donor_label": donor,
                **{k: v for k, v in timing.items() if not str(k).startswith("_t")},
            },
        )
    print(
        "Module 1 QC OK",
        donor,
        block,
        f"_compute_seconds={timing.get('_compute_seconds')}",
        f"peak_rss_mb={timing.get('peak_rss_mb')}",
    )
    for k, p in out.items():
        print(f"  {k}: {p} exists={Path(p).exists()}")
    return Path(out["h5ad"])


def ensure_module3_inputs(force: bool = False) -> None:
    from scripts.common.paths import load_config
    from scripts.nb02_qc.export import resolve_qc_h5ad

    cfg = load_config()
    for item in cfg["module2"]["inputs"]:
        out = resolve_qc_h5ad(cfg, item["block_id"])
        if out.exists() and not force:
            print(f"Module 1 QC input present: {out}")
            continue
        print(f"Running Module 1 QC for {item}")
        run_module2(block_id=item["block_id"], donor_label=item["donor_label"])


def run_module3() -> None:
    """Learner Module 2 -- integration half (CLI --module 3 for historical script numbers)."""
    from scripts.common.paths import load_config, resolve
    from scripts.common.plotting import umap_panel
    from scripts.common.runtime import measure_run, rss_checkpoint
    from scripts.nb03_integration.concatenate import concatenate_donors
    from scripts.nb03_integration.compare import (
        composition_weighted_pseudobulk,
        hubmap_marker_matrix,
        max_across_labels_pseudobulk,
        plot_hubmap_marker_heatmap,
    )
    from scripts.nb03_integration.export import save_module2_integration_outputs
    from scripts.nb03_integration.integrate import (
        donor_composition,
        ensure_cluster_labels,
        harmony_mixing_metrics,
        harmony_silhouette_by_celltype,
        leiden_resolution_sweep,
        post_harmony_label_views,
        recompute_pca_umap,
        run_harmony,
    )

    with measure_run() as timing:
        cfg = load_config()
        fig_dir = resolve(cfg, "outputs_figures")
        diag = cfg["module2"].get("diagnostics") or {}
        adata, donor_summary, _shared_ids = concatenate_donors(cfg)
        print(donor_summary)
        rss_checkpoint("after_concat", adata)
        adata = recompute_pca_umap(adata, cfg)
        rss_checkpoint("after_pca_umap", adata)
        fig_paths = {
            "umap_before": umap_panel(
                adata,
                colors=[c for c in ["donor_id", "donor_label"] if c in adata.obs.columns],
                path=fig_dir / "module2_umap_before_harmony.png",
            )
        }
        if cfg["module2"]["harmony"].get("enabled", True):
            adata = run_harmony(adata, cfg)
            adata.obsm["X_umap"] = adata.obsm["X_umap_harmony"]
            rss_checkpoint("after_harmony", adata)
        comparison_colors = post_harmony_label_views(adata, cfg)
        fig_paths["umap_after"] = umap_panel(
            adata,
            colors=[c for c in ["donor_id", "donor_label"] if c in adata.obs.columns],
            path=fig_dir / "module2_umap_after_harmony.png",
        )
        fig_paths["umap_label_comparison"] = umap_panel(
            adata,
            colors=comparison_colors,
            path=fig_dir / "module2_umap_label_comparison.png",
            ncols=2,
        )

        label_col = ensure_cluster_labels(adata, cfg)
        # Leiden sweep after neighbors exist (Harmony neighbors preferred)
        sweep = leiden_resolution_sweep(adata, cfg)
        print("Leiden sweep:\n", sweep)

        mix = None
        sil = None
        if diag.get("enabled", True):
            mix = harmony_mixing_metrics(
                adata,
                k=int(diag.get("knn_k", 30)),
                max_cells=int(diag.get("max_cells", 8000)),
                random_seed=int(
                    diag.get("random_seed")
                    if diag.get("random_seed") is not None
                    else cfg["module2"].get("random_seed", 0)
                ),
            )
            sil = harmony_silhouette_by_celltype(
                adata,
                label_col="cell_class",
                max_cells=int(diag.get("max_cells", 8000)),
                random_seed=int(
                    diag.get("random_seed")
                    if diag.get("random_seed") is not None
                    else cfg["module2"].get("random_seed", 0)
                ),
            )
            rss_checkpoint("after_diagnostics", adata)
            print("Harmony mixing:\n", mix)
        else:
            print("Skipping Harmony mixing/silhouette (module2.diagnostics.enabled=false)")
            rss_checkpoint("diagnostics_skipped", adata)

        crosstab = donor_composition(adata, label_col=label_col)
        weighted, fracs, _means = composition_weighted_pseudobulk(adata, label_col=label_col)
        max_vec = max_across_labels_pseudobulk(adata, label_col=label_col, linear=True)
        marker_mat = hubmap_marker_matrix(
            adata, cfg["module2"]["markers"], label_col=label_col, top_n_labels=12
        )
        fig_paths["marker_heatmap"] = plot_hubmap_marker_heatmap(
            marker_mat, fig_dir / "module2_hubmap_marker_heatmap.png"
        )

        from scripts.common.runtime import finalize_timing

        finalize_timing(timing)
        out = save_module2_integration_outputs(
            cfg,
            adata,
            donor_summary=donor_summary,
            donor_label_crosstab=crosstab,
            composition_weighted=weighted,
            composition_fractions=fracs,
            max_across_labels=max_vec,
            marker_matrix=marker_mat,
            harmony_mixing=mix,
            harmony_silhouette=sil,
            leiden_sweep=sweep,
            figure_paths=fig_paths,
            extras={
                "label_col": label_col,
                "cell_class_unmapped": adata.uns.get("cell_class_unmapped_labels", []),
                "cell_class_mapping_n": adata.uns.get("cell_class_mapping_n", {}),
                **{k: v for k, v in timing.items() if not str(k).startswith("_t")},
            },
        )
    print(
        "Module 2 integration OK (CLI module 3)",
        f"_compute_seconds={timing.get('_compute_seconds')}",
        f"peak_rss_mb={timing.get('peak_rss_mb')}",
    )
    for k, p in out.items():
        print(f"  {k}: {p} exists={Path(p).exists()}")


def run_module4() -> None:
    """Learner Module 2 -- DE / pathways half (CLI --module 4)."""
    import pandas as pd

    from scripts.common.paths import load_config, resolve
    from scripts.common.runtime import measure_run
    from scripts.nb04_de.enrich import run_enrichment_for_groups, run_prerank_gsea
    from scripts.nb04_de.export import save_module2_de_outputs
    from scripts.nb04_de.load import gene_symbol_series, load_integrated
    from scripts.nb04_de.markers import (
        drop_ribosomal_markers,
        filter_markers_for_enrichment,
        groups_below_min_cells,
        known_marker_presence,
        rank_genes_to_frame,
        run_azimuth_focus_de,
        run_rank_genes,
        run_rank_genes_for_prerank,
        wilcoxon_score_ranks_for_prerank,
    )
    from scripts.nb04_de.plots import (
        plot_enrichment_bar,
        plot_known_marker_dotplot,
        plot_marker_dotplot,
    )

    with measure_run() as timing:
        cfg = load_config()
        fig_dir = resolve(cfg, "outputs_figures")
        adata, input_path = load_integrated(cfg)
        print("Module 2 DE input:", input_path, adata.shape)
        groupby = cfg["module2"]["groupby"]
        min_cells = int(cfg["module2"]["de"].get("min_cells_per_group", 50))
        dropped = groups_below_min_cells(adata, groupby, min_cells)
        print(dropped)

        rank_key = run_rank_genes(adata, cfg, groupby=groupby)
        markers = drop_ribosomal_markers(rank_genes_to_frame(adata, rank_key), cfg)
        az_key, az_markers = run_azimuth_focus_de(adata, cfg)
        if az_key and not az_markers.empty:
            az_markers = drop_ribosomal_markers(az_markers, cfg)
        known = known_marker_presence(adata, cfg)

        fig_paths = {}
        p = plot_known_marker_dotplot(
            adata, cfg, groupby=groupby, path=fig_dir / "module2_known_marker_dotplot.png"
        )
        if p:
            fig_paths["known_marker_dotplot"] = p
        p = plot_marker_dotplot(
            adata, markers, groupby=groupby, path=fig_dir / "module2_marker_dotplot.png"
        )
        if p:
            fig_paths["marker_dotplot"] = p

        markers_for_enr = filter_markers_for_enrichment(markers, cfg)
        background = (
            gene_symbol_series(adata)
            .astype(str)
            .str.strip()
            .replace({"": pd.NA, "nan": pd.NA, "None": pd.NA, "NAN": pd.NA})
            .dropna()
            .unique()
            .tolist()
        )
        background = [g for g in background if not str(g).startswith("ENSG")]
        print(f"Enrichment background genes: {len(background)}")
        enrichment = run_enrichment_for_groups(markers_for_enr, cfg, background_genes=background)
        p = plot_enrichment_bar(enrichment, fig_dir / "module2_enrichment_summary.png")
        if p:
            fig_paths["enrichment_summary"] = p

        prerank_frames = []
        kept_groups = dropped.loc[dropped["kept"], "group"].astype(str).tolist()
        prerank_key = run_rank_genes_for_prerank(adata, cfg, groupby=groupby)
        for g in kept_groups:
            ranked = wilcoxon_score_ranks_for_prerank(adata, prerank_key, g, cfg=cfg)
            if ranked.empty or len(ranked) < 50:
                continue
            pre = run_prerank_gsea(ranked, cfg)
            if pre is not None and not pre.empty:
                pre = pre.copy()
                pre.insert(0, "group", g)
                prerank_frames.append(pre)
        prerank_df = pd.concat(prerank_frames, ignore_index=True) if prerank_frames else pd.DataFrame()

        from scripts.common.runtime import finalize_timing

        finalize_timing(timing)
        out = save_module2_de_outputs(
            cfg,
            marker_df=markers,
            enrichment_df=enrichment,
            known_markers_df=known,
            figure_paths=fig_paths,
            prerank_df=prerank_df,
            dropped_groups=dropped,
            extras={
                "rank_key": rank_key,
                "prerank_rank_key": prerank_key,
                "azimuth_focus_key": az_key,
                "n_azimuth_focus_markers": int(az_markers.shape[0]),
                "n_markers_for_enrichment": int(markers_for_enr.shape[0]),
                "n_background_genes": len(background),
                "n_ribosomal_removed_from_marker_table": int(markers.attrs.get("n_ribosomal_removed", 0)),
                "n_ribosomal_removed_before_enrichment": int(
                    markers_for_enr.attrs.get("n_ribosomal_removed", 0)
                ),
                "exclude_ribosomal_genes_de": bool(cfg["module2"]["de"].get("exclude_ribosomal_genes", True)),
                "exclude_ribosomal_genes": bool(
                    markers_for_enr.attrs.get("exclude_ribosomal_genes", True)
                ),
                "hugo_symbol_gtex_fill": adata.uns.get("hugo_symbol_gtex_fill"),
                "feature_space": adata.uns.get("module2_feature_space"),
                **{k: v for k, v in timing.items() if not str(k).startswith("_t")},
            },
        )
    print(
        "Module 2 DE/pathways OK (CLI module 4)",
        f"_compute_seconds={timing.get('_compute_seconds')}",
        f"peak_rss_mb={timing.get('peak_rss_mb')}",
    )
    for k, pth in out.items():
        print(f"  {k}: {pth} exists={Path(pth).exists()}")


def run_module5() -> None:
    """Learner Module 3 -- multiome / MOFA (CLI --module 5)."""
    import json
    from pathlib import Path as _Path

    from scripts.common.paths import load_config, resolve, portable_path
    from scripts.common.runtime import measure_run
    from scripts.nb05_multimodal.bridge import (
        azimuth_composition,
        barcode_overlap_from_mofa,
        barcode_overlap_table,
        cluster_concordance,
        combine_focus_correlations,
        factor_label_association,
        focus_gene_availability,
        focus_gene_bridge,
        focus_gene_correlations,
        focus_gene_correlations_by_cell_type,
        focus_from_paired_matrix,
        modality_contrast_table,
        mofa_factor_dominance,
        qc_association,
        rna_atac_gene_bridge,
    )
    from scripts.nb05_multimodal.export import save_module3_outputs
    from scripts.nb05_multimodal.inventory import inventory_module3_assets
    from scripts.nb05_multimodal.load import (
        adata_from_mofa_factors,
        assemble_nucleus_frame,
        focus_gene_matrix,
        focus_gene_matrix_from_mofa,
        load_joint_embedding,
        load_mofa,
        mofa_nondominant_r2_summary,
        multiome_paths,
        recompute_mofa_r2_total,
    )
    from scripts.nb05_multimodal.plots import (
        plot_cluster_concordance,
        plot_factor_label_heatmap,
        plot_focus_correlations,
        plot_focus_gene_bars,
        plot_mofa_umap,
        plot_mofa_variance,
        plot_qc_association,
        plot_rna_atac_scatter,
        plot_umap_panel,
        plot_wnn_umap,
    )

    from scripts.common.runtime import peak_rss_mb
    import time as _time

    cfg = load_config()
    _t0 = _time.perf_counter()
    _rss0 = peak_rss_mb()
    fig_dir = resolve(cfg, "outputs_figures")
    paths = multiome_paths(cfg)
    donor = cfg["module3"]["donor_label"]
    load_h5mu = bool(cfg["module3"].get("params", {}).get("load_h5mu", False))
    noise_delta = float((cfg["module3"].get("thresholds") or {}).get("noise_delta_r", 0.01))
    min_n = int((cfg["module3"].get("thresholds") or {}).get("min_nuclei_correlation", 50))
    max_fac = int((cfg["module3"].get("analysis") or {}).get("mofa_variance_max_factor", 14))

    inventory = inventory_module3_assets(cfg)
    print(inventory[["file_name", "exists", "size_mb", "core_or_optional"]])

    mofa = load_mofa(cfg)
    r2_recompute = recompute_mofa_r2_total(cfg, view="rna")
    print("R2 recompute:", r2_recompute.get("label"), r2_recompute)
    nondom = mofa_nondominant_r2_summary(
        mofa["variance_per_factor"], r2_recompute=r2_recompute
    )
    print(
        f"Max non-dominant R^2 = {nondom.get('max_nondominant_r2'):.5f} "
        f"({nondom.get('max_nondominant_percent'):.4f} percent; "
        f"units={nondom.get('r2_units')})"
    )

    contrast = modality_contrast_table(cfg)
    if load_h5mu and paths["secondary_h5mu"].exists():
        overlap = barcode_overlap_table(cfg)
        adata = load_joint_embedding(cfg)
        if "X_mofa" not in adata.obsm:
            adata = adata[adata.obs_names.isin(mofa["factors"].index)].copy()
            adata.obsm["X_mofa"] = mofa["factors"].loc[adata.obs_names].to_numpy()
        focus_long = focus_gene_matrix(cfg, mofa)
        layer_meta = {
            "rna_layer": focus_long.attrs.get("rna_layer"),
            "rna_layer_encoding": focus_long.attrs.get("rna_layer_encoding"),
            "rna_layer_layout": focus_long.attrs.get("rna_layer_layout"),
            "available_rna_layers": focus_long.attrs.get("available_rna_layers"),
            "source": focus_long.attrs.get("source"),
        }
    else:
        overlap = barcode_overlap_from_mofa(mofa)
        adata = adata_from_mofa_factors(mofa)
        focus_long = focus_gene_matrix_from_mofa(cfg, mofa)
        layer_meta = {
            "rna_layer": focus_long.attrs.get("rna_layer"),
            "source": focus_long.attrs.get("source"),
            "note": focus_long.attrs.get("note"),
            "load_h5mu": False,
        }

    nucleus = assemble_nucleus_frame(cfg)
    for col in ("azimuth_label", "leiden_wnn", "rna_leiden", "atac_leiden"):
        if col in nucleus.columns:
            adata.obs[col] = nucleus.reindex(adata.obs_names)[col].to_numpy()
    if not focus_long.empty and "azimuth_label" in nucleus.columns and "azimuth_label" not in focus_long.columns:
        focus_long = focus_long.merge(
            nucleus.reset_index()[["barcode", "azimuth_label"]],
            on="barcode",
            how="left",
        )
    conc = cluster_concordance(nucleus)
    qc_df = qc_association(
        nucleus,
        label_cols=["atac_leiden", "rna_leiden"],
        qc_cols=["nFrags", "ReadsInTSS", "DoubletEnrichment", "PromoterRatio"],
    )
    assoc = factor_label_association(
        mofa["factors"],
        nucleus,
        label_cols=["azimuth_label", "leiden_wnn", "rna_leiden", "atac_leiden"],
    )
    comp = azimuth_composition(nucleus)

    # Attach recorded layer probe if present (dev-machine Q1)
    probe_path = resolve(cfg, "outputs_reports") / "module3_rna_layer_probe.json"
    layer_probe = None
    if probe_path.exists():
        layer_probe = json.loads(probe_path.read_text())

    print(overlap)
    fig_paths = {}
    if load_h5mu:
        p = plot_umap_panel(adata, fig_dir / "module3_umap_modalities.png")
        if p:
            fig_paths["umap_modalities"] = p
    p = plot_mofa_umap(adata, fig_dir / "module3_mofa_latent.png", donor_label=donor)
    if p:
        fig_paths["mofa_latent"] = p
    p = plot_mofa_variance(
        mofa["variance_total"],
        mofa["variance_per_factor"],
        fig_dir / "module3_mofa_variance.png",
        n_factors=max_fac,
        r2_recompute=r2_recompute,
    )
    if p:
        fig_paths["mofa_variance"] = p
    p = plot_qc_association(qc_df, fig_dir / "module3_qc_association.png")
    if p:
        fig_paths["qc_association"] = p
    p = plot_factor_label_heatmap(
        assoc, fig_dir / "module3_factor_label_heatmap.png", n_factors=max_fac
    )
    if p:
        fig_paths["factor_label_heatmap"] = p
    p = plot_wnn_umap(nucleus, fig_dir / "module3_wnn_umap.png", donor_label=donor)
    if p:
        fig_paths["wnn_umap"] = p
    p = plot_cluster_concordance(nucleus, fig_dir / "module3_cluster_concordance.png")
    if p:
        fig_paths["cluster_concordance"] = p

    bridge = rna_atac_gene_bridge(mofa)
    focus = focus_from_paired_matrix(focus_long)
    if focus.empty:
        focus = focus_gene_bridge(bridge, cfg)
    corr_global = focus_gene_correlations(focus_long, min_nuclei=min_n)
    corr_cell = focus_gene_correlations_by_cell_type(focus_long, cfg)
    corr = combine_focus_correlations(corr_global, corr_cell)
    availability = focus_gene_availability(cfg, focus_long, bridge)
    factor_dom = mofa_factor_dominance(mofa["variance_per_factor"], n_factors=max_fac)
    print(
        "shared genes",
        bridge.shape[0],
        "metric",
        bridge.attrs.get("bridge_metric"),
        "r",
        bridge.attrs.get("pearson_log1p"),
        "rho",
        bridge.attrs.get("spearman_log1p"),
    )
    p = plot_rna_atac_scatter(
        bridge, fig_dir / "module3_rna_atac_scatter.png", focus_df=focus, donor_label=donor
    )
    if p:
        fig_paths["rna_atac_scatter"] = p
    p = plot_focus_gene_bars(focus, fig_dir / "module3_focus_gene_bars.png")
    if p:
        fig_paths["focus_gene_bars"] = p
    p = plot_focus_correlations(
        corr, fig_dir / "module3_focus_gene_correlations.png", noise_delta=noise_delta
    )
    if p:
        fig_paths["focus_correlations"] = p

    for stale_name in (
        "module5_umap_rna_vs_atac.png",
        "module5_mofa_variance.png",
        "module5_focus_gene_correlations.png",
    ):
        stale = fig_dir / stale_name
        if stale.exists():
            stale.unlink()

    out = save_module3_outputs(
        cfg,
        inventory_df=inventory,
        contrast_df=contrast,
        overlap_df=overlap,
        bridge_df=bridge,
        focus_df=focus,
        corr_df=corr,
        variance_total=mofa["variance_total"],
        variance_per_factor=mofa["variance_per_factor"],
        figure_paths=fig_paths,
        extras={
            "multiome_dir": portable_path(cfg, paths["dir"]),
            "secondary_h5mu": str(paths["secondary_h5mu"]),
            "mofa_hdf5": portable_path(cfg, paths["mofa_hdf5"]),
            "n_nuclei": int(adata.n_obs),
            "n_mofa_factors": int(mofa["factors"].shape[1]),
            "portal_url": cfg["module3"].get("portal_url"),
            "load_h5mu": load_h5mu,
            "layer_meta": layer_meta,
            "layer_probe_path": str(probe_path) if probe_path.exists() else None,
            "layer_probe_summary": {
                "NKG7_layer_sums": (layer_probe or {}).get("NKG7_layer_sums"),
                "rna_layers": (layer_probe or {}).get("rna_layers"),
                "layout_spliced": (layer_probe or {}).get("layer_spliced_inferred_layout"),
            }
            if layer_probe
            else None,
            "_compute_seconds": round(_time.perf_counter() - _t0, 3),
            "peak_rss_mb": peak_rss_mb(),
            "peak_rss_mb_start": _rss0,
        },
        availability_df=availability,
        factor_dom_df=factor_dom,
        cluster_concordance_df=conc,
        qc_association_df=qc_df,
        factor_label_association_df=assoc,
        azimuth_composition_df=comp,
        r2_recompute=r2_recompute,
        nondom_summary=nondom,
    )
    print("Module 3 OK (CLI module 5)")
    for k, pth in out.items():
        print(f"  {k}: {pth} exists={_Path(pth).exists()}")


def run_module6() -> None:
    """Learner Module 4 -- cross-ecosystem bulk rebuild (CLI --module 6)."""
    from pathlib import Path as _Path

    import time as _time

    import pandas as pd
    from scipy.stats import spearmanr

    from scripts.common.paths import load_config, module_root, resolve
    from scripts.common.runtime import peak_rss_mb, rss_checkpoint
    from scripts.nb06_crossres.compare import (
        compare_configured_nes_pairs,
        hubmap_gtex_bridge_table,
        is_ercc_symbol,
        nes_sensitivity_table,
        plot_nes_scatter,
        plot_nes_sensitivity,
        plot_pseudobulk_vs_gtex,
        plot_volcano,
        prerank_contrast,
        resolve_nes_sensitivity,
    )
    from scripts.nb06_crossres.contrasts import (
        contrast_h_aging_gtex,
        contrast_h_aging_hubmap,
        contrast_h_disease,
        contrast_m_flight,
        gtex_age_bracket_inventory,
        logfc_series,
    )
    from scripts.nb06_crossres.ercc_calib import (
        ercc_calibration_frame,
        ercc_calibration_summary,
        plot_ercc_calibration,
    )
    from scripts.nb06_crossres.export import save_module4_outputs
    from scripts.nb06_crossres.genelab_subset import extract_genelab_isst_de_subset
    from scripts.nb06_crossres.load import native_levels_panel, sample_counts_table
    from scripts.nb06_crossres.orthologs import load_ortholog_table, ortholog_loss_record

    cfg = load_config()
    _t0 = _time.perf_counter()
    _rss0 = peak_rss_mb()
    _rss_checkpoints: list[dict] = []
    def _rss(label: str) -> None:
        _rss_checkpoints.append(rss_checkpoint(label))

    _rss("m4_start")
    fig_dir = resolve(cfg, "outputs_figures")
    panel = list(cfg["module4"]["params"]["gene_panel"])
    highlight = list(cfg["module4"]["params"].get("highlight_pathways") or [])

    print("=== GTEx age bracket inventory (declare before DE) ===")
    gtex_inv = gtex_age_bracket_inventory(cfg)
    print(gtex_inv)

    print("=== sample counts ===")
    sample_counts = sample_counts_table(cfg)
    print(sample_counts)

    print("=== native levels (deliberately incomparable) ===")
    native = native_levels_panel(cfg)
    print(native[["gene_symbol_human", "gtex_lung_mean_tpm", "gse150910_mean_counts_ipf"]].head())
    _rss("after_native_levels")

    print("=== H_DISEASE (pydeseq2) ===")
    de_h = contrast_h_disease(cfg)
    print(de_h.sort_values("padj").head(5)[["gene_symbol", "logFC", "padj"]])
    _rss("after_h_disease")

    print("=== H_AGING_GTEX (pydeseq2 on gene_reads; predeclared 20-39 vs 60-79) ===")
    de_age_gtex = contrast_h_aging_gtex(cfg)
    print(de_age_gtex.sort_values("padj").head(5)[["gene_symbol", "logFC", "padj"]])
    print("design", getattr(de_age_gtex, "attrs", {}).get("design"))
    _rss("after_h_aging_gtex")

    print("=== M_FLIGHT (pydeseq2; formerly M_UNLOAD) ===")
    from scripts.nb06_crossres.ercc_calib import (
        ercc_read_fraction_arm_summary,
        ercc_read_fraction_by_sample,
    )

    ercc_frac = ercc_read_fraction_by_sample(cfg)
    ercc_frac_path = resolve(cfg, "outputs_tables") / "module4_ercc_read_fraction_by_sample.tsv"
    ercc_frac.to_csv(ercc_frac_path, sep="\t", index=False)
    ercc_frac_summary = ercc_read_fraction_arm_summary(ercc_frac)
    print("ERCC read fraction by arm:", ercc_frac_summary)

    de_m = contrast_m_flight(cfg, drop_low_qa=False)
    de_m_sens = contrast_m_flight(cfg, drop_low_qa=True)
    print(
        "size_factor polarity extras:",
        (getattr(de_m, "attrs", {}) or {}).get("polarity"),
    )
    print(de_m.sort_values("padj").head(5)[["gene_symbol", "logFC", "padj", "mouse_ensembl"]])
    de_m_no_ercc = de_m[~de_m["gene_symbol"].map(is_ercc_symbol)]
    print("=== M_FLIGHT top genes excluding ERCC- spike-ins ===")
    print(de_m_no_ercc.sort_values("padj").head(5)[["gene_symbol", "logFC", "padj", "mouse_ensembl"]])
    _rss("after_m_flight")

    print("=== H_AGING_HUBMAP (2 vs 2 illustration; summed/rounded counts layer) ===")
    de_age_hub = contrast_h_aging_hubmap(cfg)
    n_sig = int((pd.to_numeric(de_age_hub["padj"], errors="coerce") < 0.05).sum())
    print(de_age_hub.sort_values("padj").head(5)[["gene_symbol", "logFC", "padj"]])
    print("padj<0.05:", n_sig, "(expect few/none; do not loosen thresholds)")
    print("design", getattr(de_age_hub, "attrs", {}).get("design"))
    _rss("after_h_aging_hubmap")

    print("=== ERCC ExFold calibration (flight=Mix1, GC=Mix2) ===")
    ercc_cal = ercc_calibration_frame(de_m, cfg)
    ercc_stats = ercc_calibration_summary(ercc_cal)
    print(ercc_stats)

    print("=== GeneLab ISS-T DE subset (committed; not the 111 MB file) ===")
    _subset_rel = ((cfg.get("module4") or {}).get("data") or {}).get(
        "osd248_genelab_isst_de_subset",
        "outputs/tables/module4_osd248_genelab_isst_de_subset.tsv",
    )
    subset_path = module_root(cfg) / _subset_rel
    if not subset_path.exists():
        subset_path = extract_genelab_isst_de_subset(cfg)
    genelab_subset = pd.read_csv(subset_path, sep="\t")
    print(subset_path, "rows", len(genelab_subset))
    # Cross-check log2FC on shared gene symbols (not ERCC spike-ins)
    gl = genelab_subset.copy()
    gl["sym"] = gl["SYMBOL"].astype(str).str.upper()
    mjoin = de_m.copy()
    mjoin["sym"] = mjoin["gene_symbol"].astype(str).str.upper()
    j = mjoin.merge(gl[["sym", "log2fc_isst_flight_vs_gc"]], on="sym", how="inner")
    j = j.dropna(subset=["logFC", "log2fc_isst_flight_vs_gc"])
    j = j[~j["sym"].str.startswith("ERCC-")]
    genelab_rho = (
        float(spearmanr(j["logFC"], j["log2fc_isst_flight_vs_gc"]).correlation)
        if len(j) >= 3
        else float("nan")
    )
    print(f"GeneLab vs pydeseq2 log2FC Spearman rho={genelab_rho:.3f} n={len(j)}")

    print("=== ortholog loss (M_FLIGHT only) ===")
    ortho = load_ortholog_table(cfg)
    loss = ortholog_loss_record(
        de_h["gene_symbol"].astype(str),
        de_m["human_symbol"].dropna().astype(str),
        ortho,
    )
    print(loss)

    print("=== prerank GSEA for four contrasts (Hallmark + Reactome) ===")
    m_for_gsea = de_m.dropna(subset=["human_symbol"]).copy()
    m_for_gsea["gene_symbol"] = m_for_gsea["human_symbol"].astype(str).str.upper()
    nes_by_id = {
        "H_DISEASE": prerank_contrast(logfc_series(de_h), cfg, contrast_id="H_DISEASE"),
        "H_AGING_GTEX": prerank_contrast(
            logfc_series(de_age_gtex), cfg, contrast_id="H_AGING_GTEX"
        ),
        "M_FLIGHT": prerank_contrast(
            logfc_series(m_for_gsea), cfg, contrast_id="M_FLIGHT"
        ),
        "H_AGING_HUBMAP": prerank_contrast(
            logfc_series(de_age_hub), cfg, contrast_id="H_AGING_HUBMAP"
        ),
    }
    for cid, nes in nes_by_id.items():
        print(f"  {cid}: {len(nes)} pathways")
    _rss("after_nes_prerank")

    print("=== four NES pairs ===")
    nes_pair_summary, nes_pair_details = compare_configured_nes_pairs(nes_by_id, cfg)
    print(nes_pair_summary)

    headline = list(
        (cfg.get("module4", {}).get("params") or {}).get("nes_headline_pair")
        or ["H_DISEASE", "H_AGING_GTEX"]
    )
    headline_key = f"{headline[0]}__{headline[1]}"
    nes_joined, nes_stats = nes_pair_details[headline_key]

    print("=== NES sensitivity on headline pair (includes GO BP) ===")
    from scripts.nb06_crossres.compare import sensitivity_enabled

    _de_map = {
        "H_AGING_GTEX": de_age_gtex,
        "M_FLIGHT": m_for_gsea,
        "H_AGING_HUBMAP": de_age_hub,
        "H_DISEASE": de_h,
    }
    if sensitivity_enabled(cfg):
        print("sensitivity.enabled=true: re-preranking Hallmark+Reactome+GO BP (high RAM)")
        nes_sens, sens_meta = resolve_nes_sensitivity(
            cfg,
            logfc_x=logfc_series(_de_map[headline[0]]),
            logfc_y=logfc_series(_de_map[headline[1]]),
            contrast_x=headline[0],
            contrast_y=headline[1],
        )
    else:
        print("sensitivity.enabled=false: loading committed module4_nes_sensitivity.tsv")
        nes_sens, sens_meta = resolve_nes_sensitivity(cfg)
    print(sens_meta)
    print(nes_sens)
    _rss("after_nes_sensitivity")

    print("=== HuBMAP<->GTEx level bridge ===")
    bridge = hubmap_gtex_bridge_table(cfg)
    _rss("after_bridge")

    fig_paths: dict[str, _Path] = {}
    for de, name, title in (
        (
            de_h,
            "volcano_h_disease",
            "H_DISEASE: GSE150910 IPF vs control (pydeseq2, ~ condition only)",
        ),
        (
            de_age_gtex,
            "volcano_h_aging_gtex",
            "H_AGING_GTEX: older vs younger lung (gene_reads; 20-39 vs 60-79)",
        ),
        (
            de_m,
            "volcano_m_flight",
            "M_FLIGHT: OSD-248 ISS-T flight vs GC (n=10 vs 10)",
        ),
        (
            de_age_hub,
            "volcano_h_aging_hubmap",
            "H_AGING_HUBMAP: older vs younger donors (2 vs 2 illustration only)",
        ),
    ):
        p = plot_volcano(de, fig_dir / f"module4_{name}.png", title=title, panel=panel)
        if p:
            fig_paths[name] = p

    p = plot_ercc_calibration(
        ercc_cal, fig_dir / "module4_ercc_exfold_calibration.png", summary=ercc_stats
    )
    if p:
        fig_paths["ercc_calibration"] = p

    for key, (joined, stats) in nes_pair_details.items():
        safe = key.replace("__", "_vs_")
        p = plot_nes_scatter(
            joined,
            fig_dir / f"module4_nes_scatter_{safe}.png",
            stats=stats,
            highlight=highlight,
            pair_label=stats.get("pair_label"),
        )
        if p:
            fig_paths[f"nes_scatter_{safe}"] = p
        if key == headline_key:
            p2 = plot_nes_scatter(
                joined,
                fig_dir / "module4_nes_scatter.png",
                stats=stats,
                highlight=highlight,
                pair_label=stats.get("pair_label"),
            )
            if p2:
                fig_paths["nes_scatter"] = p2

    p = plot_nes_sensitivity(nes_sens, fig_dir / "module4_nes_sensitivity.png")
    if p:
        fig_paths["nes_sensitivity"] = p
    p = plot_pseudobulk_vs_gtex(
        bridge, fig_dir / "module4_pseudobulk_vs_gtex.png", panel=panel
    )
    if p:
        fig_paths["pseudobulk_vs_gtex"] = p
    _rss("after_figures")

    out = save_module4_outputs(
        cfg,
        sample_counts=sample_counts,
        native_levels=native,
        contrast_h=de_h,
        contrast_m=de_m,
        contrast_m_sens=de_m_sens,
        contrast_aging_gtex=de_age_gtex,
        contrast_aging_hubmap=de_age_hub,
        ortholog_loss=loss,
        nes_by_id=nes_by_id,
        nes_pair_summary=nes_pair_summary,
        nes_pair_details=nes_pair_details,
        bridge_df=bridge,
        figure_paths=fig_paths,
        nes_sensitivity=nes_sens,
        ercc_calibration=ercc_cal,
        genelab_subset=genelab_subset,
        gtex_bracket_inventory=gtex_inv,
        extras={
            "n_h_disease": int(de_h.shape[0]),
            "n_h_aging_gtex": int(de_age_gtex.shape[0]),
            "n_m_flight": int(de_m.shape[0]),
            "n_m_flight_qa_filtered": int(de_m_sens.shape[0]),
            "n_h_aging_hubmap": int(de_age_hub.shape[0]),
            "n_h_aging_hubmap_padj_lt_05": n_sig,
            "flight_low_qa_sample": "Mmus_C57-6T_LNG_FLT_ISS-T_Rep7_F7",
            "ercc_calibration": ercc_stats,
            "ercc_read_fraction": ercc_frac_summary,
            "m_flight_polarity": getattr(de_m, "attrs", {}).get("polarity"),
            "rss_checkpoints": _rss_checkpoints,
            "nes_sensitivity_meta": sens_meta,
            "genelab_isst_de_subset": str(subset_path),
            "n_genelab_subset_rows": int(len(genelab_subset)),
            "genelab_logfc_spearman": genelab_rho,
            "n_genelab_joined": int(len(j)),
            "h_aging_gtex_design": getattr(de_age_gtex, "attrs", {}).get("design"),
            "h_aging_hubmap_design": getattr(de_age_hub, "attrs", {}).get("design"),
            "_compute_seconds": round(_time.perf_counter() - _t0, 3),
            "peak_rss_mb": peak_rss_mb(),
            "peak_rss_mb_start": _rss0,
        },
    )

    print("=== shared-mechanism concordance (headline; uses existing NES tables) ===")
    from scripts.nb06_crossres.export import append_mechanism_section
    from scripts.nb06_crossres.mechanisms import run_mechanism_analysis, save_mechanism_outputs

    mech = run_mechanism_analysis(cfg)
    mech_paths = save_mechanism_outputs(cfg, mech)
    append_mechanism_section(
        out["methods"],
        null_df=mech["concordance_null"],
        summary=mech["mechanism_summary"],
        mechanisms=mech["mechanisms"],
    )
    print(mech["concordance_null"][
        ["set_id", "observed_concordant", "observed_concordant_unflagged", "perm_mean", "perm_p_ge_obs"]
    ].to_string(index=False))
    for k, pth in mech_paths.items():
        print(f"  {k}: {pth} exists={_Path(pth).exists()}")
    _rss("after_mechanisms")

    # Stamp mid-run RSS profile onto run_params (save happened before mechanisms).
    rp = _Path(out["run_params"])
    if rp.exists():
        import json as _json

        payload = _json.loads(rp.read_text())
        extras = payload.setdefault("extras", {})
        extras["rss_checkpoints"] = _rss_checkpoints
        extras["peak_rss_mb"] = peak_rss_mb()
        payload["peak_rss_mb"] = peak_rss_mb()
        rp.write_text(_json.dumps(payload, indent=2))

    print("Module 4 rebuild OK (CLI module 6)")
    for k, pth in out.items():
        print(f"  {k}: {pth} exists={_Path(pth).exists()}")
    print("RSS checkpoints:")
    for row in _rss_checkpoints:
        print(f"  {row['label']}: {row['peak_rss_mb']}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--module",
        choices=["1", "2", "3", "4", "5", "6", "all"],
        default="all",
        help="CLI numbers: 3/4=Module2 halves, 5=Module3, 6=Module4 crossres",
    )
    parser.add_argument("--block-id", default=None)
    parser.add_argument("--donor-label", default=None)
    parser.add_argument(
        "--force-module2",
        action="store_true",
        help="Re-run Module 2 inputs even if processed h5ads exist",
    )
    args = parser.parse_args()

    if args.module in {"1", "all"}:
        run_module1()
    if args.module == "2":
        if args.block_id or args.donor_label:
            run_module2(block_id=args.block_id, donor_label=args.donor_label)
        else:
            # Default: all teaching donors (Module 1 QC half)
            from scripts.common.paths import load_config

            cfg = load_config()
            for item in cfg["module2"]["inputs"]:
                run_module2(
                    block_id=item.get("block_id") or item.get("primary_id"),
                    donor_label=item["donor_label"],
                )
    if args.module in {"3", "all"}:
        ensure_module3_inputs(force=args.force_module2)
        run_module3()
    if args.module in {"4", "all"}:
        run_module4()
    if args.module in {"5", "all"}:
        run_module5()
    if args.module in {"6", "all"}:
        run_module6()


if __name__ == "__main__":
    main()
