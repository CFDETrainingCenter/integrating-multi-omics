"""Methods / parameter reporting for reproducible notebook write-ups."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

import pandas as pd

from scripts.common.io import environment_versions


def deep_update(base: dict[str, Any], overrides: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge overrides into a copy of base."""
    out = deepcopy(base)
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_update(out[key], value)
        else:
            out[key] = deepcopy(value)
    return out


def apply_analysis_overrides(cfg: dict[str, Any], overrides: dict[str, Any] | None) -> dict[str, Any]:
    """
    Merge notebook-local overrides into config.

    Typical Module 3 usage::

        OVERRIDES = {"module2": {"embedding": {"umap_min_dist": 0.2, "n_neighbors": 30}}}
        cfg = apply_analysis_overrides(cfg, OVERRIDES)
    """
    if not overrides:
        return cfg
    return deep_update(cfg, overrides)


def module1_embedding_params(cfg: dict[str, Any]) -> dict[str, Any]:
    """QC / single-block preprocess knobs (learner Module 1)."""
    return dict(cfg["module1"]["preprocess"])


def module3_embedding_params(cfg: dict[str, Any]) -> dict[str, Any]:
    """Integration embedding knobs (learner Module 2; name kept for import stability)."""
    base = module1_embedding_params(cfg)
    emb = dict(cfg.get("module2", {}).get("embedding") or {})
    merged = {**base, **emb}
    merged.setdefault("umap_min_dist", 0.3)
    merged.setdefault("umap_spread", 1.0)
    merged.setdefault("scale_max_value", 10)
    merged.setdefault("leiden_resolution", 0.8)
    merged.setdefault("n_pcs_neighbors", merged.get("n_pcs"))
    return merged


def methods_table_module2(cfg: dict[str, Any]) -> pd.DataFrame:
    qc = cfg["module1"]["qc"]
    pp = module1_embedding_params(cfg)
    rows = [
        ("Input matrix", "HuBMAP Salmon `raw_expr.h5ad` (counts-like)", "anndata / scanpy"),
        ("QC metrics", "scanpy.pp.calculate_qc_metrics (mt, ribo)", "scanpy"),
        ("Gene filter", f"min_cells={qc['min_cells_per_gene']}", "scanpy.pp.filter_genes"),
        ("Cell filter: complexity", f"min_genes={qc['min_genes']}", "scanpy.pp.filter_cells"),
        ("Cell filter: depth", f"min_counts={qc['min_counts']}", "obs threshold"),
        ("Cell filter: mito", f"max_mito_percent={qc['max_mito_percent']}", "obs threshold"),
        ("Cell filter: max genes", f"max_genes={qc['max_genes']}", "obs threshold"),
        (
            "Doublet filtering",
            "Not applied in core path (HuBMAP processed nuclei; optional Scrublet deferred)",
            "n/a",
        ),
        ("Normalization", f"normalize_total target_sum={pp['target_sum']}", "scanpy.pp.normalize_total"),
        ("Log transform", "log1p", "scanpy.pp.log1p"),
        ("Feature selection", f"Seurat HVG n_top_genes={pp['n_top_genes']}", "scanpy.pp.highly_variable_genes"),
        ("Scaling", f"max_value={pp.get('scale_max_value', 10)} (HVG subset)", "scanpy.pp.scale"),
        ("PCA", f"n_comps={pp['n_pcs']}, svd_solver=arpack", "scanpy.tl.pca"),
        ("Neighbors graph", f"n_neighbors={pp['n_neighbors']}, n_pcs={pp['n_pcs']}", "scanpy.pp.neighbors"),
        (
            "UMAP",
            f"min_dist={pp.get('umap_min_dist', 0.5)}, spread={pp.get('umap_spread', 1.0)}",
            "scanpy.tl.umap (umap-learn)",
        ),
        ("Teaching subsample", f"max_nuclei={cfg['module1'].get('max_nuclei')}", "numpy RNG"),
        ("Random seed", str(cfg["module1"].get("random_seed", 0)), "config"),
    ]
    return pd.DataFrame(rows, columns=["step", "parameters", "package_or_function"])


def methods_table_module3(cfg: dict[str, Any]) -> pd.DataFrame:
    """Integration methods table (learner Module 2; function name kept for imports)."""
    emb = module3_embedding_params(cfg)
    m2 = cfg["module2"]
    harm = m2.get("harmony", {})
    n_pcs_neighbors = emb.get("n_pcs_neighbors") or emb["n_pcs"]
    rows = [
        ("Inputs", "Module 1 QC/preprocessed h5ads (counts layer retained)", "anndata"),
        ("Gene ID harmonization", "strip Ensembl version; intersect genes", "pandas"),
        ("Concatenation", "inner gene intersect across donors (join='inner')", "anndata.concat"),
        ("HVG (joint)", f"Seurat flavor n_top_genes={emb['n_top_genes']}", "scanpy.pp.highly_variable_genes"),
        ("Scaling", f"max_value={emb.get('scale_max_value', 10)}", "scanpy.pp.scale"),
        ("PCA", f"n_comps={emb['n_pcs']}, svd_solver=arpack", "scanpy.tl.pca"),
        (
            "Neighbors (uncorrected)",
            f"use_rep=X_pca, n_neighbors={emb['n_neighbors']}, n_pcs={n_pcs_neighbors}",
            "scanpy.pp.neighbors",
        ),
        (
            "UMAP (uncorrected)",
            f"min_dist={emb['umap_min_dist']}, spread={emb['umap_spread']}",
            "scanpy.tl.umap",
        ),
        (
            "Batch correction",
            (
                f"Harmony on X_pca; key={harm.get('key')}; "
                f"max_iter_harmony={harm.get('max_iter_harmony', 20)}; "
                f"random_state={m2.get('random_seed', 0)}; "
                f"theta={harm.get('theta', 'default')}"
            ),
            "harmonypy.run_harmony",
        ),
        (
            "Neighbors (Harmony)",
            f"use_rep=X_pca_harmony, n_neighbors={emb['n_neighbors']}, n_pcs=all Harmony PCs",
            "scanpy.pp.neighbors",
        ),
        (
            "UMAP (Harmony)",
            f"min_dist={emb['umap_min_dist']}, spread={emb['umap_spread']}",
            "scanpy.tl.umap",
        ),
        (
            "Leiden (diagnostic)",
            (
                f"Always computed for label-comparison panels "
                f"(resolution={emb.get('leiden_resolution', 0.8)}); "
                "Azimuth preferred for biological interpretation"
            ),
            "scanpy.tl.leiden (leidenalg)",
        ),
        ("Labels for pseudobulk", m2.get("label_column", "azimuth_label"), "Azimuth join from Module 1"),
        ("GTEx reference", "Deferred to Module 4 (cross-ecosystem)", "n/a in Module 2"),
        ("Random seed", str(m2.get("random_seed", 0)), "config"),
    ]
    return pd.DataFrame(rows, columns=["step", "parameters", "package_or_function"])


def methods_narrative_module2(cfg: dict[str, Any]) -> str:
    pp = module1_embedding_params(cfg)
    qc = cfg["module1"]["qc"]
    vers = environment_versions()
    return (
        "HuBMAP lung snRNA-seq nuclei were loaded from Salmon-processed `raw_expr.h5ad` "
        f"(scanpy {vers.get('scanpy')}, anndata {vers.get('anndata')}). "
        "QC used scanpy.pp.calculate_qc_metrics for gene counts, UMI totals, and mitochondrial/"
        f"ribosomal percentages, then filtered to min_genes={qc['min_genes']}, "
        f"min_counts={qc['min_counts']}, max_mito_percent={qc['max_mito_percent']}, "
        f"max_genes={qc['max_genes']}, and genes in >={qc['min_cells_per_gene']} nuclei. "
        "Doublet detection was not applied in the core teaching path. "
        f"Counts were library-size normalized (target_sum={pp['target_sum']}) and log1p-transformed "
        "with raw counts retained in `.layers['counts']`. "
        f"Highly variable genes (Seurat flavor, n={pp['n_top_genes']}) were scaled "
        f"(clip={pp.get('scale_max_value', 10)}) and embedded with PCA "
        f"({pp['n_pcs']} components, arpack). A kNN graph (n_neighbors={pp['n_neighbors']}) "
        f"and UMAP (min_dist={pp.get('umap_min_dist')}, spread={pp.get('umap_spread', 1.0)}) "
        "were computed for visualization."
    )


def methods_narrative_module3(cfg: dict[str, Any]) -> str:
    emb = module3_embedding_params(cfg)
    m2 = cfg["module2"]
    harm = m2.get("harmony", {})
    vers = environment_versions()
    try:
        import harmonypy as hm  # noqa: F401

        harm_ver = getattr(hm, "__version__", "unknown")
    except Exception:  # noqa: BLE001
        harm_ver = "unknown"
    return (
        "QC-preprocessed Module 1 objects were concatenated on versionless Ensembl IDs "
        "using anndata.concat with join='inner'. "
        f"Joint HVGs (n={emb['n_top_genes']}) entered PCA ({emb['n_pcs']} PCs). "
        f"Donor effects were corrected with Harmony (harmonypy {harm_ver}) on X_pca "
        f"using batch key `{harm.get('key')}`, max_iter_harmony={harm.get('max_iter_harmony', 20)}, "
        f"and random_state={m2.get('random_seed', 0)}. "
        f"Neighbors (n_neighbors={emb['n_neighbors']}) and UMAP "
        f"(min_dist={emb['umap_min_dist']}, spread={emb['umap_spread']}) were recomputed "
        "on the Harmony embedding for visualization. "
        f"Cell labels prefer `{m2.get('label_column', 'azimuth_label')}` for biological "
        f"interpretation; Leiden (resolution={emb.get('leiden_resolution', 0.8)}) is always computed "
        "for diagnostic label-comparison panels, not only as a fallback. "
        "GTEx comparison is deferred to Module 4. "
        f"Software: scanpy {vers.get('scanpy')}, anndata {vers.get('anndata')}."
    )


def _effective_exclude_ribosomal_enrichment(cfg: dict[str, Any]) -> bool:
    de_flag = bool(cfg["module2"]["de"].get("exclude_ribosomal_genes", True))
    enr_flag = cfg["module2"]["enrichment"].get("exclude_ribosomal_genes", None)
    if enr_flag is None:
        return de_flag
    return bool(enr_flag)


def methods_table_module4(cfg: dict[str, Any]) -> pd.DataFrame:
    """DE / pathway methods (learner Module 2 second half; name kept for imports)."""
    m4 = cfg["module2"]
    de = m4["de"]
    enr = m4["enrichment"]
    exclude_ribo_enr = _effective_exclude_ribosomal_enrichment(cfg)
    rows = [
        ("Input", m4["input_h5ad"], "Module 2 integrated AnnData"),
        ("Comparison groups", f"groupby={m4['groupby']}", "obs categorical"),
        (
            "Optional fine DE",
            f"azimuth_focus={m4.get('azimuth_focus')}",
            "scanpy.tl.rank_genes_groups",
        ),
        (
            "Marker DE",
            (
                f"method={de.get('method')}, n_genes={de.get('n_genes')}, "
                f"min_cells_per_group={de.get('min_cells_per_group')}, reference=rest"
            ),
            "scanpy.tl.rank_genes_groups (Wilcoxon)",
        ),
        (
            "Expression matrix for DE",
            "log-normalized `.X` (counts retained in `.layers['counts']`)",
            "scanpy / anndata",
        ),
        (
            "Gene symbol resolution",
            "use adata.var hugo_symbol from HuBMAP object (no GTEx read when gtex_context=false)",
            "scripts.nb04_de.load.gene_symbol_series",
        ),
        (
            "Marker DE table filter",
            (
                f"exclude_ribosomal_genes={de.get('exclude_ribosomal_genes', True)} "
                "(RPS/RPL/MRPS/MRPL plus named exceptions UBA52, FAU removed from "
                "reported marker tables)"
            ),
            "pandas filter after rank_genes_groups",
        ),
        (
            "Enrichment gene filter",
            (
                f"padj<={de.get('max_padj')}, logFC>={de.get('min_logfoldchange')}, "
                f"top={enr.get('top_genes_per_group')}, "
                f"exclude_ribosomal_genes={exclude_ribo_enr} "
                f"(config null follows de.exclude_ribosomal_genes={de.get('exclude_ribosomal_genes', True)})"
            ),
            "pandas filters (drop RPS/RPL/MRPS/MRPL and UBA52/FAU)",
        ),
        (
            "Pathway enrichment (ORA)",
            (
                f"mode={enr.get('mode', 'local_gmt')}, gene_sets={enr.get('gene_sets')}, "
                f"background=measured gene symbols, cutoff_padj={enr.get('cutoff_padj')}"
            ),
            "gseapy.enrich against local GMT (Enrichr online opt-in only)",
        ),
        (
            "Preranked GSEA",
            (
                f"enabled={(enr.get('prerank') or {}).get('enabled', True)}, "
                f"permutation_num={(enr.get('prerank') or {}).get('permutation_num', 100)}"
            ),
            "gseapy.prerank on Wilcoxon scores (NES)",
        ),
        (
            "Known marker validation",
            "curated lung marker panel via HUGO when present, else Ensembl ID in var_names",
            "dotplot",
        ),
        (
            "GTEx context",
            f"enabled={m4.get('gtex_context')} (Module 2 default false; Module 4 only)",
            "not used in Module 2",
        ),
        ("Random seed", str(m4.get("random_seed", 0)), "config"),
    ]
    return pd.DataFrame(rows, columns=["step", "parameters", "package_or_function"])


def methods_narrative_module4(cfg: dict[str, Any]) -> str:
    m4 = cfg["module2"]
    de = m4["de"]
    enr = m4["enrichment"]
    vers = environment_versions()
    exclude_ribo_enr = _effective_exclude_ribosomal_enrichment(cfg)
    try:
        import gseapy as gp

        gseapy_ver = getattr(gp, "__version__", "unknown")
    except Exception:  # noqa: BLE001
        gseapy_ver = "unavailable"
    mode = str(enr.get("mode", "local_gmt"))
    return (
        "Marker genes were identified on the Module 2 integrated HuBMAP lung object "
        f"using scanpy.tl.rank_genes_groups ({de.get('method', 'wilcoxon')}) grouped by "
        f"`{m4['groupby']}` with reference='rest' on log-normalized expression "
        f"(scanpy {vers.get('scanpy')}). "
        "Gene symbols use `adata.var['hugo_symbol']` from the HuBMAP object "
        "(no GTEx read when `gtex_context=false`). "
        + (
            "Ribosomal protein genes (RPS/RPL/MRPS/MRPL, plus UBA52 and FAU which "
            "encode ribosomal proteins under non-RP names) were removed from reported "
            "marker tables so top-gene lists are not dominated by translation housekeeping. "
            if de.get("exclude_ribosomal_genes", True)
            else ""
        )
        + f"Genes with adjusted p<={de.get('max_padj')} and log-fold change>={de.get('min_logfoldchange')} "
        f"were scored with over-representation enrichment (gseapy {gseapy_ver}, mode={mode}) "
        f"against local GMT libraries {enr.get('gene_sets')} using the object's measured gene "
        "symbols as background"
        + (
            ", again excluding ribosomal protein genes from enrichment inputs"
            if exclude_ribo_enr
            else ""
        )
        + ". Preranked GSEA (NES) uses Wilcoxon scores against the same GMTs when enabled. "
        "Cross-resource GTEx context is deferred to Module 4. Enrichment results are "
        "hypothesis-generating; this pipeline does not assert biological consistency."
    )


def methods_table_module5(cfg: dict[str, Any]) -> pd.DataFrame:
    """Multiome methods (learner Module 3; function name kept for imports)."""
    m5 = cfg["module3"]
    files = m5.get("files", {})
    analysis = m5.get("analysis", {})
    bridge_metric = analysis.get("bridge_metric", "feature_std_on_mofa_inputs")
    load_h5mu = bool((m5.get("params") or {}).get("load_h5mu", False))
    rows = [
        ("Donor", f"{m5['donor_id']} ({m5['donor_label']})", "true multiome teaching donor"),
        ("Dataset", m5.get("dataset_id", ""), m5.get("assay", "SNARE-seq2 [Salmon + ArchR + Muon]")),
        (
            "Portal URL",
            m5.get("portal_url", "https://portal.hubmapconsortium.org"),
            "HuBMAP native dataset page",
        ),
        (
            "Default path",
            "interpretation_only from multiome_mofa.hdf5"
            if not load_h5mu
            else "full MuData + MOFA path",
            "Phase 0/3 download decision",
        ),
        (
            "MuData input",
            f"{m5.get('multiome_dir')}/{files.get('secondary_h5mu', 'secondary_analysis.h5mu')}",
            "optional unless load_h5mu=true",
        ),
        (
            "MOFA input",
            f"{m5.get('multiome_dir')}/{files.get('mofa_hdf5', 'multiome_mofa.hdf5')}",
            "joint factors + rna/atac_cbg feature matrices",
        ),
        (
            "Gene bridge",
            f"{bridge_metric}: log1p(RNA feature SD) vs log1p(ATAC feature SD) on MOFA inputs (HUGO join)",
            "MOFA inputs are centered; SD captures shared variability (not mean expression)",
        ),
        (
            "Focus markers",
            (
                f"{analysis.get('focus_genes')}; "
                f"rna_layer={analysis.get('rna_layer', 'spliced_unspliced_sum')} when MuData loaded"
            ),
            "global + within-cell-type Pearson; flag low-n and rna_mean==0",
        ),
        (
            "MOFA honesty checks",
            "recompute R^2 from Z/W/Y; max non-dominant R^2; truncate dead trailing factors",
            "scripts.nb05_multimodal.load.recompute_mofa_r2_total",
        ),
        (
            "Markers not in both spaces",
            str(analysis.get("focus_genes_not_in_both_modalities", [])),
            "documented omissions",
        ),
        ("Random seed", str(analysis.get("random_seed", 0)), "config"),
    ]
    return pd.DataFrame(rows, columns=["step", "parameters", "package_or_function"])


def methods_narrative_module5(cfg: dict[str, Any]) -> str:
    m5 = cfg["module3"]
    vers = environment_versions()
    bridge_metric = m5.get("analysis", {}).get("bridge_metric", "feature_std_on_mofa_inputs")
    load_h5mu = bool((m5.get("params") or {}).get("load_h5mu", False))
    return (
        f"Module 3 examined true multiome context for HuBMAP donor {m5['donor_id']} "
        f"using SNARE-seq2 Salmon+ArchR+Muon dataset {m5.get('dataset_id')}. "
        + (
            "The default teaching path interprets the fitted MOFA model (`multiome_mofa.hdf5`) "
            "without requiring the full MuData download. "
            if not load_h5mu
            else "Analysis used MuData embeddings plus the fitted MOFA model. "
        )
        + "RNA and ATAC modalities occupy the same barcode space when MuData is loaded. "
        "Joint structure is summarized with MOFA factors (rna + atac_cbg views). "
        f"A genome-wide gene bridge uses {bridge_metric} on MOFA inputs (joined by HUGO symbol). "
        "Claims about shared latent structure must match the variance table (re-measured in Phase 3). "
        f"Software: scanpy {vers.get('scanpy')}, anndata {vers.get('anndata')}."
    )
