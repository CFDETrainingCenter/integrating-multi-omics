"""Export Module 2 DE / pathway tables, figures, methods, and interpretation report."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from scripts.common.io import environment_versions
from scripts.common.methods import methods_narrative_module4, methods_table_module4
from scripts.common.paths import ensure_output_dirs, resolve, portable_path
from scripts.nb04_de.enrich import enrichment_failed_groups


def _input_gene_counts(enrichment_df: pd.DataFrame) -> dict[str, int]:
    if enrichment_df is None or enrichment_df.empty:
        return {}
    if "group" not in enrichment_df.columns or "n_input_genes" not in enrichment_df.columns:
        return {}
    out: dict[str, int] = {}
    for g, sub in enrichment_df.groupby("group", observed=True):
        out[str(g)] = int(sub["n_input_genes"].iloc[0])
    return out


def write_interpretation_markdown(
    cfg: dict[str, Any],
    marker_df: pd.DataFrame,
    enrichment_df: pd.DataFrame,
    path: Path,
    *,
    prerank_df: pd.DataFrame | None = None,
) -> Path:
    groupby = cfg["module2"]["groupby"]
    failed = enrichment_failed_groups(enrichment_df)
    n_inputs = _input_gene_counts(enrichment_df)
    enr = cfg["module2"]["enrichment"]
    mode = str(enr.get("mode", "local_gmt"))

    lines = [
        "# Module 2 -- Biological interpretation summary (DE + pathways)",
        "",
        f"_Auto-drafted from run at {datetime.now(timezone.utc).isoformat()}. "
        "Edit with learner conclusions -- this file does not assert biological consistency._",
        "",
        "## Comparison design",
        "",
        f"- Grouping variable: `{groupby}`",
        f"- DE method: `{cfg['module2']['de'].get('method', 'wilcoxon')}` (scanpy.tl.rank_genes_groups)",
        f"- Enrichment mode: `{mode}` against `{enr.get('gene_sets')}`",
        (
            f"- Ribosomal filter: "
            f"`de.exclude_ribosomal_genes={cfg['module2']['de'].get('exclude_ribosomal_genes', True)}` "
            "removes RPS/RPL/MRPS/MRPL plus named exceptions UBA52 and FAU from reported "
            "marker tables, ORA inputs, and Wilcoxon "
            "prerank ranked lists (same effective enrichment flag)"
        ),
        "- GTEx context: off in Module 2 (cross-resource comparison is Module 4)",
        "",
        "## Top markers (per group)",
        "",
    ]
    if marker_df.empty:
        lines.append("_No markers exported._")
    else:
        top = marker_df.copy()
        if "pvals_adj" in top.columns:
            top = top.sort_values(["group", "pvals_adj"]).groupby("group", observed=True).head(15)
        top = top[top["gene_symbol"].notna()]
        top = top[~top["gene_symbol"].astype(str).isin(["", "nan", "None"])]
        top = top.groupby("group", observed=True).head(5)
        for group, sub in top.groupby("group", observed=True):
            genes = ", ".join(sub["gene_symbol"].astype(str).tolist())
            lines.append(f"- **{group}**: {genes}")
        lines.append("")

    lines.extend(["## Over-representation enrichment (per group)", ""])
    if enrichment_df is None or enrichment_df.empty:
        lines.append("_No enrichment results (disabled, failed, or insufficient genes)._")
    else:
        if n_inputs:
            lines.append("### Input list sizes (`n_input_genes`)")
            lines.append("")
            for g, n in sorted(n_inputs.items()):
                lines.append(f"- **{g}**: {n} genes submitted")
            lines.append("")
            lines.append(
                "Small lists can produce extreme or unexpected terms; treat ORA as "
                "hypothesis-generating and check gene membership before claiming biology."
            )
            lines.append("")

        if failed:
            lines.append("### Failed / skipped enrichment groups")
            lines.append("")
            for g in failed:
                err = ""
                sub = enrichment_df[
                    (enrichment_df["group"].astype(str) == g)
                    & (
                        enrichment_df["Term"].astype(str).str.startswith("ENRICHMENT_FAILED")
                        | enrichment_df["Term"].astype(str).str.startswith("ENRICHMENT_SKIPPED")
                    )
                ]
                if not sub.empty and "error" in sub.columns:
                    err = f" -- {sub.iloc[0]['error']}"
                elif not sub.empty:
                    err = f" -- {sub.iloc[0]['Term']}"
                lines.append(f"- **{g}**: enrichment incomplete{err}")
            lines.append("")

        df = enrichment_df.copy()
        if "Term" in df.columns:
            df = df[
                ~df["Term"].astype(str).str.startswith("ENRICHMENT_FAILED")
                & ~df["Term"].astype(str).str.startswith("ENRICHMENT_SKIPPED")
            ]
        padj = "Adjusted P-value" if "Adjusted P-value" in df.columns else None
        if padj is None and "FDR q-val" in df.columns:
            padj = "FDR q-val"
        if padj:
            df = df[df[padj].astype(float) <= 0.05]
            df = df.sort_values(["group", padj]).groupby("group", observed=True).head(5)
        if df.empty:
            lines.append("_No significant enriched terms after filtering._")
        else:
            for group, sub in df.groupby("group", observed=True):
                terms = "; ".join(sub["Term"].astype(str).tolist())
                lines.append(f"- **{group}**: {terms}")
        lines.append("")

    lines.extend(["## Preranked GSEA (NES)", ""])
    if prerank_df is None or prerank_df.empty:
        lines.append("_No prerank results in this run (disabled or empty)._")
    else:
        lines.append(
            "NES columns come from `gseapy.prerank` on Wilcoxon scores against local GMTs. "
            "Terms are listed for inspection; biological plausibility is the learner's task. "
            "Even after ribosomal-protein filtering, translation-adjacent Reactome/GO terms can "
            "still rank high for small secretory compartments -- treat marker genes as the "
            "identity check before claiming a pathway story."
        )
        lines.append("")
        # Prefer FDR / NES columns if present
        show = prerank_df.copy()
        fdr_col = next((c for c in ("FDR q-val", "FDR", "fdr") if c in show.columns), None)
        nes_col = next((c for c in ("NES", "nes") if c in show.columns), None)
        term_col = next((c for c in ("Term", "term", "Name") if c in show.columns), None)
        if fdr_col and term_col:
            show = show.sort_values(fdr_col).head(15)
            for _, row in show.iterrows():
                nes = f", NES={row[nes_col]:.3f}" if nes_col and pd.notna(row.get(nes_col)) else ""
                grp = f" [{row['group']}]" if "group" in show.columns else ""
                lines.append(f"- {row[term_col]}{grp}{nes}")
        else:
            lines.append(f"_Prerank table has {len(prerank_df)} rows (see TSV)._")
        lines.append("")

    lines.extend(
        [
            "## Interpretation prompts (learner task)",
            "",
            "1. Which cell population or cluster was analyzed?",
            "2. What are the top marker genes, and do they match known lung identities?",
            "3. What pathways are enriched, and are they biologically coherent for that population?",
            "4. What biological interpretation is supported by markers + pathways together?",
            "5. What limitations remain (label transfer uncertainty, enrichment bias, small input lists)?",
            "",
            "## Three-sentence results draft (edit me -- no auto biology claim)",
            "",
            (
                f"Using the Module 2 integrated HuBMAP lung snRNA-seq object, we identified marker genes "
                f"for `{groupby}` groups with Wilcoxon rank-sum testing (scanpy) and scored pathways "
                f"with gseapy (`{mode}` ORA"
                + (" + prerank GSEA" if prerank_df is not None and not prerank_df.empty else "")
                + ") against local Hallmark / Reactome / GO BP gene sets. "
                "Enriched terms are listed above with per-group `n_input_genes`; they have **not** been "
                "checked for biological plausibility in this auto-draft -- that check is the learner's task. "
                "Because labels are reference-mapped and enrichment is hypothesis-generating, "
                "follow-up validation should use independent cohorts or orthogonal assays."
            ),
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))
    return path


def save_module2_de_outputs(
    cfg: dict[str, Any],
    marker_df: pd.DataFrame,
    enrichment_df: pd.DataFrame,
    known_markers_df: pd.DataFrame,
    figure_paths: dict[str, Path],
    extras: dict[str, Any] | None = None,
    *,
    prerank_df: pd.DataFrame | None = None,
    dropped_groups: pd.DataFrame | None = None,
) -> dict[str, Path]:
    ensure_output_dirs(cfg)
    tables = resolve(cfg, "outputs_tables")
    reports = resolve(cfg, "outputs_reports")

    paths: dict[str, Path] = {
        "markers": tables / "module2_marker_genes.tsv",
        "enrichment": tables / "module2_pathway_enrichment.tsv",
        "prerank": tables / "module2_prerank_gsea.tsv",
        "known_markers": tables / "module2_known_marker_presence.tsv",
        "methods_parameters": tables / "module2_de_methods_parameters.tsv",
        "interpretation": reports / "module2_biological_interpretation.md",
        "run_params": reports / "module2_de_run_params.json",
    }
    if dropped_groups is not None:
        paths["dropped_groups"] = tables / "module2_groups_below_min_cells.tsv"
        dropped_groups.to_csv(paths["dropped_groups"], sep="\t", index=False)

    marker_df.to_csv(paths["markers"], sep="\t", index=False)
    # gseapy returns Genes in nondeterministic order; sort so re-runs diff cleanly.
    if "Genes" in enrichment_df.columns:
        enrichment_df = enrichment_df.copy()
        enrichment_df["Genes"] = (
            enrichment_df["Genes"]
            .fillna("")
            .map(lambda s: ";".join(sorted(str(s).split(";"))) if s else s)
        )
    enrichment_df.to_csv(paths["enrichment"], sep="\t", index=False)
    (prerank_df if prerank_df is not None else pd.DataFrame()).to_csv(
        paths["prerank"], sep="\t", index=False
    )
    known_markers_df.to_csv(paths["known_markers"], sep="\t", index=False)
    methods_table_module4(cfg).to_csv(paths["methods_parameters"], sep="\t", index=False)
    write_interpretation_markdown(
        cfg, marker_df, enrichment_df, paths["interpretation"], prerank_df=prerank_df
    )

    params = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "module2_de": cfg["module2"],
        "n_markers": int(marker_df.shape[0]),
        "n_enrichment_rows": int(enrichment_df.shape[0]) if enrichment_df is not None else 0,
        "n_prerank_rows": int(prerank_df.shape[0]) if prerank_df is not None else 0,
        "enrichment_failed_groups": enrichment_failed_groups(enrichment_df),
        "n_input_genes_by_group": _input_gene_counts(enrichment_df),
        "methods_narrative": methods_narrative_module4(cfg),
        "methods_table": methods_table_module4(cfg).to_dict(orient="records"),
        "environment": environment_versions(),
        "figures": {k: portable_path(cfg, v) for k, v in figure_paths.items()},
        "extras": extras or {},
        "_compute_seconds": (extras or {}).get("_compute_seconds"),
        "peak_rss_mb": (extras or {}).get("peak_rss_mb"),
    }
    try:
        import gseapy as gp

        params["environment"]["gseapy"] = getattr(gp, "__version__", "unknown")
    except Exception:  # noqa: BLE001
        params["environment"]["gseapy"] = "unavailable"

    paths["run_params"].write_text(json.dumps(params, indent=2, default=str))
    for k, v in figure_paths.items():
        paths[f"figure_{k}"] = Path(v)
    return paths


def save_module4_outputs(*args, **kwargs):
    raise RuntimeError(
        "save_module4_outputs was renamed to save_module2_de_outputs "
        "(learner Module 2 DE half; no GTEx context table)."
    )
