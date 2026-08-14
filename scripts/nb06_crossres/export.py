"""Export Module 4 cross-ecosystem tables, figures, and methods narrative."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from scripts.common.io import environment_versions
from scripts.common.paths import ensure_output_dirs, portable_path, resolve


def write_methods_narrative(
    cfg: dict[str, Any],
    path: Path,
    *,
    sample_counts: pd.DataFrame,
    ortholog_loss: pd.DataFrame,
    nes_stats: dict[str, Any],
    extras: dict[str, Any] | None = None,
    nes_pair_summary: pd.DataFrame | None = None,
) -> Path:
    m4 = cfg["module4"]
    ometa = m4.get("orthologs_meta") or {}
    extras = extras or {}
    ercc = extras.get("ercc_calibration") or {}
    lines = [
        "# Module 4 -- Cross-ecosystem comparison (methods + results draft)",
        "",
        f"_Auto-drafted at {datetime.now(timezone.utc).isoformat()}. Edit with learner conclusions._",
        "",
        "## Framing",
        "",
        "Module 4 is **cross-ecosystem**, not a merge. Four resources, four within-resource "
        "contrasts, plus a HuBMAP<->GTEx level bridge. Counts are **never combined** across species "
        "or studies. What travels is the direction and relative magnitude of a change "
        "(logFC / NES), never the native value.",
        "",
        "GEO (GSE150910) is NIH but **not** CFDE. OSDR is NASA. GTEx and HuBMAP are CFDE core.",
        "",
        "Methods assets can succeed while biology pairs stay near-null or uninterpretable: "
        "ERCC calibration and the GeneLab cross-check validate the pipeline; global NES is "
        "supporting only; shared-mechanism concordance is the headline.",
        "",
        "## Resources and contrasts",
        "",
        "| Ecosystem | Resource | Contrast | Native unit |",
        "|---|---|---|---|",
        "| NIH (not CFDE) | GSE150910 | H_DISEASE: IPF vs control | bulk counts |",
        "| CFDE | GTEx v10 lung | H_AGING_GTEX: older vs younger (gene_reads) | bulk counts |",
        "| NASA | OSD-248 | M_FLIGHT: ISS-T flight vs GC | bulk counts (mouse; RSEM rounded) |",
        "| CFDE | HuBMAP | H_AGING_HUBMAP: older vs younger donors (2 vs 2 illustration) | summed counts-layer (rounded) |",
        "| CFDE | HuBMAP / GTEx | level bridge only | linear CP10K vs TPM |",
        "",
        "## Sample counts (from data)",
        "",
    ]
    for _, row in sample_counts.iterrows():
        lines.append(
            f"- **{row['resource']} / {row['contrast_id']} / {row['arm']}**: n={int(row['n'])} "
            f"({row.get('notes', '')})"
        )
    lines.extend(
        [
            "",
            "## Ortholog loss (M_FLIGHT only)",
            "",
            f"- Source: {ometa.get('source')} ({ometa.get('download_date')})",
            "",
        ]
    )
    for _, row in ortholog_loss.iterrows():
        lines.append(f"- {row['metric']}: {int(row['n'])}")

    lines.extend(["", "## NES pairs (Hallmark + Reactome)", ""])
    if nes_pair_summary is not None and not nes_pair_summary.empty:
        lines.append("| Pair | Frame | n | Spearman rho | Sign agree |")
        lines.append("|---|---|---:|---:|---:|")
        for _, row in nes_pair_summary.iterrows():
            try:
                rho = f"{float(row['spearman_rho']):.3f}"
            except (TypeError, ValueError):
                rho = str(row.get("spearman_rho"))
            try:
                agree = f"{100 * float(row['sign_agreement']):.1f}%"
            except (TypeError, ValueError):
                agree = str(row.get("sign_agreement"))
            lines.append(
                f"| {row['contrast_x']} vs {row['contrast_y']} | {row.get('frame')} | "
                f"{int(row['n_shared_pathways'])} | {rho} | {agree} |"
            )
        lines.append("")
        lines.append(
            "`H_DISEASE` vs `H_AGING_GTEX` is the **best-case** within-species pair "
            "(same tissue, both powered) -- not a guaranteed positive control. "
            "IPF is age-associated; it is not aging."
        )
        lines.append(
            "`H_AGING_HUBMAP` pairs only with `H_AGING_GTEX` (method question). "
            "State age-range and power differences before reading that result."
        )
    else:
        r_nes = nes_stats.get("spearman_nes")
        try:
            r_nes_fmt = f"{float(r_nes):.3f}"
        except (TypeError, ValueError):
            r_nes_fmt = str(r_nes)
        lines.append(
            f"- Headline pair stats available: rho~{r_nes_fmt}, "
            f"n={nes_stats.get('n_shared_pathways')}"
        )

    lines.extend(
        [
            "",
            "- ``spearman_p`` is not headline evidence: overlapping gene sets violate independence.",
            "- Gene-set sensitivity table (Hallmark / Reactome / GO BP slices) ships as committed "
            "`module4_nes_sensitivity.tsv`. Default `module4.sensitivity.enabled: false` loads it; "
            "`true` re-preranks (high RAM).",
            "",
            "## Pipeline calibration",
            "",
        ]
    )
    if ercc:
        lines.append(
            f"- ERCC ExFold (flight=Mix1, GC=Mix2): n_obs={ercc.get('n_ercc_observed', 'NA')}/92; "
            f"OLS slope={ercc.get('ols_slope', float('nan'))}, "
            f"intercept={ercc.get('ols_intercept', float('nan'))}, "
            f"R^2={ercc.get('ols_r2', float('nan'))} "
            f"(Spearman secondary rho~{ercc.get('spearman_obs_vs_exp', float('nan'))})."
        )
        lines.append(
            "- Mix is aligned with arm, so spike-ins calibrate and cannot normalize this contrast."
        )
        lines.append(
            "- `ERCC-` rows remain in the DE table for ExFold calibration; size factors are "
            "estimated on endogenous genes only (`exclude_ercc_from_size_factors`). "
            "ERCC read fraction is recorded by sample; excluding spike-ins from size factors "
            "does not move pathway passer `p_up` off zero."
        )
    if extras.get("genelab_logfc_spearman") is not None:
        lines.append(
            f"- GeneLab cross-check: Spearman rho~{extras.get('genelab_logfc_spearman')} "
            f"across {extras.get('n_genelab_joined', 'NA')} shared gene symbols "
            "(spike-in IDs are not in the GeneLab gene table)."
        )
    lines.extend(
        [
            "",
            "## Confounds (one clause each)",
            "",
            "- GTEx is postmortem: age is entangled with cause of death and ischemic time.",
            "- OSD-248 is female C57BL/6NTac at 36 weeks: nothing here estimates an age effect.",
            "- GSE150910 controls are not matched to GTEx donors.",
            "- ISS-T flight carcasses were transported from orbit; ground carcasses were not "
            "(both held at -80; dissection months after euthanasia). Flight euthanasia dates "
            "06-09 Feb vs GC 03-06 Feb; dissection dates identical; **no times of day recorded** "
            "(circadian genes often dominate non-ERCC ranks).",
            "- H_DISEASE uses `~ condition` only; the published analysis adjusted for clinical covariates.",
            "- OSD-248 ERCC spike-ins often dominate padj ranks because flight and GC used different "
            "ExFold mixes; volcanoes exclude ERCC by default.",
            "- H_AGING_HUBMAP is 2 vs 2 illustration only; HuBMAP older ages sit inside GTEx middle brackets.",
            "- HuBMAP integrated `counts` layer is not integer UMIs; donor sums are rounded before pydeseq2.",
            "- No CFDE program supplies fibrotic human lung; that is why Module 4 goes to GEO (NIH, not CFDE).",
            "- Reference levels are set explicitly (`Treatment(reference)` + DeseqStats contrast). "
            "Alphabetical factor order would have inverted three of four contrasts "
            "(`IPF` before `control`; `older` before `younger`).",
            "",
            "## What this module teaches",
            "",
            "1. Levels never travel (four resources, four unit systems).",
            "2. Global NES correlations show these contrasts are not broadly equivalent (supporting line).",
            "3. Shared-mechanism concordance names pathways that pass thresholds in the same direction, "
            "against chance; that is the headline.",
            "4. Validate method before interpreting biology (ERCC, GeneLab; explicit DE polarity).",
            "5. An underpowered contrast looks like nothing; recognizing that is a skill.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join([ln for ln in lines if ln is not None]))
    return path


def save_module4_outputs(
    cfg: dict[str, Any],
    *,
    sample_counts: pd.DataFrame,
    native_levels: pd.DataFrame,
    contrast_h: pd.DataFrame,
    contrast_m: pd.DataFrame,
    contrast_m_sens: pd.DataFrame | None,
    contrast_aging_gtex: pd.DataFrame | None,
    contrast_aging_hubmap: pd.DataFrame | None,
    ortholog_loss: pd.DataFrame,
    nes_by_id: dict[str, pd.DataFrame],
    nes_pair_summary: pd.DataFrame,
    nes_pair_details: dict[str, tuple[pd.DataFrame, dict[str, Any]]],
    bridge_df: pd.DataFrame,
    figure_paths: dict[str, Path],
    extras: dict[str, Any] | None = None,
    nes_sensitivity: pd.DataFrame | None = None,
    ercc_calibration: pd.DataFrame | None = None,
    genelab_subset: pd.DataFrame | None = None,
    gtex_bracket_inventory: pd.DataFrame | None = None,
    # Legacy kwargs accepted for transitional callers
    contrast_aging: pd.DataFrame | None = None,
    nes_human: pd.DataFrame | None = None,
    nes_mouse: pd.DataFrame | None = None,
    nes_joined: pd.DataFrame | None = None,
    nes_stats: dict[str, Any] | None = None,
) -> dict[str, Path]:
    ensure_output_dirs(cfg)
    tables = resolve(cfg, "outputs_tables")
    reports = resolve(cfg, "outputs_reports")

    if contrast_aging_gtex is None and contrast_aging is not None:
        contrast_aging_gtex = contrast_aging

    headline = list(
        (cfg.get("module4", {}).get("params") or {}).get("nes_headline_pair")
        or ["H_DISEASE", "H_AGING_GTEX"]
    )
    headline_key = f"{headline[0]}__{headline[1]}"
    if nes_stats is None:
        nes_stats = (nes_pair_details.get(headline_key) or (None, {}))[1] or {}
    if nes_joined is None and headline_key in nes_pair_details:
        nes_joined = nes_pair_details[headline_key][0]

    paths: dict[str, Path] = {
        "sample_counts": tables / "module4_sample_counts.tsv",
        "native_levels": tables / "module4_native_levels_panel.tsv",
        "contrast_h_disease": tables / "module4_contrast_H_DISEASE.tsv",
        "contrast_m_flight": tables / "module4_contrast_M_FLIGHT.tsv",
        "ortholog_loss": tables / "module4_ortholog_loss.tsv",
        "nes_pair_summary": tables / "module4_nes_pair_summary.tsv",
        "nes_sensitivity": tables / "module4_nes_sensitivity.tsv",
        "ercc_calibration": tables / "module4_ercc_exfold_calibration.tsv",
        "bridge": tables / "module4_hubmap_gtex_bridge.tsv",
        "methods": reports / "module4_methods_narrative.md",
        "run_params": reports / "module4_run_params.json",
    }

    sample_counts.to_csv(paths["sample_counts"], sep="\t", index=False)
    native_levels.to_csv(paths["native_levels"], sep="\t", index=False)
    contrast_h.to_csv(paths["contrast_h_disease"], sep="\t", index=False)
    contrast_m.to_csv(paths["contrast_m_flight"], sep="\t", index=False)
    ortholog_loss.to_csv(paths["ortholog_loss"], sep="\t", index=False)
    nes_pair_summary.to_csv(paths["nes_pair_summary"], sep="\t", index=False)
    bridge_df.to_csv(paths["bridge"], sep="\t", index=False)

    if contrast_m_sens is not None and not contrast_m_sens.empty:
        paths["contrast_m_flight_qa_filtered"] = (
            tables / "module4_contrast_M_FLIGHT_qa_filtered.tsv"
        )
        contrast_m_sens.to_csv(paths["contrast_m_flight_qa_filtered"], sep="\t", index=False)
    if contrast_aging_gtex is not None and not contrast_aging_gtex.empty:
        paths["contrast_h_aging_gtex"] = tables / "module4_contrast_H_AGING_GTEX.tsv"
        contrast_aging_gtex.to_csv(paths["contrast_h_aging_gtex"], sep="\t", index=False)
    if contrast_aging_hubmap is not None and not contrast_aging_hubmap.empty:
        paths["contrast_h_aging_hubmap"] = tables / "module4_contrast_H_AGING_HUBMAP.tsv"
        contrast_aging_hubmap.to_csv(paths["contrast_h_aging_hubmap"], sep="\t", index=False)
    if gtex_bracket_inventory is not None and not gtex_bracket_inventory.empty:
        paths["gtex_age_brackets"] = tables / "module4_gtex_age_bracket_inventory.tsv"
        gtex_bracket_inventory.to_csv(paths["gtex_age_brackets"], sep="\t", index=False)

    for cid, nes_df in nes_by_id.items():
        p = tables / f"module4_nes_{cid}.tsv"
        paths[f"nes_{cid}"] = p
        nes_df.to_csv(p, sep="\t", index=False)

    for key, (joined, _stats) in nes_pair_details.items():
        p = tables / f"module4_nes_comparison_{key}.tsv"
        paths[f"nes_comparison_{key}"] = p
        joined.to_csv(p, sep="\t", index=False)

    # Keep a convenience copy of the headline pair under the old filename
    if nes_joined is not None and not nes_joined.empty:
        paths["nes_comparison"] = tables / "module4_nes_comparison.tsv"
        nes_joined.to_csv(paths["nes_comparison"], sep="\t", index=False)

    if nes_sensitivity is not None and not nes_sensitivity.empty:
        nes_sensitivity.to_csv(paths["nes_sensitivity"], sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(paths["nes_sensitivity"], sep="\t", index=False)
    if ercc_calibration is not None and not ercc_calibration.empty:
        ercc_calibration.to_csv(paths["ercc_calibration"], sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(paths["ercc_calibration"], sep="\t", index=False)
    if genelab_subset is not None and not genelab_subset.empty:
        subset_path = tables / "module4_osd248_genelab_isst_de_subset.tsv"
        paths["genelab_isst_de_subset"] = subset_path
        if not subset_path.exists():
            genelab_subset.to_csv(subset_path, sep="\t", index=False)

    write_methods_narrative(
        cfg,
        paths["methods"],
        sample_counts=sample_counts,
        ortholog_loss=ortholog_loss,
        nes_stats=nes_stats,
        extras=extras,
        nes_pair_summary=nes_pair_summary,
    )

    params = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "module4": cfg["module4"],
        "nes_stats_headline": nes_stats,
        "nes_pair_summary": nes_pair_summary.to_dict(orient="records"),
        "environment": environment_versions(),
        "figures": {k: portable_path(cfg, v) for k, v in figure_paths.items()},
        "extras": extras or {},
        "_compute_seconds": (extras or {}).get("_compute_seconds"),
        "peak_rss_mb": (extras or {}).get("peak_rss_mb"),
        "honesty": {
            "counts_merged_across_studies": False,
            "ortholog_loss_reported": True,
            "geo_is_cfde": False,
            "nes_headline_gene_sets": (cfg.get("module4", {}).get("enrichment") or {}).get(
                "comparison_gene_sets"
            ),
            "nes_headline_pair": headline,
            "osd248_counts_pipeline": "RSEM_Unnormalized_Counts_rRNArm",
            "hubmap_donor_counts_rounded": True,
            "h_aging_hubmap_illustration_only": True,
        },
    }
    # Prefer repo-relative paths inside extras when present.
    ex = params["extras"]
    for key in ("sensitivity_table", "genelab_isst_de_subset", "table", "path"):
        if key in ex and ex[key]:
            ex[key] = portable_path(cfg, ex[key])
    meta = ex.get("nes_sensitivity_meta")
    if isinstance(meta, dict):
        for key, val in list(meta.items()):
            if isinstance(val, (str, Path)) and ("/" in str(val) or str(val).endswith(".tsv")):
                meta[key] = portable_path(cfg, val)
    try:
        import pydeseq2

        params["environment"]["pydeseq2"] = getattr(pydeseq2, "__version__", "unknown")
    except Exception:  # noqa: BLE001
        params["environment"]["pydeseq2"] = "unavailable"

    paths["run_params"].write_text(json.dumps(params, indent=2, default=str))
    for k, v in figure_paths.items():
        paths[f"figure_{k}"] = Path(v)
    return paths


def append_mechanism_section(
    methods_path: Path,
    *,
    null_df: pd.DataFrame,
    summary: pd.DataFrame,
    mechanisms: pd.DataFrame,
) -> None:
    """Append mechanism-concordance close to an existing methods narrative."""
    lines = [
        "",
        "## Shared-mechanism concordance (headline)",
        "",
        "Global NES Spearman is supporting only. The headline is pathways that pass "
        "FDR q<0.25 and |NES|>=1.5 in each contrast independently, with agreeing sign, "
        "compared to analytic and permutation nulls. All pathways ship as tier 3 "
        "(concordance only). Leading-edge Jaccard on this run spans roughly 0.21 to 0.48 "
        "with no natural tier 1/2 cut; do not invent one.",
        "",
        "Suspect ribosome / translation / OxPhos terms are flagged, not excluded.",
        "",
    ]
    if null_df is not None and not null_df.empty:
        lines.append("| Set | Obs | Unflagged | Perm mean | Obs/perm | perm p |")
        lines.append("|---|---:|---:|---:|---:|---:|")
        for _, r in null_df.iterrows():
            lines.append(
                f"| {r['set_id']} | {int(r['observed_concordant'])} | "
                f"{int(r['observed_concordant_unflagged'])} | "
                f"{float(r['perm_mean']):.2f} | "
                f"{float(r['perm_ratio_obs_over_mean']) if r['perm_ratio_obs_over_mean']==r['perm_ratio_obs_over_mean'] else float('nan'):.2f} | "
                f"{float(r['perm_p_ge_obs']):.3f} |"
            )
        lines.append("")
    if summary is not None and not summary.empty:
        for _, r in summary.iterrows():
            n_u = int(r["n_mechanisms_unflagged"])
            n_t = int(r["n_mechanisms_total"])
            lines.append(
                f"- **{r['set_id']}**: {n_u} mechanisms unflagged "
                f"({n_t} including flagged); "
                f"{int(r['n_concordant_pathways'])} concordant pathways "
                f"(perm p={float(r['perm_p_ge_obs']):.3f})."
            )
        lines.append("")
    if mechanisms is not None and not mechanisms.empty:
        flagged = mechanisms.loc[mechanisms["suspect_flag_any"], "representative_term"]
        if len(flagged):
            lines.append(
                "Flagged representatives: " + "; ".join(flagged.astype(str).tolist()) + "."
            )
            lines.append("")
    text = methods_path.read_text() if methods_path.exists() else ""
    methods_path.write_text(text.rstrip() + "\n" + "\n".join(lines) + "\n")
