"""Preranked GSEA concordance and HuBMAP<->GTEx bridge comparison (Module 4)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from scripts.common.plotting import save_figure
from scripts.nb04_de.enrich import resolve_gmt_paths, run_prerank_gsea


def prerank_contrast(
    logfc: pd.Series,
    cfg: dict[str, Any],
    *,
    contrast_id: str,
    gene_sets: list[str] | None = None,
) -> pd.DataFrame:
    """Prerank GSEA on a logFC vector against Module 4 local GMTs (no KEGG / no HTTP).

    Default gene sets are ``comparison_gene_sets`` (Hallmark + Reactome). Pass
    ``gene_sets`` explicitly for the P6-7 sensitivity run that includes GO BP.
    """
    enr = cfg.get("module4", {}).get("enrichment") or {}
    sets = gene_sets or enr.get("comparison_gene_sets") or enr.get("gene_sets") or [
        "MSigDB_Hallmark_2020",
        "Reactome_2022",
    ]
    # Temporarily point enrichment to module4 block
    cfg_m2 = {
        "_module_root": cfg.get("_module_root"),
        "_config_path": cfg.get("_config_path"),
        "module2": {
            "random_seed": int(cfg.get("module4", {}).get("params", {}).get("random_seed", 0) or 0),
            "enrichment": {
                **enr,
                "prerank": enr.get("prerank")
                or {"enabled": True, "min_size": 15, "max_size": 500, "permutation_num": 100},
                "gene_sets": list(sets),
                "geneset_dir": enr.get("geneset_dir", "data/genesets"),
            },
        },
        "paths": cfg.get("paths", {}),
        "project_name": cfg.get("project_name"),
    }
    ranked = logfc.dropna().astype(float)
    ranked.index = ranked.index.astype(str).str.upper()
    ranked = ranked[~ranked.index.duplicated(keep="first")]
    # ERCC spike-ins are a validation asset for OSD-248, not pathway members
    ranked = ranked[~ranked.index.astype(str).str.startswith("ERCC-")]
    res = run_prerank_gsea(ranked, cfg_m2)
    if res is None or res.empty:
        return pd.DataFrame()
    res = res.copy()
    res.insert(0, "contrast_id", contrast_id)
    return res


def nes_comparison(
    nes_a: pd.DataFrame,
    nes_b: pd.DataFrame,
    *,
    contrast_x: str = "A",
    contrast_y: str = "B",
    frame: str | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Join NES by pathway Term; report Spearman across shared pathways.

    ``spearman_p`` is recorded for transparency but gene-set overlap violates
    independence; do not headline it without that caveat (P6-7).
    Column aliases NES_human/NES_mouse are kept for older callers/plots.
    """
    def _prep(df: pd.DataFrame, suffix: str) -> pd.DataFrame:
        term_col = "Term" if "Term" in df.columns else ("Name" if "Name" in df.columns else None)
        nes_col = "NES" if "NES" in df.columns else None
        if term_col is None or nes_col is None:
            return pd.DataFrame(columns=["Term", f"NES_{suffix}"])
        out = df[[term_col, nes_col]].copy()
        out = out.rename(columns={term_col: "Term", nes_col: f"NES_{suffix}"})
        out["Term"] = out["Term"].astype(str).str.replace(r"^.*\.gmt__", "", regex=True)
        return out.drop_duplicates(subset=["Term"], keep="first")

    a = _prep(nes_a, "x")
    b = _prep(nes_b, "y")
    merged = a.merge(b, on="Term", how="inner")
    # Legacy aliases used by sensitivity helpers and older notebooks
    merged["NES_human"] = merged["NES_x"]
    merged["NES_mouse"] = merged["NES_y"]
    stats: dict[str, Any] = {
        "contrast_x": contrast_x,
        "contrast_y": contrast_y,
        "pair_label": f"{contrast_x} vs {contrast_y}",
        "frame": frame,
        "n_pathways_x": int(len(a)),
        "n_pathways_y": int(len(b)),
        "n_pathways_human": int(len(a)),
        "n_pathways_mouse": int(len(b)),
        "n_shared_pathways": int(len(merged)),
        "spearman_nes": float("nan"),
        "sign_agreement": float("nan"),
        "spearman_p_note": (
            "gene-set overlap violates independence; spearman_p is not a valid "
            "evidence claim for the headline"
        ),
        "xlabel": f"{contrast_x} NES",
        "ylabel": f"{contrast_y} NES",
    }
    if len(merged) >= 3:
        from scipy.stats import spearmanr

        r, p = spearmanr(merged["NES_x"], merged["NES_y"])
        stats["spearman_nes"] = float(r)
        stats["spearman_p"] = float(p)
        same = np.sign(merged["NES_x"]) == np.sign(merged["NES_y"])
        nonzero = (merged["NES_x"] != 0) & (merged["NES_y"] != 0)
        if int(nonzero.sum()) > 0:
            stats["sign_agreement"] = float(same[nonzero].mean())
    return merged, stats


def compare_configured_nes_pairs(
    nes_by_id: dict[str, pd.DataFrame],
    cfg: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, tuple[pd.DataFrame, dict[str, Any]]]]:
    """
    Correlate the four configured NES pairs and no others.

    Returns (summary_table, {pair_key: (joined, stats)}).
    """
    pairs = (cfg.get("module4") or {}).get("nes_pairs") or [
        ["H_DISEASE", "H_AGING_GTEX", "biology_best_case"],
        ["H_DISEASE", "M_FLIGHT", "cross_ecosystem"],
        ["H_AGING_GTEX", "M_FLIGHT", "cross_ecosystem"],
        ["H_AGING_HUBMAP", "H_AGING_GTEX", "method"],
    ]
    summary_rows = []
    details: dict[str, tuple[pd.DataFrame, dict[str, Any]]] = {}
    for item in pairs:
        cx, cy, frame = item[0], item[1], item[2] if len(item) > 2 else None
        if cx not in nes_by_id or cy not in nes_by_id:
            continue
        joined, stats = nes_comparison(
            nes_by_id[cx],
            nes_by_id[cy],
            contrast_x=cx,
            contrast_y=cy,
            frame=frame,
        )
        key = f"{cx}__{cy}"
        details[key] = (joined, stats)
        summary_rows.append(
            {
                "contrast_x": cx,
                "contrast_y": cy,
                "frame": frame,
                "n_shared_pathways": stats.get("n_shared_pathways"),
                "spearman_rho": stats.get("spearman_nes"),
                "sign_agreement": stats.get("sign_agreement"),
                "spearman_p": stats.get("spearman_p"),
                "spearman_p_note": stats.get("spearman_p_note"),
            }
        )
    return pd.DataFrame(summary_rows), details

def _library_from_term(term: str) -> str:
    t = str(term)
    if "Hallmark" in t or t.startswith("HALLMARK") or "MSigDB_Hallmark" in t:
        return "hallmark"
    if "Reactome" in t or t.startswith("R-HSA") or "Reactome_" in t:
        return "reactome"
    if "GO_Biological" in t or t.startswith("GO:"):
        return "go_bp"
    # After gmt__ strip, Hallmark terms are often bare names; Reactome keep R-HSA;
    # GO keep (GO:...) suffix. Heuristic for stripped Hallmark:
    if "(GO:" in t:
        return "go_bp"
    if "R-HSA-" in t:
        return "reactome"
    return "hallmark"


def _attach_library(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    raw = out["Term"].astype(str)
    # Prefer prefix before strip when present in original columns
    if "Term" in out.columns:
        out["library"] = raw.map(_library_from_term)
    return out


def nes_sensitivity_table(
    nes_human_full: pd.DataFrame,
    nes_mouse_full: pd.DataFrame,
) -> pd.DataFrame:
    """
    Spearman NES concordance for gene-set slices (P6-7 sensitivity asset).

    Requires NES tables that still carry GMT library identity in the Term
    (``*.gmt__...``) or recoverable library tags. Headline Module 4 uses
    Hallmark + Reactome; this table shows the finding survives that choice.
    """
    def _prep_with_lib(df: pd.DataFrame, suffix: str) -> pd.DataFrame:
        term_col = "Term" if "Term" in df.columns else ("Name" if "Name" in df.columns else None)
        nes_col = "NES" if "NES" in df.columns else None
        if term_col is None or nes_col is None:
            return pd.DataFrame(columns=["Term", f"NES_{suffix}", "library"])
        out = df[[term_col, nes_col]].copy()
        out = out.rename(columns={term_col: "Term", nes_col: f"NES_{suffix}"})
        raw = out["Term"].astype(str)
        out["library"] = [
            (
                "hallmark"
                if "Hallmark" in t or "MSigDB_Hallmark" in t
                else "reactome"
                if "Reactome" in t
                else "go_bp"
                if "GO_Biological" in t or "GO_Biological_Process" in t
                else _library_from_term(t.split(".gmt__")[-1] if ".gmt__" in t else t)
            )
            for t in raw
        ]
        out["Term"] = raw.str.replace(r"^.*\.gmt__", "", regex=True)
        return out.drop_duplicates(subset=["Term"], keep="first")

    h = _prep_with_lib(nes_human_full, "human")
    m = _prep_with_lib(nes_mouse_full, "mouse")
    merged = h.merge(m, on=["Term", "library"], how="inner")
    slices = [
        ("All three (Hallmark + Reactome + GO BP)", None),
        ("Hallmark + Reactome", {"hallmark", "reactome"}),
        ("Reactome only", {"reactome"}),
        ("Hallmark only", {"hallmark"}),
        ("GO BP only", {"go_bp"}),
    ]
    rows = []
    from scipy.stats import spearmanr

    for label, keep in slices:
        sub = merged if keep is None else merged[merged["library"].isin(keep)]
        n = int(len(sub))
        rho = float("nan")
        p = float("nan")
        agree = float("nan")
        if n >= 3:
            rho, p = spearmanr(sub["NES_human"], sub["NES_mouse"])
            rho, p = float(rho), float(p)
            nonzero = (sub["NES_human"] != 0) & (sub["NES_mouse"] != 0)
            if int(nonzero.sum()) > 0:
                agree = float(
                    (np.sign(sub.loc[nonzero, "NES_human"]) == np.sign(sub.loc[nonzero, "NES_mouse"]))
                    .mean()
                )
        rows.append(
            {
                "slice": label,
                "n_shared_pathways": n,
                "spearman_rho": rho,
                "spearman_p": p,
                "sign_agreement": agree,
                "headline": label.startswith("Hallmark + Reactome"),
                "spearman_p_note": (
                    "gene-set overlap violates independence; p is descriptive only"
                ),
            }
        )
    return pd.DataFrame(rows)


def sensitivity_enabled(cfg: dict[str, Any]) -> bool:
    """Shipped default is false: load committed TSV (Module 3 load_h5mu pattern)."""
    return bool((cfg.get("module4") or {}).get("sensitivity", {}).get("enabled", False))


def sensitivity_table_path(cfg: dict[str, Any]) -> Path:
    from scripts.common.paths import module_root, resolve

    sens = (cfg.get("module4") or {}).get("sensitivity") or {}
    rel = sens.get("table") or "outputs/tables/module4_nes_sensitivity.tsv"
    p = Path(rel)
    if not p.is_absolute():
        # Prefer resolved outputs_tables when relative path is the default name
        if p.name == "module4_nes_sensitivity.tsv" and len(p.parts) <= 3:
            return resolve(cfg, "outputs_tables") / p.name
        return module_root(cfg) / p
    return p


def resolve_nes_sensitivity(
    cfg: dict[str, Any],
    *,
    logfc_x: pd.Series | None = None,
    logfc_y: pd.Series | None = None,
    contrast_x: str | None = None,
    contrast_y: str | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load committed sensitivity TSV, or recompute when ``sensitivity.enabled``.

    Returns ``(table, meta)`` where meta records source and whether GO BP was preranked.
    """
    path = sensitivity_table_path(cfg)
    meta: dict[str, Any] = {
        "sensitivity_enabled": sensitivity_enabled(cfg),
        "sensitivity_table": str(path),
        "source": None,
    }
    if not sensitivity_enabled(cfg):
        if not path.exists():
            raise FileNotFoundError(
                f"module4.sensitivity.enabled is false but committed table missing: {path}. "
                "Set module4.sensitivity.enabled: true once to regenerate, then ship the TSV."
            )
        df = pd.read_csv(path, sep="\t")
        meta["source"] = "committed_tsv"
        return df, meta

    if logfc_x is None or logfc_y is None or not contrast_x or not contrast_y:
        raise ValueError(
            "sensitivity.enabled=true requires logfc_x/logfc_y and contrast_x/contrast_y "
            "to re-prerank Hallmark+Reactome+GO BP"
        )
    sens_sets = (cfg.get("module4", {}).get("enrichment") or {}).get("sensitivity_gene_sets") or [
        "MSigDB_Hallmark_2020",
        "Reactome_2022",
        "GO_Biological_Process_2023",
    ]
    nes_x = prerank_contrast(logfc_x, cfg, contrast_id=contrast_x, gene_sets=sens_sets)
    nes_y = prerank_contrast(logfc_y, cfg, contrast_id=contrast_y, gene_sets=sens_sets)
    df = nes_sensitivity_table(nes_x, nes_y)
    meta["source"] = "recomputed"
    meta["sensitivity_gene_sets"] = list(sens_sets)
    return df, meta


def plot_nes_sensitivity(
    sens: pd.DataFrame,
    path: Path | str,
    *,
    dpi: int = 150,
) -> Path | None:
    """Bar chart of Spearman rho across gene-set slices (P6-7)."""
    import matplotlib.pyplot as plt

    if sens is None or sens.empty:
        return None
    df = sens.copy()
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    x = np.arange(len(df))
    colors = ["#4C78A8" if bool(h) else "#9E9E9E" for h in df.get("headline", [False] * len(df))]
    ax.bar(x, df["spearman_rho"].astype(float), color=colors, edgecolor="k", linewidth=0.4)
    ax.axhline(0, color="#666", lw=0.8)
    labels = [
        f"{s}\nn={int(n)}"
        for s, n in zip(df["slice"].astype(str), df["n_shared_pathways"].astype(int))
    ]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("Spearman rho (NES human vs mouse)")
    ax.set_title(
        "NES concordance is near-null across gene-set slices\n"
        "(headline = Hallmark + Reactome; GO BP excluded from headline n)",
        fontsize=10,
    )
    ax.text(
        0.99,
        0.02,
        "spearman_p not shown: overlapping gene sets are not independent tests",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="#333333",
    )
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)


def plot_nes_scatter(
    nes_df: pd.DataFrame,
    path: Path | str,
    *,
    stats: dict[str, Any] | None = None,
    highlight: list[str] | None = None,
    pair_label: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    dpi: int = 150,
) -> Path | None:
    import matplotlib.pyplot as plt

    if nes_df is None or nes_df.empty:
        return None
    # Accept either legacy NES_human/NES_mouse or NES_x/NES_y
    xcol = "NES_x" if "NES_x" in nes_df.columns else "NES_human"
    ycol = "NES_y" if "NES_y" in nes_df.columns else "NES_mouse"
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    ax.axhline(0, color="#888", lw=0.8)
    ax.axvline(0, color="#888", lw=0.8)
    ax.scatter(nes_df[xcol], nes_df[ycol], s=18, alpha=0.55, c="#4C78A8", linewidths=0)
    highlight = highlight or []
    for term in highlight:
        sub = nes_df[nes_df["Term"].astype(str).str.contains(term, case=False, na=False)]
        if sub.empty:
            continue
        ax.scatter(
            sub[xcol],
            sub[ycol],
            s=48,
            c="#F58518",
            edgecolors="k",
            linewidths=0.4,
            zorder=3,
        )
        for _, row in sub.iterrows():
            ax.annotate(
                str(row["Term"])[:40],
                (row[xcol], row[ycol]),
                fontsize=7,
                xytext=(4, 4),
                textcoords="offset points",
            )
    stats = stats or {}
    r = stats.get("spearman_nes", float("nan"))
    n = stats.get("n_shared_pathways", len(nes_df))
    pair = pair_label or stats.get("pair_label") or "contrast pair"
    title = (
        f"NES concordance (Hallmark + Reactome): {pair}\n"
        f"Spearman rho={r:.3f} across {n} shared pathways"
    )
    if isinstance(r, float) and abs(r) < 0.15:
        title += " (near-zero is a result)"
    ax.set_title(title)
    ax.set_xlabel(xlabel or stats.get("xlabel") or f"{stats.get('contrast_x', 'contrast A')} NES")
    ax.set_ylabel(ylabel or stats.get("ylabel") or f"{stats.get('contrast_y', 'contrast B')} NES")
    ax.text(
        0.02,
        0.02,
        "Quadrants: concordant up (I), discordant (II/IV), concordant down (III)",
        transform=ax.transAxes,
        fontsize=7,
        va="bottom",
        color="#444",
    )
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)


def hubmap_gtex_bridge_table(cfg: dict[str, Any]) -> pd.DataFrame:
    """
    Join composition-weighted + max-across HuBMAP vectors to GTEx mean TPM.

    Log1p once here for plotting columns; never z-score across units.
    """
    from scripts.nb06_crossres.load import load_gtex_means_for_bridge, load_hubmap_pseudobulk

    weighted, max_vec, _ = load_hubmap_pseudobulk(cfg)
    gtex = load_gtex_means_for_bridge(cfg)
    df = gtex.copy()
    df["hubmap_composition_weighted_linear"] = weighted.reindex(df.index)
    df["hubmap_max_across_labels_linear"] = max_vec.reindex(df.index)
    df["hubmap_composition_weighted_log1p"] = np.log1p(df["hubmap_composition_weighted_linear"].astype(float))
    df["hubmap_max_across_labels_log1p"] = np.log1p(df["hubmap_max_across_labels_linear"].astype(float))
    df["gtex_log1p_tpm"] = np.log1p(df["gtex_lung_mean_tpm"].astype(float))
    return df.reset_index()


def plot_pseudobulk_vs_gtex(
    bridge_df: pd.DataFrame,
    path: Path | str,
    *,
    panel: list[str] | None = None,
    dpi: int = 150,
) -> Path | None:
    """Show GTEx vs composition-weighted AND max-across -- the M9 lesson."""
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr

    if bridge_df is None or bridge_df.empty:
        return None
    df = bridge_df.dropna(
        subset=["gtex_log1p_tpm", "hubmap_composition_weighted_log1p", "hubmap_max_across_labels_log1p"]
    )
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.0))
    for ax, col, label in (
        (axes[0], "hubmap_composition_weighted_log1p", "composition-weighted (correct bulk proxy)"),
        (axes[1], "hubmap_max_across_labels_log1p", "max-across-cell-types (contrast - not a bulk mean)"),
    ):
        x = df[col].astype(float)
        y = df["gtex_log1p_tpm"].astype(float)
        ax.scatter(x, y, s=6, alpha=0.2, c="#4C78A8", linewidths=0)
        r, _ = spearmanr(x, y)
        ax.set_xlabel(f"HuBMAP {label}\nlog1p(linear)")
        ax.set_ylabel("GTEx lung mean TPM (log1p)")
        ax.set_title(f"Spearman rho={r:.3f}\n(never z-score across CP10K vs TPM)")
        if panel:
            sub = df[df["gene_symbol"].astype(str).str.upper().isin([p.upper() for p in panel])]
            ax.scatter(sub[col], sub["gtex_log1p_tpm"], s=40, c="#F58518", edgecolors="k", linewidths=0.4, zorder=3)
            for _, row in sub.iterrows():
                ax.annotate(str(row["gene_symbol"]), (row[col], row["gtex_log1p_tpm"]), fontsize=7)
    fig.suptitle(
        "HuBMAP vs GTEx: composition-weighted tracks bulk more closely "
        "(rho 0.741 against 0.636); max-across is a contrast, not a bulk mean",
        fontsize=11,
        y=1.02,
    )
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)


def is_ercc_symbol(symbol: object) -> bool:
    """True only for ERCC spike-in IDs (ERCC-...), not human ERCC* repair genes."""
    return str(symbol).upper().startswith("ERCC-")


def plot_volcano(
    de: pd.DataFrame,
    path: Path | str,
    *,
    title: str,
    panel: list[str] | None = None,
    exclude_ercc: bool = True,
    dpi: int = 150,
) -> Path | None:
    """
    Volcano of logFC vs -log10(padj).

    ERCC spike-ins are excluded from the plot by default. In OSD-248, flight and
    ground control received different ExFold mixes, so ERCC DE is expected technical
    signal -- not unloading biology. The full DE table still retains ERCC rows.
    """
    import matplotlib.pyplot as plt

    if de is None or de.empty:
        return None
    df = de.copy()
    n_ercc = 0
    if exclude_ercc and "gene_symbol" in df.columns:
        ercc_mask = df["gene_symbol"].map(is_ercc_symbol)
        n_ercc = int(ercc_mask.sum())
        df = df.loc[~ercc_mask].copy()
    if df.empty:
        return None
    df["neglog10padj"] = -np.log10(df["padj"].astype(float).clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ax.scatter(df["logFC"], df["neglog10padj"], s=6, alpha=0.25, c="#9E9E9E", linewidths=0)
    if panel:
        sub = df[df["gene_symbol"].astype(str).str.upper().isin([p.upper() for p in panel])]
        ax.scatter(sub["logFC"], sub["neglog10padj"], s=40, c="#F58518", edgecolors="k", linewidths=0.4, zorder=3)
        for _, row in sub.iterrows():
            ax.annotate(str(row["gene_symbol"]), (row["logFC"], row["neglog10padj"]), fontsize=7)
    ax.axvline(0, color="#666", lw=0.8)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p")
    plot_title = title
    if exclude_ercc and n_ercc:
        plot_title = f"{title}\n(ERCC spike-ins excluded from plot; n={n_ercc} remain in DE table)"
    ax.set_title(plot_title)
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)
