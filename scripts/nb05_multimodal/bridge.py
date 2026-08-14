"""True-multiome summaries: barcode pairing, modality contrast, RNA<->ATAC gene bridge."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from scripts.nb03_integration.gtex import strip_ensembl_version
from scripts.nb05_multimodal.load import barcode_sets, modality_shapes


def barcode_overlap_table(cfg: dict[str, Any]) -> pd.DataFrame:
    """Confirm shared-nucleus pairing across MuData modalities."""
    sets = barcode_sets(cfg)
    rna = sets.get("rna", set())
    atac = sets.get("atac_cbg", sets.get("atac_cbb", set()))
    inter = rna & atac
    n_inter = len(inter)
    paired = n_inter == len(rna) == len(atac) and n_inter > 0
    return pd.DataFrame(
        [
            {
                "rna_n_barcodes": len(rna),
                "atac_n_barcodes": len(atac),
                "n_exact_overlap": n_inter,
                "overlap_fraction_of_rna": n_inter / max(len(rna), 1),
                "overlap_fraction_of_atac": n_inter / max(len(atac), 1),
                "paired_multiome": paired,
                "interpretation": (
                    "True multiome: RNA and ATAC barcode sets are identical (shared nuclei)."
                    if paired
                    else "Barcode spaces differ -- inspect before claiming shared-nucleus multiome."
                ),
            }
        ]
    )


def modality_contrast_table(cfg: dict[str, Any]) -> pd.DataFrame:
    m5 = cfg["module3"]
    load_h5mu = bool((m5.get("params") or {}).get("load_h5mu", False))
    rows = [
        {
            "aspect": "Donor",
            "value": f"{m5['donor_id']} ({m5['donor_label']})",
        },
        {
            "aspect": "Dataset",
            "value": m5["dataset_id"],
        },
        {
            "aspect": "Portal URL",
            "value": m5.get("portal_url", "https://portal.hubmapconsortium.org"),
        },
        {
            "aspect": "Assay / pipeline",
            "value": m5.get("assay", "SNARE-seq2 [Salmon + ArchR + Muon]"),
        },
        {
            "aspect": "Design",
            "value": "True multiome (shared-nucleus RNA + ATAC)",
        },
        {
            "aspect": "Relation to Modules 1-2 / 4",
            "value": (
                "Those modules use 10x snRNA-seq; Module 3 is the configured "
                f"donor SNARE true multiome ({m5['donor_label']})"
            ),
        },
    ]
    if load_h5mu:
        try:
            shapes = modality_shapes(cfg).set_index("modality")
            for mod in ("rna", "atac_cbg", "atac_cbb"):
                if mod not in shapes.index:
                    continue
                r = shapes.loc[mod]
                rows.append(
                    {
                        "aspect": f"{mod} shape",
                        "value": (
                            f"{int(r['n_nuclei'])} nuclei x {int(r['n_features'])} features "
                            f"({r['matrix_type']})"
                        ),
                    }
                )
        except FileNotFoundError as exc:
            rows.append({"aspect": "MuData shapes", "value": f"unavailable: {exc}"})
    else:
        rows.append(
            {
                "aspect": "MuData shapes",
                "value": "skipped (load_h5mu=false); see MOFA n_nuclei in run_params",
            }
        )
    rows.append(
        {
            "aspect": "MOFA views",
            "value": "rna + atac_cbg (30 factors in multiome_mofa.hdf5; trailing factors may be dead)",
        }
    )
    rows.append(
        {
            "aspect": "Gene bridge metric",
            "value": m5.get("analysis", {}).get(
                "bridge_metric", "feature_std_on_mofa_inputs"
            ),
        }
    )
    rows.append(
        {
            "aspect": "Teaching note",
            "value": (
                "Default path interprets multiome_mofa.hdf5 plus a committed "
                "Azimuth/WNN label cache; secondary_analysis.h5mu is optional "
                "(load_h5mu) and is not a learner download"
            ),
        }
    )
    return pd.DataFrame(rows)

def rna_atac_gene_bridge(mofa: dict[str, Any]) -> pd.DataFrame:
    """
    Gene bridge on MOFA feature space.

    MOFA RNA features are Ensembl IDs (+ hugo metadata); ATAC features are HUGO
    symbols. Inputs are centered, so we compare **feature standard deviations**
    (variability present in both views) rather than means.
    """
    rna_ids = strip_ensembl_version(pd.Index(mofa["features"]["rna"].astype(str))).astype(str)
    rna_df = pd.DataFrame(
        {
            "gene_id_no_version": rna_ids.to_numpy(),
            "gene_symbol": np.asarray(mofa["rna_hugo"].astype(str)),
            "rna_std": np.asarray(mofa["rna_std"], dtype=float),
            "rna_mean_expression": np.asarray(mofa["rna_means"], dtype=float),
        }
    )
    rna_df["_sym"] = rna_df["gene_symbol"].str.upper()
    rna_df = rna_df[rna_df["_sym"].ne("NAN") & rna_df["_sym"].ne("")]
    rna_df = (
        rna_df.sort_values("rna_std", ascending=False)
        .groupby("_sym", as_index=False, sort=False)
        .first()
    )

    atac_df = pd.DataFrame(
        {
            "gene_symbol_atac": np.asarray(mofa["features"]["atac_cbg"].astype(str)),
            "atac_std": np.asarray(mofa["atac_std"], dtype=float),
            "atac_mean_gene_activity": np.asarray(mofa["atac_means"], dtype=float),
        }
    )
    atac_df["_sym"] = atac_df["gene_symbol_atac"].str.upper()
    atac_df = atac_df[atac_df["_sym"].ne("NAN") & atac_df["_sym"].ne("")]
    atac_df = (
        atac_df.sort_values("atac_std", ascending=False)
        .groupby("_sym", as_index=False, sort=False)
        .first()
    )

    merged = rna_df.merge(atac_df[["_sym", "atac_std", "atac_mean_gene_activity"]], on="_sym", how="inner")
    df = pd.DataFrame(
        {
            "gene_symbol": merged["gene_symbol"].astype(str).values,
            "gene_id_no_version": merged["gene_id_no_version"].astype(str).values,
            "rna_std": merged["rna_std"].astype(float).values,
            "atac_std": merged["atac_std"].astype(float).values,
            # keep mean columns for focus-bar compatibility / inspection
            "rna_mean_expression": merged["rna_mean_expression"].astype(float).values,
            "atac_mean_gene_activity": merged["atac_mean_gene_activity"].astype(float).values,
        }
    )

    x = np.log1p(df["rna_std"].to_numpy())
    y = np.log1p(df["atac_std"].to_numpy())
    if len(df) >= 3 and np.std(x) > 0 and np.std(y) > 0:
        df.attrs["pearson_log1p"] = float(np.corrcoef(x, y)[0, 1])
        from scipy.stats import spearmanr

        df.attrs["spearman_log1p"] = float(spearmanr(x, y).correlation)
    else:
        df.attrs["pearson_log1p"] = float("nan")
        df.attrs["spearman_log1p"] = float("nan")
    df.attrs["bridge_metric"] = "feature_std_on_mofa_inputs"
    return df


def focus_gene_bridge(bridge_df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    focus = [g.upper() for g in (cfg["module3"]["analysis"].get("focus_genes") or [])]
    if not focus or bridge_df.empty:
        return pd.DataFrame()
    df = bridge_df.copy()
    df["_sym"] = df["gene_symbol"].astype(str).str.upper()
    out = df[df["_sym"].isin(focus)].drop(columns=["_sym"]).copy()
    order = {g: i for i, g in enumerate(focus)}
    out["_ord"] = out["gene_symbol"].astype(str).str.upper().map(order)
    return out.sort_values("_ord").drop(columns=["_ord"]).reset_index(drop=True)


def focus_gene_availability(
    cfg: dict[str, Any],
    focus_long: pd.DataFrame,
    bridge_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Document which configured focus genes were recovered in both modalities."""
    analysis = cfg["module3"]["analysis"]
    focus = [g.upper() for g in (analysis.get("focus_genes") or [])]
    documented_missing = {
        g.upper() for g in (analysis.get("focus_genes_not_in_both_modalities") or [])
    }
    recovered = set()
    if focus_long is not None and not focus_long.empty:
        recovered |= set(focus_long["gene_symbol"].astype(str).str.upper())
    if bridge_df is not None and not bridge_df.empty:
        recovered |= set(bridge_df["gene_symbol"].astype(str).str.upper())

    rows = []
    for g in focus:
        in_both = g in recovered
        rows.append(
            {
                "gene_symbol": g,
                "in_focus_panel": True,
                "recovered_in_both_modalities": in_both,
                "notes": (
                    "Recovered in paired MuData RNA + ATAC gene activity"
                    if in_both
                    else "Not recovered in both modalities used by Module 3"
                ),
            }
        )
    for g in sorted(documented_missing - set(focus)):
        rows.append(
            {
                "gene_symbol": g,
                "in_focus_panel": False,
                "recovered_in_both_modalities": False,
                "notes": (
                    "Documented omission: not present in both MOFA/MuData ATAC gene spaces "
                    "(e.g. RNA-only epithelial marker)"
                ),
            }
        )
    return pd.DataFrame(rows)


def focus_from_paired_matrix(focus_long: pd.DataFrame) -> pd.DataFrame:
    """Pseudobulk means from paired per-nucleus matrix.

    On the MOFA path, values are centered inputs: keep them for SD / correlation
    bridging only (column names reflect that) and flag ``means_are_centered``.
    """
    if focus_long is None or focus_long.empty:
        return pd.DataFrame()
    centered = bool(getattr(focus_long, "attrs", {}).get("means_are_centered", False))
    source = str(getattr(focus_long, "attrs", {}).get("source", ""))
    if "mofa" in source.lower():
        centered = True
    rows = []
    for sym, sub in focus_long.groupby("gene_symbol", sort=False):
        rna_m = float(sub["rna_value"].mean())
        atac_m = float(sub["atac_value"].mean())
        if centered:
            rows.append(
                {
                    "gene_symbol": sym,
                    "gene_id_no_version": sub["gene_id_no_version"].iloc[0],
                    # Kept for SD bridge / correlation only - not expression means
                    "rna_centered_mean": rna_m,
                    "atac_centered_mean": atac_m,
                    "rna_mean_expression": rna_m,
                    "atac_mean_gene_activity": atac_m,
                }
            )
        else:
            rows.append(
                {
                    "gene_symbol": sym,
                    "gene_id_no_version": sub["gene_id_no_version"].iloc[0],
                    "rna_mean_expression": rna_m,
                    "atac_mean_gene_activity": atac_m,
                }
            )
    out = pd.DataFrame(rows)
    out.attrs["means_are_centered"] = centered
    out.attrs["source"] = getattr(focus_long, "attrs", {}).get("source", "")
    return out


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _correlation_note(
    *,
    n_nuclei: int,
    rna_mean: float,
    pearson_r: float,
    min_n: int,
) -> str:
    notes = []
    if n_nuclei < int(min_n):
        notes.append(f"n_nuclei={n_nuclei} below threshold {min_n}; correlation not interpretable")
    if rna_mean == 0.0 or (isinstance(rna_mean, float) and rna_mean == 0):
        notes.append("rna_mean==0 (zero is a finding -- marker absent in this RNA layer/subset)")
    if pearson_r != pearson_r:  # nan
        if np.std is not None:
            notes.append("pearson_r undefined (zero variance in RNA and/or ATAC)")
    return "; ".join(notes) if notes else ""


def focus_gene_correlations(
    focus_long: pd.DataFrame,
    *,
    min_nuclei: int = 50,
) -> pd.DataFrame:
    """Per-gene Pearson correlation across all paired nuclei (expected to be weak)."""
    if focus_long is None or focus_long.empty:
        return pd.DataFrame()
    rows = []
    for sym, sub in focus_long.groupby("gene_symbol", sort=False):
        x = sub["rna_value"].to_numpy(dtype=float)
        y = sub["atac_value"].to_numpy(dtype=float)
        r = _pearson(x, y)
        rna_mean = float(np.mean(x))
        n = int(len(sub))
        rows.append(
            {
                "gene_symbol": sym,
                "gene_id_no_version": sub["gene_id_no_version"].iloc[0],
                "scope": "global",
                "cell_type_filter": "all_nuclei",
                "n_nuclei": n,
                "pearson_r": r,
                "rna_mean": rna_mean,
                "atac_mean": float(np.mean(y)),
                "note": _correlation_note(
                    n_nuclei=n, rna_mean=rna_mean, pearson_r=r, min_n=min_nuclei
                ),
                "interpretation": (
                    "Weak global r is expected: gene activity != expression; "
                    "marker signal diluted across cell types"
                ),
            }
        )
    return pd.DataFrame(rows)


def focus_gene_correlations_by_cell_type(
    focus_long: pd.DataFrame,
    cfg: dict[str, Any],
    label_col: str = "azimuth_label",
) -> pd.DataFrame:
    """
    Within-cell-type RNA<->ATAC Pearson for focus genes.

    Emits explicit ``note`` when n is low or rna_mean==0 (a zero is a finding).
    """
    if focus_long is None or focus_long.empty:
        return pd.DataFrame()
    if label_col not in focus_long.columns:
        return pd.DataFrame()

    mapping = cfg["module3"]["analysis"].get("focus_gene_cell_types") or {}
    min_n = int((cfg["module3"].get("thresholds") or {}).get("min_nuclei_correlation", 50))
    rows = []
    for sym, sub in focus_long.groupby("gene_symbol", sort=False):
        key = str(sym).upper()
        types = mapping.get(key) or mapping.get(str(sym)) or []
        if not types:
            continue
        mask = sub[label_col].astype(str).isin([str(t) for t in types])
        restricted = sub.loc[mask]
        x = restricted["rna_value"].to_numpy(dtype=float)
        y = restricted["atac_value"].to_numpy(dtype=float)
        rna_mean = float(np.mean(x)) if len(x) else float("nan")
        atac_mean = float(np.mean(y)) if len(y) else float("nan")
        r = _pearson(x, y)
        n = int(len(restricted))
        rows.append(
            {
                "gene_symbol": sym,
                "gene_id_no_version": sub["gene_id_no_version"].iloc[0],
                "scope": "within_cell_type",
                "cell_type_filter": ";".join(map(str, types)),
                "n_nuclei": n,
                "pearson_r": r,
                "rna_mean": rna_mean,
                "atac_mean": atac_mean,
                "note": _correlation_note(
                    n_nuclei=n,
                    rna_mean=0.0 if rna_mean != rna_mean else rna_mean,
                    pearson_r=r,
                    min_n=min_n,
                ),
                "interpretation": (
                    "Restricting to one cell type removes between-type variance; "
                    "within-type r is usually smaller and is a stricter test"
                ),
            }
        )
    return pd.DataFrame(rows)


def combine_focus_correlations(
    global_df: pd.DataFrame,
    celltype_df: pd.DataFrame,
) -> pd.DataFrame:
    frames = [f for f in (global_df, celltype_df) if f is not None and not f.empty]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def mofa_factor_dominance(
    variance_per_factor: pd.DataFrame,
    n_factors: int = 14,
) -> pd.DataFrame:
    """Summarize modality dominance; default truncate at Factor 14 (dead beyond)."""
    if variance_per_factor is None or variance_per_factor.empty:
        return pd.DataFrame()
    sub = variance_per_factor.iloc[:n_factors].copy()
    views = list(sub.columns)
    rows = []
    for factor, row in sub.iterrows():
        vals = {v: float(row[v]) for v in views}
        dominant = max(vals, key=vals.get)
        nondom = sorted(vals.values(), reverse=True)
        nondom_r2 = float(nondom[1]) if len(nondom) > 1 else 0.0
        rows.append(
            {
                "factor": factor,
                **vals,
                "dominant_view": dominant,
                "nondominant_r2": nondom_r2,
                "rna_minus_atac": float(vals.get("rna", 0.0) - vals.get("atac_cbg", 0.0)),
            }
        )
    return pd.DataFrame(rows)


DEFAULT_CONCORDANCE_PAIRS: list[tuple[str, str]] = [
    ("atac_leiden", "azimuth_label"),
    ("atac_clusters", "azimuth_label"),
    ("leiden_wnn", "rna_leiden"),
    ("leiden_wnn", "atac_leiden"),
    ("rna_leiden", "azimuth_label"),
    ("leiden_wnn", "azimuth_label"),
    ("atac_leiden", "rna_leiden"),
]


def cluster_concordance(
    frame: pd.DataFrame,
    pairs: list[tuple[str, str]] | None = None,
) -> pd.DataFrame:
    """Pairwise ARI / NMI for named label columns."""
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    if frame is None or frame.empty:
        return pd.DataFrame()
    pairs = pairs or DEFAULT_CONCORDANCE_PAIRS
    rows = []
    for col_a, col_b in pairs:
        if col_a not in frame.columns or col_b not in frame.columns:
            rows.append(
                {
                    "pair": f"{col_a} vs {col_b}",
                    "col_a": col_a,
                    "col_b": col_b,
                    "n_nuclei": 0,
                    "k_a": 0,
                    "k_b": 0,
                    "ari": float("nan"),
                    "nmi": float("nan"),
                    "note": "missing column",
                }
            )
            continue
        sub = frame[[col_a, col_b]].dropna()
        a = sub[col_a].astype(str).to_numpy()
        b = sub[col_b].astype(str).to_numpy()
        rows.append(
            {
                "pair": f"{col_a} vs {col_b}",
                "col_a": col_a,
                "col_b": col_b,
                "n_nuclei": int(len(sub)),
                "k_a": int(pd.Series(a).nunique()),
                "k_b": int(pd.Series(b).nunique()),
                "ari": float(adjusted_rand_score(a, b)) if len(sub) else float("nan"),
                "nmi": float(normalized_mutual_info_score(a, b)) if len(sub) else float("nan"),
                "note": "",
            }
        )
    return pd.DataFrame(rows)


def _eta_squared(values: np.ndarray, groups: np.ndarray) -> float:
    mask = np.isfinite(values) & pd.notna(groups)
    y = np.asarray(values[mask], dtype=float)
    g = np.asarray(groups[mask], dtype=object)
    if len(y) < 3:
        return float("nan")
    grand = float(np.mean(y))
    ss_tot = float(np.sum((y - grand) ** 2))
    if ss_tot <= 0:
        return 0.0
    ss_between = 0.0
    for label in pd.unique(g):
        yg = y[g == label]
        ss_between += float(len(yg)) * (float(np.mean(yg)) - grand) ** 2
    return ss_between / ss_tot


_LOG10_QC = frozenset({"nFrags", "ReadsInTSS", "ReadsInPromoter"})


def qc_association(
    frame: pd.DataFrame,
    label_cols: list[str],
    qc_cols: list[str],
) -> pd.DataFrame:
    """Eta squared of each QC metric under each clustering."""
    if frame is None or frame.empty:
        return pd.DataFrame()
    rows = []
    for qc in qc_cols:
        if qc not in frame.columns:
            continue
        raw = pd.to_numeric(frame[qc], errors="coerce").to_numpy(dtype=float)
        if qc in _LOG10_QC:
            metric = np.log10(np.clip(raw, a_min=1.0, a_max=None))
            transform = "log10"
            metric_name = f"log10 {qc}"
        else:
            metric = raw
            transform = "none"
            metric_name = qc
        for lab in label_cols:
            if lab not in frame.columns:
                continue
            groups = frame[lab].to_numpy()
            rows.append(
                {
                    "qc_metric": metric_name,
                    "qc_column": qc,
                    "labelling": lab,
                    "eta_squared": _eta_squared(metric, groups),
                    "transform": transform,
                    "n_nuclei": int(np.isfinite(metric).sum()),
                    "n_groups": int(pd.Series(groups).nunique(dropna=True)),
                }
            )
    return pd.DataFrame(rows)


def factor_label_association(
    Z: pd.DataFrame,
    frame: pd.DataFrame,
    label_cols: list[str],
) -> pd.DataFrame:
    """Eta squared of every MOFA factor score under every labelling."""
    if Z is None or Z.empty or frame is None or frame.empty:
        return pd.DataFrame()
    aligned = Z.join(frame[label_cols], how="inner")
    rows = []
    for factor in Z.columns:
        y = aligned[factor].to_numpy(dtype=float)
        for lab in label_cols:
            groups = aligned[lab].to_numpy()
            rows.append(
                {
                    "factor": factor,
                    "labelling": lab,
                    "eta_squared": _eta_squared(y, groups),
                    "n_groups": int(pd.Series(groups).nunique(dropna=True)),
                    "n_nuclei": int(len(aligned)),
                }
            )
    return pd.DataFrame(rows)


def azimuth_composition(frame: pd.DataFrame, label_col: str = "azimuth_label") -> pd.DataFrame:
    """Azimuth label counts and fractions for the cross-assay comparison."""
    if frame is None or frame.empty or label_col not in frame.columns:
        return pd.DataFrame()
    counts = frame[label_col].astype(str).value_counts(dropna=False)
    out = counts.rename("n_nuclei").rename_axis(label_col).reset_index()
    out["fraction"] = out["n_nuclei"] / out["n_nuclei"].sum()
    return out


def barcode_overlap_from_mofa(mofa: dict[str, Any]) -> pd.DataFrame:
    """When MuData is absent: MOFA samples are the paired nucleus set by construction."""
    n = int(mofa["factors"].shape[0])
    return pd.DataFrame(
        [
            {
                "rna_n_barcodes": n,
                "atac_n_barcodes": n,
                "n_exact_overlap": n,
                "overlap_fraction_of_rna": 1.0,
                "overlap_fraction_of_atac": 1.0,
                "paired_multiome": True,
                "interpretation": (
                    "NOT independent pairing evidence from a barcode crosstab. "
                    "Pairing is asserted by MOFA file construction: one shared sample "
                    "barcode list is written once for rna + atac_cbg views. "
                    "n_exact_overlap equals that list length by definition "
                    "(MuData barcode crosstab skipped because load_h5mu=false)."
                ),
            }
        ]
    )
