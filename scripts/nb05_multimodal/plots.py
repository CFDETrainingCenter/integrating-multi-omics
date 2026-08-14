"""Module 3 plotting helpers for paired SNARE multiome (MOFA interpretation)."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.common.plotting import save_figure


def _scatter_embedding(ax, xy, color_vals, title: str, cmap: str = "tab20"):
    if color_vals is None:
        ax.scatter(xy[:, 0], xy[:, 1], s=3, c="#4C78A8", linewidths=0)
    else:
        cats = pd.Series(color_vals).astype("category")
        codes = cats.cat.codes.to_numpy()
        ax.scatter(xy[:, 0], xy[:, 1], c=codes, s=3, cmap=cmap, linewidths=0)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_xticks([])
    ax.set_yticks([])


def plot_umap_panel(
    adata,
    path: Path | str,
    dpi: int = 150,
    *,
    donor_label: str | None = None,
) -> Path | None:
    """Joint / RNA / ATAC UMAPs with consistent coloring (requires MuData embeddings)."""
    import matplotlib.pyplot as plt

    emb_keys = []
    if "X_umap" in adata.obsm:
        emb_keys.append(("Joint MuData", "X_umap"))
    if "X_umap_rna" in adata.obsm:
        emb_keys.append(("RNA", "X_umap_rna"))
    if "X_umap_atac" in adata.obsm:
        emb_keys.append(("ATAC gene-activity", "X_umap_atac"))
    if not emb_keys:
        return None

    label_rows = []
    if "azimuth_label" in adata.obs.columns:
        label_rows.append(("azimuth_label", "Azimuth label"))
    if "leiden" in adata.obs.columns:
        label_rows.append(("leiden", "RNA Leiden"))
    if not label_rows:
        label_rows.append((None, "no label"))

    n_rows = len(label_rows)
    n_cols = len(emb_keys)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.0 * n_cols, 3.6 * n_rows), squeeze=False)
    for r, (color_key, color_name) in enumerate(label_rows):
        vals = adata.obs[color_key] if color_key and color_key in adata.obs.columns else None
        for c, (title, key) in enumerate(emb_keys):
            xy = np.asarray(adata.obsm[key])
            _scatter_embedding(
                axes[r, c],
                xy,
                vals,
                f"{title}\n({color_name})",
            )
    donor = donor_label or "configured donor"
    fig.suptitle(
        f"{donor} SNARE multiome embeddings (same label across modalities)",
        y=1.01,
        fontsize=11,
    )
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)


def plot_mofa_umap(
    adata,
    path: Path | str,
    color: str = "azimuth_label",
    dpi: int = 150,
    *,
    donor_label: str | None = None,
) -> Path | None:
    """UMAP/PCA view of MOFA factors -- do not claim a joint shared space."""
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA

    if "X_mofa" not in adata.obsm:
        return None
    z = np.asarray(adata.obsm["X_mofa"], dtype=np.float64)
    if "X_umap" in adata.obsm:
        xy = np.asarray(adata.obsm["X_umap"])
        title = "Joint UMAP (MuData) colored by label"
    else:
        xy = PCA(n_components=2, random_state=0).fit_transform(z)
        title = "PCA of MOFA factors"

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4))
    vals = adata.obs[color] if color in adata.obs.columns else None
    _scatter_embedding(axes[0], xy, vals, f"{title}\n({color})")
    axes[1].scatter(z[:, 0], z[:, 1], s=3, c="#54A24B", linewidths=0, alpha=0.5)
    axes[1].set_xlabel("MOFA Factor1 (RNA-private)")
    axes[1].set_ylabel("MOFA Factor2 (RNA-private)")
    axes[1].set_title("Factor1 vs Factor2\n(modality-private, not shared)")
    sc = axes[2].scatter(xy[:, 0], xy[:, 1], c=z[:, 0], s=3, cmap="viridis", linewidths=0)
    axes[2].set_title("Embedding\n(colored by MOFA Factor1)")
    axes[2].set_xlabel("Dim1")
    axes[2].set_ylabel("Dim2")
    axes[2].set_xticks([])
    axes[2].set_yticks([])
    fig.colorbar(sc, ax=axes[2], fraction=0.046, pad=0.04, label="Factor1")
    donor = donor_label or "configured donor"
    fig.suptitle(
        f"{donor} MOFA factors (stacked modality-private spaces -- see variance figure)",
        y=1.02,
        fontsize=11,
    )
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)


def plot_mofa_variance(
    variance_total: pd.DataFrame,
    variance_per_factor: pd.DataFrame,
    path: Path | str,
    n_factors: int = 14,
    dpi: int = 150,
    *,
    r2_recompute: dict | None = None,
) -> Path | None:
    """
    Left: total R^2 by view (with recomputation note if provided).
    Right: grouped bars -- both views per factor on a shared axis, truncated at Factor 14.
    """
    import matplotlib.pyplot as plt

    if variance_total is None or variance_total.empty:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    colors = ["#4C78A8", "#F58518"][: len(variance_total)]
    axes[0].bar(variance_total["view"], variance_total["r2_total"], color=colors)
    axes[0].set_ylabel("Total variance explained (percent)")
    axes[0].set_title("MOFA variance explained by view")
    ymax = float(variance_total["r2_total"].max()) * 1.25
    axes[0].set_ylim(0, ymax)
    if r2_recompute and r2_recompute.get("ok"):
        note = (
            f"Recomputed RNA R^2={r2_recompute.get('rna_r2_recomputed'):.4f} "
            f"(stored={r2_recompute.get('rna_r2_stored'):.4f}). "
            f"{r2_recompute.get('label', '')}"
        )
    else:
        note = (
            "Values are percent, not fractions: MOFA explains about 1 percent of RNA "
            "variance here. See the recomputation in run_params."
        )
    if note:
        axes[0].text(
            0.5,
            0.97,
            textwrap.fill(note, 62),
            transform=axes[0].transAxes,
            ha="center",
            va="top",
            fontsize=6.5,
            color="#444444",
        )

    sub = variance_per_factor.iloc[: int(n_factors)].copy()
    x = np.arange(len(sub))
    width = 0.38
    view_cols = list(sub.columns)
    for i, col in enumerate(view_cols):
        axes[1].bar(x + (i - 0.5) * width, sub[col].to_numpy(), width, label=str(col))
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(sub.index, rotation=45, ha="right")
    # Same percent scale as the totals panel: 0.49 here means 0.49 percent.
    axes[1].set_ylabel("Variance explained per factor (percent)")
    axes[1].set_title(
        f"Factors 1-{n_factors}: both views (shared axis)\n"
        f"Factors {n_factors + 1}-30 carry ~0 variance (dropped)"
    )
    axes[1].legend(frameon=False, fontsize=8)
    # annotate max non-dominant (units follow scale_branch from R^2 recompute)
    if len(view_cols) >= 2:
        arr = sub.to_numpy(dtype=float)
        nondom = []
        for row in arr:
            order = np.argsort(row)
            nondom.append(float(row[order[-2]]) if len(row) > 1 else 0.0)
        max_nd = max(nondom) if nondom else 0.0
        branch = (r2_recompute or {}).get("scale_branch")
        if branch == "both":
            stored = (r2_recompute or {}).get("rna_r2_stored")
            branch = "percent" if stored is not None and float(stored) > 1.0 else "direct"
        if branch == "percent":
            ann = (
                f"Max non-dominant R^2 = {max_nd:.5f}%\n"
                "(percent-scale HDF5 -- not a joint space)"
            )
        elif branch == "direct":
            ann = (
                f"Max non-dominant R^2 = {max_nd:.5f}\n"
                f"({100.0 * max_nd:.4f}% -- not a joint space)"
            )
        else:
            ann = (
                f"Max non-dominant R^2 = {max_nd:.5f}\n"
                "(units: see scale_branch -- not a joint space)"
            )
        axes[1].text(
            0.98,
            0.98,
            ann,
            transform=axes[1].transAxes,
            ha="right",
            va="top",
            fontsize=7.5,
            color="#333333",
        )
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)


def plot_rna_atac_scatter(
    bridge_df: pd.DataFrame,
    path: Path | str,
    focus_df: pd.DataFrame | None = None,
    dpi: int = 150,
    *,
    donor_label: str | None = None,
) -> Path | None:
    import matplotlib.pyplot as plt

    if bridge_df is None or bridge_df.empty:
        return None
    if "rna_std" in bridge_df.columns and "atac_std" in bridge_df.columns:
        x = np.log1p(bridge_df["rna_std"].astype(float))
        y = np.log1p(bridge_df["atac_std"].astype(float))
        xlabel = "RNA feature SD (log1p, MOFA inputs)"
        ylabel = "ATAC feature SD (log1p, MOFA inputs)"
        metric_note = "Bridge metric: feature SD on centered MOFA inputs (not mean expression)"
    else:
        x = np.log1p(bridge_df["rna_mean_expression"].astype(float))
        y = np.log1p(bridge_df["atac_mean_gene_activity"].astype(float))
        xlabel = "RNA mean (log1p)"
        ylabel = "ATAC gene activity mean (log1p)"
        metric_note = "Bridge metric: mean values"
    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    ax.scatter(x, y, s=4, alpha=0.25, c="#4C78A8", linewidths=0)
    if focus_df is not None and not focus_df.empty:
        focus_syms = set(focus_df["gene_symbol"].astype(str).str.upper())
        sub = bridge_df[bridge_df["gene_symbol"].astype(str).str.upper().isin(focus_syms)]
        if not sub.empty and "rna_std" in sub.columns:
            ax.scatter(
                np.log1p(sub["rna_std"].astype(float)),
                np.log1p(sub["atac_std"].astype(float)),
                s=36,
                c="#F58518",
                edgecolors="k",
                linewidths=0.4,
                zorder=3,
            )
            for _, row in sub.iterrows():
                sym = str(row.get("gene_symbol") or row["gene_id_no_version"])
                ax.annotate(
                    sym,
                    (np.log1p(float(row["rna_std"])), np.log1p(float(row["atac_std"]))),
                    fontsize=7,
                    xytext=(4, 4),
                    textcoords="offset points",
                )
    r = bridge_df.attrs.get("pearson_log1p")
    rho = bridge_df.attrs.get("spearman_log1p")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    title = f"{donor_label or 'configured donor'} RNA<->ATAC gene bridge (supporting)"
    if r == r:
        title += f"\nPearson r={r:.3f}, Spearman rho={rho:.3f}"
    ax.set_title(title)
    ax.text(
        0.02,
        0.02,
        metric_note + "\nWeak genome-wide concordance is expected",
        transform=ax.transAxes,
        fontsize=7,
        va="bottom",
        color="#444444",
    )
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)


def plot_focus_gene_bars(
    focus_df: pd.DataFrame,
    path: Path | str,
    dpi: int = 150,
) -> Path | None:
    """Focus-marker panel with dual scales (concordance, not equal expression).

    Returns None on the MOFA interpretation path where values are centered residuals
    (not expression means), or when ``means_are_centered`` is set / negatives appear.
    """
    import matplotlib.pyplot as plt

    if focus_df is None or focus_df.empty:
        return None
    if bool(getattr(focus_df, "attrs", {}).get("means_are_centered", False)):
        return None
    if "rna_mean_expression" not in focus_df.columns:
        return None
    df = focus_df.copy()
    rna = df["rna_mean_expression"].astype(float).to_numpy()
    atac = df["atac_mean_gene_activity"].astype(float).to_numpy()
    # Centered MOFA inputs: near-zero mean and/or negatives are not expression
    if (np.nanmean(np.abs(rna)) < 1e-3 and np.nanmin(rna) < 0) or np.nanmin(rna) < -1e-6:
        return None
    if np.nanmin(atac) < -1e-6:
        return None
    labels = df["gene_symbol"].replace("", pd.NA).fillna(df["gene_id_no_version"]).astype(str)
    x = np.arange(len(df))
    width = 0.38

    fig, axes = plt.subplots(1, 2, figsize=(max(10, 0.85 * len(df) + 4), 4.0))

    ax = axes[0]
    ax.bar(x - width / 2, np.log1p(rna), width, label="RNA (log1p)", color="#4C78A8")
    ax.set_ylabel("RNA log1p mean", color="#4C78A8")
    ax.tick_params(axis="y", labelcolor="#4C78A8")
    ax2 = ax.twinx()
    ax2.bar(x + width / 2, np.log1p(atac), width, label="ATAC activity (log1p)", color="#F58518")
    ax2.set_ylabel("ATAC log1p mean", color="#F58518")
    ax2.tick_params(axis="y", labelcolor="#F58518")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_title("Focus means (different units)")
    lines1, labs1 = ax.get_legend_handles_labels()
    lines2, labs2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labs1 + labs2, frameon=False, fontsize=8, loc="upper right")

    ax = axes[1]
    rna_z = (rna - rna.mean()) / (rna.std() if rna.std() > 0 else 1.0)
    atac_z = (atac - atac.mean()) / (atac.std() if atac.std() > 0 else 1.0)
    ax.bar(x - width / 2, rna_z, width, label="RNA z", color="#4C78A8")
    ax.bar(x + width / 2, atac_z, width, label="ATAC z", color="#F58518")
    ax.axhline(0, color="#666666", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Within-panel z-score")
    ax.set_title("Concordance (z-scored means)")
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle(
        "Focus markers: co-detection / concordance (not equal expression)",
        fontsize=11,
        y=1.02,
    )
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)


def _correlation_title_stats(corr_df: pd.DataFrame, noise_delta: float = 0.01) -> str:
    """Build an honest title from global vs within-type rows."""
    if corr_df is None or corr_df.empty or "scope" not in corr_df.columns:
        return "Focus-marker RNA<->ATAC correlations"
    g = corr_df[corr_df["scope"] == "global"].copy()
    c = corr_df[corr_df["scope"] == "within_cell_type"].copy()
    if g.empty:
        return "Focus-marker RNA<->ATAC correlations"
    if c.empty:
        return "Focus-marker RNA<->ATAC correlations (global scope only)"
    gmap = {str(s).upper(): float(r) for s, r in zip(g["gene_symbol"], g["pearson_r"])}
    higher = 0
    meaningfully_higher = 0
    comparable = 0
    for _, row in c.iterrows():
        sym = str(row["gene_symbol"]).upper()
        if sym not in gmap:
            continue
        gv, cv = gmap[sym], float(row["pearson_r"])
        if not (cv == cv) or not (gv == gv):  # nan
            continue
        comparable += 1
        if cv > gv:
            higher += 1
            if (cv - gv) > float(noise_delta):
                meaningfully_higher += 1
    return (
        f"Focus correlations: within-type exceeds global in {higher} of {comparable} "
        f"({meaningfully_higher} by >{noise_delta:g}); "
        "restricting usually shrinks r"
    )


def plot_focus_correlations(
    corr_df: pd.DataFrame,
    path: Path | str,
    dpi: int = 150,
    *,
    noise_delta: float = 0.01,
) -> Path | None:
    """Compare global vs within-cell-type Pearson r; title matches the bars."""
    import matplotlib.pyplot as plt

    if corr_df is None or corr_df.empty or "scope" not in corr_df.columns:
        return None
    global_df = corr_df[corr_df["scope"] == "global"].copy()
    cell_df = corr_df[corr_df["scope"] == "within_cell_type"].copy()
    if global_df.empty:
        return None

    genes = global_df["gene_symbol"].astype(str).tolist()
    gmap = {g.upper(): float(r) for g, r in zip(global_df["gene_symbol"], global_df["pearson_r"])}
    cmap = {}
    notes = {}
    has_within = not cell_df.empty
    if has_within:
        cmap = {g.upper(): float(r) for g, r in zip(cell_df["gene_symbol"], cell_df["pearson_r"])}
        if "note" in cell_df.columns:
            notes = {
                g.upper(): str(n)
                for g, n in zip(cell_df["gene_symbol"], cell_df["note"])
                if str(n) not in {"", "nan", "None"}
            }

    x = np.arange(len(genes))
    width = 0.38
    g_vals = [gmap.get(g.upper(), np.nan) for g in genes]
    c_vals = [cmap.get(g.upper(), np.nan) for g in genes] if has_within else None

    fig, ax = plt.subplots(figsize=(max(7.5, 0.75 * len(genes)), 4.4))
    if has_within:
        ax.bar(x - width / 2, g_vals, width, label="Global (all nuclei)", color="#9E9E9E")
        ax.bar(x + width / 2, c_vals, width, label="Within expected cell type(s)", color="#54A24B")
    else:
        # Global-only MOFA path: do not invent a within-cell-type legend series
        ax.bar(x, g_vals, width * 1.5, label="Global (all nuclei)", color="#9E9E9E")
    ax.axhline(0, color="#666666", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=45, ha="right")
    ax.set_ylabel("Pearson r (RNA vs ATAC activity)")
    ax.set_title(_correlation_title_stats(corr_df, noise_delta=noise_delta))
    # Flag low-n / zero-RNA notes under bars
    ymin = ax.get_ylim()[0]
    for i, g in enumerate(genes):
        note = notes.get(g.upper(), "")
        if note:
            ax.text(i, ymin, "!", ha="center", va="bottom", fontsize=9, color="#B22222")
    ax.legend(frameon=False, fontsize=8)
    ax.text(
        0.01,
        0.01,
        "! = flagged (n below threshold, rna_mean=0, or undefined r). "
        "Near-zero r everywhere is the main pattern.",
        transform=ax.transAxes,
        fontsize=7,
        va="bottom",
        color="#444444",
    )
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)


def plot_qc_association(
    qc_df: pd.DataFrame,
    path: Path | str,
    dpi: int = 150,
) -> Path | None:
    """Grouped bars: eta squared by QC metric, one series per clustering (R2)."""
    import matplotlib.pyplot as plt

    if qc_df is None or qc_df.empty:
        return None
    metrics = list(qc_df["qc_metric"].astype(str).unique())
    labels = list(qc_df["labelling"].astype(str).unique())
    x = np.arange(len(metrics))
    width = 0.8 / max(len(labels), 1)
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    for i, lab in enumerate(labels):
        sub = qc_df[qc_df["labelling"].astype(str) == lab]
        vals = [
            float(sub.loc[sub["qc_metric"] == m, "eta_squared"].iloc[0])
            if m in set(sub["qc_metric"].astype(str))
            else np.nan
            for m in metrics
        ]
        ax.bar(x + (i - (len(labels) - 1) / 2) * width, vals, width, label=str(lab))
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, rotation=20, ha="right")
    ax.set_ylabel("Eta squared")
    ax.set_title("ATAC QC vs clustering: depth explains ATAC Leiden")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)


def plot_factor_label_heatmap(
    assoc_df: pd.DataFrame,
    path: Path | str,
    dpi: int = 150,
    n_factors: int = 18,
) -> Path | None:
    """Factors vs labellings, eta squared as fill (R3)."""
    import matplotlib.pyplot as plt

    if assoc_df is None or assoc_df.empty:
        return None
    sub = assoc_df.copy()
    factor_order = [f"Factor{i}" for i in range(1, int(n_factors) + 1)]
    present = [f for f in factor_order if f in set(sub["factor"].astype(str))]
    if not present:
        present = list(sub["factor"].astype(str).unique())
    pivot = (
        sub.pivot_table(index="factor", columns="labelling", values="eta_squared", aggfunc="mean")
        .reindex(present)
    )
    fig, ax = plt.subplots(figsize=(7.2, max(4.5, 0.28 * len(present) + 1.5)))
    im = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="viridis", vmin=0, vmax=0.7)
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=30, ha="right")
    ax.set_title("MOFA factor scores vs cell identity (eta squared)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="eta squared")
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)


def plot_wnn_umap(
    cache_df: pd.DataFrame,
    path: Path | str,
    dpi: int = 150,
    *,
    donor_label: str | None = None,
    atac_leiden: pd.Series | None = None,
) -> Path | None:
    """WNN UMAP colored by Azimuth label and by ATAC Leiden (R4)."""
    import matplotlib.pyplot as plt

    if cache_df is None or cache_df.empty:
        return None
    if "wnn_umap_1" not in cache_df.columns or "wnn_umap_2" not in cache_df.columns:
        return None
    xy = cache_df[["wnn_umap_1", "wnn_umap_2"]].to_numpy(dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    azimuth = cache_df["azimuth_label"] if "azimuth_label" in cache_df.columns else None
    _scatter_embedding(axes[0], xy, azimuth, "WNN UMAP (Azimuth label)")
    axes[0].set_xlabel("WNN UMAP1")
    axes[0].set_ylabel("WNN UMAP2")
    atac = atac_leiden
    if atac is None and "atac_leiden" in cache_df.columns:
        atac = cache_df["atac_leiden"]
    _scatter_embedding(axes[1], xy, atac, "WNN UMAP (ATAC Leiden)")
    axes[1].set_xlabel("WNN UMAP1")
    axes[1].set_ylabel("WNN UMAP2")
    donor = donor_label or "configured donor"
    fig.suptitle(
        f"{donor} WNN embedding follows RNA identity more than ATAC",
        y=1.02,
        fontsize=11,
    )
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)


def plot_cluster_concordance(
    frame: pd.DataFrame,
    path: Path | str,
    dpi: int = 150,
    *,
    row_col: str = "rna_leiden",
    col_col: str = "atac_leiden",
) -> Path | None:
    """RNA Leiden vs ATAC Leiden contingency, row-normalised (R1)."""
    import matplotlib.pyplot as plt

    if frame is None or frame.empty:
        return None
    if row_col not in frame.columns or col_col not in frame.columns:
        return None
    ct = pd.crosstab(frame[row_col].astype(str), frame[col_col].astype(str), normalize="index")
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    im = ax.imshow(ct.to_numpy(dtype=float), aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_yticks(np.arange(len(ct.index)))
    ax.set_yticklabels(ct.index, fontsize=7)
    ax.set_xticks(np.arange(len(ct.columns)))
    ax.set_xticklabels(ct.columns, fontsize=7, rotation=90)
    ax.set_ylabel("RNA Leiden")
    ax.set_xlabel("ATAC Leiden")
    ax.set_title("RNA Leiden vs ATAC Leiden (row-normalised)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="row fraction")
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)
