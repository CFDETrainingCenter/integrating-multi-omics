"""HuBMAP-only pseudobulk summaries and teaching heatmaps (Module 2).

GTEx comparison lives in Module 4 -- do not read GTEx here.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd


def mean_expression_by_label(
    adata,
    label_col: str = "azimuth_label",
    *,
    linear: bool = False,
) -> pd.DataFrame:
    """
    Mean expression per label (columns=labels, index=genes).

    If ``linear=True``, recover a linear scale with ``expm1`` assuming ``adata.X``
    is log1p(CP10K). Otherwise use ``adata.X`` as-is (typically log1p).

    Stays sparse end-to-end: ``expm1(0)=0`` so zeros never densify. Peak memory is
    sparse float64 ``.data`` (~0.4 GB here) plus one dense 1xn_genes mean vector
    per label (~0.3 MB), not a full dense 18kx42k matrix.
    """
    import scipy.sparse as sp

    if label_col not in adata.obs.columns:
        raise KeyError(f"{label_col} not in adata.obs")
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    else:
        X = X.tocsr()
    # Cast on sparse .data only (0.22 GB float32 -> ~0.44 GB float64), never densify.
    X = X.astype(np.float64, copy=True)
    if linear:
        X = X.expm1()
    labels = adata.obs[label_col].astype(str).fillna("NA")
    genes = adata.var_names.astype(str)
    out = {}
    for lab in sorted(labels.unique()):
        idx = np.where(labels.values == lab)[0]
        if len(idx) == 0:
            continue
        # Sparse mean includes implicit zeros; same arithmetic as the dense path.
        out[lab] = np.asarray(X[idx].mean(axis=0)).ravel()
    pb = pd.DataFrame(out, index=genes)
    pb.index.name = "gene_id_no_version"
    return pb


# Backward-compatible name used by older call sites
def pseudobulk_by_label(adata, label_col: str = "azimuth_label") -> pd.DataFrame:
    return mean_expression_by_label(adata, label_col=label_col, linear=False)


def label_fractions(adata, label_col: str = "azimuth_label") -> pd.Series:
    labs = adata.obs[label_col].astype(str)
    return labs.value_counts(normalize=True).sort_index()


def composition_weighted_pseudobulk(
    adata,
    label_col: str = "azimuth_label",
) -> tuple[pd.Series, pd.Series, pd.DataFrame]:
    """
    Composition-weighted bulk reconstruction on the **linear** scale.

    For each gene: sum_t (fraction_t * mean_linear_expression_t).
    Returns (weighted_vector, fractions, per_label_linear_means).
    """
    means = mean_expression_by_label(adata, label_col=label_col, linear=True)
    fracs = label_fractions(adata, label_col=label_col)
    # align
    shared = [c for c in means.columns if c in fracs.index]
    means = means[shared]
    fracs = fracs.loc[shared]
    weighted = means.mul(fracs, axis=1).sum(axis=1)
    weighted.name = "composition_weighted_linear"
    return weighted, fracs, means


def max_across_labels_pseudobulk(
    adata,
    label_col: str = "azimuth_label",
    *,
    linear: bool = True,
) -> pd.Series:
    """Max across cell-type means (contrast to composition-weighted; not a bulk average)."""
    means = mean_expression_by_label(adata, label_col=label_col, linear=linear)
    s = means.max(axis=1)
    s.name = "max_across_labels_linear" if linear else "max_across_labels"
    return s


def hubmap_marker_matrix(
    adata,
    markers: Sequence[str],
    label_col: str = "azimuth_label",
    *,
    top_n_labels: int = 12,
) -> pd.DataFrame:
    """
    Marker x cell-type mean matrix (log1p X) restricted to the most abundant labels.
    Index = gene symbol when resolvable via ``hugo_symbol``, else gene id.
    """
    counts = adata.obs[label_col].astype(str).value_counts()
    keep = list(counts.head(int(top_n_labels)).index)
    pb = mean_expression_by_label(adata, label_col=label_col, linear=False)[keep]

    # Resolve markers to var index
    hugo = None
    if "hugo_symbol" in adata.var.columns:
        hugo = (
            adata.var["hugo_symbol"]
            .astype(str)
            .str.upper()
            .replace({"NAN": pd.NA, "NONE": pd.NA, "": pd.NA})
        )
    rows = []
    for sym in markers:
        sym_u = str(sym).upper()
        gid = None
        if hugo is not None:
            hits = hugo[hugo == sym_u]
            if len(hits):
                gid = hits.index[0]
        if gid is None and sym in adata.var_names:
            gid = sym
        if gid is None or gid not in pb.index:
            continue
        row = {"gene_symbol": sym_u, "gene_id_no_version": gid}
        for lab in keep:
            row[lab] = float(pb.loc[gid, lab])
        rows.append(row)
    return pd.DataFrame(rows)


def plot_hubmap_marker_heatmap(
    marker_df: pd.DataFrame,
    path: Path,
    *,
    dpi: int = 150,
    title: str = "HuBMAP markers (row z-score within resource)",
) -> Path:
    """Legible heatmap: row-z-score across cell types inside HuBMAP only."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    from scripts.common.plotting import apply_figure_style

    apply_figure_style()
    if marker_df.empty or "gene_symbol" not in marker_df.columns:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, "No markers resolved", ha="center")
        ax.axis("off")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return path

    meta = {"gene_symbol", "gene_id_no_version"}
    cols = [c for c in marker_df.columns if c not in meta]
    mat = marker_df.set_index("gene_symbol")[cols].astype(float)
    mat_z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, np.nan), axis=0)
    # Readable sizing: ~0.55 in per row/col
    fig_w = max(7.0, 0.55 * len(cols) + 2.0)
    fig_h = max(4.0, 0.55 * len(mat_z) + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(mat_z.fillna(0), cmap="vlag", center=0, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Cell type (top by abundance)")
    ax.set_ylabel("Gene")
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    fig.tight_layout()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


# --- Deprecated GTEx helpers kept as stubs so old imports fail clearly ---
def marker_comparison_table(*_a, **_k):
    raise RuntimeError("GTEx left Module 2; use Module 4 cross-ecosystem comparison.")


def hubmap_vs_gtex_scatter_df(*_a, **_k):
    raise RuntimeError("GTEx left Module 2; use Module 4 cross-ecosystem comparison.")


def plot_marker_heatmap(*_a, **_k):
    raise RuntimeError("Use plot_hubmap_marker_heatmap (HuBMAP-only).")


def plot_hubmap_gtex_scatter(*_a, **_k):
    raise RuntimeError("GTEx left Module 2; use Module 4 cross-ecosystem comparison.")
