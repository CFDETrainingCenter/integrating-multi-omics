"""Small plotting helpers that always write to disk."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt

_STYLE_APPLIED = False


def apply_figure_style(cfg=None) -> None:
    """Load config/figure_style.mplstyle if present."""
    from scripts.common.paths import MODULE_ROOT

    style = MODULE_ROOT / "config" / "figure_style.mplstyle"
    if style.exists():
        plt.style.use(str(style))


def _ensure_figure_style() -> None:
    global _STYLE_APPLIED
    if not _STYLE_APPLIED:
        apply_figure_style()
        _STYLE_APPLIED = True


def save_figure(path: Path | str, dpi: int = 150, close: bool = True) -> Path:
    """Write the current matplotlib figure to disk; display inline under IPython/Jupyter."""
    _ensure_figure_style()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    try:
        from IPython import get_ipython
        from IPython.display import Image, display

        if get_ipython() is not None:
            display(Image(filename=str(path)))
    except Exception:  # noqa: BLE001 -- never fail a write because display is unavailable
        pass
    if close:
        plt.close()
    return path


def qc_violin_panel(adata, keys: Sequence[str], path: Path | str, dpi: int = 150) -> Path:
    import scanpy as sc

    n = len(keys)
    fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 3.5))
    if n == 1:
        axes = [axes]
    for ax, key in zip(axes, keys):
        sc.pl.violin(adata, keys=key, ax=ax, show=False, stripplot=False)
        ax.set_title(key)
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)


def qc_scatter(
    adata,
    x: str,
    y: str,
    path: Path | str,
    color: str | None = None,
    dpi: int = 150,
) -> Path:
    import scanpy as sc

    sc.pl.scatter(adata, x=x, y=y, color=color, show=False)
    return save_figure(path, dpi=dpi, close=True)


def umap_panel(
    adata,
    colors: Sequence[str],
    path: Path | str,
    dpi: int = 150,
    ncols: int | None = None,
    wspace: float = 0.35,
    legend_loc: str = "right margin",
    legend_fontsize: int = 7,
    size: float | None = None,
) -> Path:
    import scanpy as sc

    color_list = list(colors)
    n = len(color_list)
    cols = ncols if ncols is not None else min(2, n) if n >= 4 else min(3, max(n, 1))
    sc.pl.umap(
        adata,
        color=color_list,
        ncols=cols,
        wspace=wspace,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize,
        size=size,
        show=False,
    )
    return save_figure(path, dpi=dpi, close=True)
