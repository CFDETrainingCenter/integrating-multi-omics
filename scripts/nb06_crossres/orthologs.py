"""One-to-one human-mouse ortholog mapping with loss accounting."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from scripts.common.paths import module_root


def load_ortholog_table(cfg: dict[str, Any]) -> pd.DataFrame:
    rel = cfg["module4"]["data"]["orthologs"]
    path = module_root(cfg) / rel
    # skip comment header lines
    df = pd.read_csv(path, sep="\t", comment="#")
    need = {"human_symbol", "mouse_symbol"}
    if not need.issubset(df.columns):
        raise ValueError(f"Ortholog table missing columns {need}: {path}")
    df["human_symbol"] = df["human_symbol"].astype(str).str.upper()
    # mouse symbols stay title-case from Ensembl (Sftpc); also store upper for joins
    df["mouse_symbol"] = df["mouse_symbol"].astype(str)
    df["mouse_symbol_upper"] = df["mouse_symbol"].str.upper()
    return df


def map_human_to_mouse(symbols: list[str], ortho: pd.DataFrame) -> pd.DataFrame:
    panel = pd.DataFrame({"human_symbol": [s.upper() for s in symbols]})
    return panel.merge(ortho[["human_symbol", "mouse_symbol"]], on="human_symbol", how="left")


def ortholog_loss_record(
    human_genes: list[str] | pd.Index,
    mouse_genes: list[str] | pd.Index,
    ortho: pd.DataFrame,
) -> pd.DataFrame:
    """
    Report intersection loss -- every ortholog-using step must surface this.
    """
    h = {str(g).upper() for g in human_genes if str(g).strip()}
    # mouse genes may be title-case or upper
    m_raw = [str(g) for g in mouse_genes if str(g).strip()]
    m_upper = {g.upper() for g in m_raw}
    o_h = set(ortho["human_symbol"].astype(str).str.upper())
    o_m = set(ortho["mouse_symbol_upper"].astype(str))
    # one2one pairs where both sides present in the contrast gene universes
    pairs = ortho[
        ortho["human_symbol"].isin(h) & ortho["mouse_symbol_upper"].isin(m_upper)
    ]
    n_h, n_m, n_o, n_r = len(h), len(m_upper), len(ortho), len(pairs)
    rows = [
        {"metric": "genes_in_human_contrast", "n": n_h},
        {"metric": "genes_in_mouse_contrast", "n": n_m},
        {"metric": "one2one_orthologs_available", "n": n_o},
        {"metric": "retained_after_intersection", "n": n_r},
        {"metric": "dropped_human_side", "n": n_h - n_r},
        {"metric": "dropped_mouse_side", "n": n_m - n_r},
    ]
    return pd.DataFrame(rows)


def align_logfc_by_ortholog(
    human_logfc: pd.Series,
    mouse_logfc: pd.Series,
    ortho: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Return (aligned dataframe with human_symbol, mouse_symbol, logFC_h, logFC_m,
    loss_record).
    """
    h = human_logfc.copy()
    h.index = h.index.astype(str).str.upper()
    m = mouse_logfc.copy()
    # index mouse by upper for join, keep display symbol
    m_df = m.rename("logFC_mouse").to_frame()
    m_df["mouse_symbol_upper"] = m_df.index.astype(str).str.upper()
    m_df["mouse_symbol"] = m_df.index.astype(str)

    loss = ortholog_loss_record(h.index, m.index, ortho)
    pairs = ortho[
        ortho["human_symbol"].isin(h.index) & ortho["mouse_symbol_upper"].isin(m_df["mouse_symbol_upper"])
    ][["human_symbol", "mouse_symbol", "mouse_symbol_upper"]].drop_duplicates()
    aligned = pairs.merge(h.rename("logFC_human"), left_on="human_symbol", right_index=True, how="inner")
    aligned = aligned.merge(
        m_df[["mouse_symbol_upper", "logFC_mouse"]],
        on="mouse_symbol_upper",
        how="inner",
    )
    return aligned.drop(columns=["mouse_symbol_upper"]), loss
