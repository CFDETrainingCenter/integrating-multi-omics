"""GTEx GCT loading and gene ID harmonization helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from scripts.common.paths import module_root


def strip_ensembl_version(ids: pd.Series | pd.Index) -> pd.Series:
    s = pd.Series(ids.astype(str))
    return s.str.replace(r"\.\d+$", "", regex=True)


def read_gct_gz(path: Path | str) -> pd.DataFrame:
    """
    Read a GCT/GCT.gz into a DataFrame indexed by versionless Ensembl ID.

    Returns columns: gene_symbol + sample TPM/count columns.
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t", compression="gzip", skiprows=2)
    if "Name" not in df.columns or "Description" not in df.columns:
        raise ValueError(f"Unexpected GCT columns in {path}: {df.columns[:5].tolist()}")
    df = df.rename(columns={"Name": "gene_id", "Description": "gene_symbol"})
    df["gene_id_no_version"] = strip_ensembl_version(df["gene_id"])
    df = df.drop_duplicates(subset=["gene_id_no_version"], keep="first")
    df = df.set_index("gene_id_no_version")
    return df


def _gtex_cfg(cfg: dict[str, Any]) -> dict[str, Any]:
    """Resolve GTEx paths. Prefer Module 4 (authoritative); Module 2 may omit gtex."""
    if cfg.get("module2", {}).get("gtex"):
        return dict(cfg["module2"]["gtex"])
    m4 = cfg.get("module4") or {}
    tpm = (m4.get("gtex") or {}).get("tpm_file") or (m4.get("data") or {}).get("gtex_tpm")
    if not tpm:
        raise KeyError("GTEx tpm_file not found under module4.data / module4.gtex")
    return {"tpm_file": tpm, "max_subjects": (m4.get("gtex") or {}).get("max_subjects", 100)}


def load_gtex_lung_tpm(cfg: dict[str, Any]) -> pd.DataFrame:
    rel = _gtex_cfg(cfg)["tpm_file"]
    path = module_root(cfg) / rel
    return read_gct_gz(path)


def gtex_subject_means(gtex: pd.DataFrame, cfg: dict[str, Any]) -> pd.Series:
    """Mean TPM across (optionally subsampled) subjects for each gene."""
    sample_cols = [c for c in gtex.columns if c not in {"gene_id", "gene_symbol"}]
    gcfg = _gtex_cfg(cfg)
    max_n = gcfg.get("max_subjects")
    if max_n is not None and len(sample_cols) > int(max_n):
        rng = np.random.default_rng(int(cfg.get("module2", {}).get("random_seed", 0)))
        sample_cols = list(rng.choice(sample_cols, size=int(max_n), replace=False))
    means = gtex[sample_cols].astype(float).mean(axis=1)
    means.name = "gtex_lung_mean_tpm"
    return means


def shared_gene_table(hubmap_genes: list[str], gtex: pd.DataFrame) -> pd.DataFrame:
    hub = pd.Index(hubmap_genes.astype(str) if hasattr(hubmap_genes, "astype") else list(map(str, hubmap_genes)))
    hub = strip_ensembl_version(hub)
    shared = sorted(set(hub) & set(gtex.index.astype(str)))
    out = pd.DataFrame({"gene_id_no_version": shared})
    out["gene_symbol"] = gtex.loc[shared, "gene_symbol"].astype(str).values
    out["in_hubmap"] = True
    out["in_gtex"] = True
    return out
