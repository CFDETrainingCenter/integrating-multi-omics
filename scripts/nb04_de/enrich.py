"""Pathway enrichment: local GMT (default) or optional Enrichr online."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import pandas as pd

from scripts.common.paths import module_root, resolve


def _geneset_dir(cfg: dict[str, Any]) -> Path:
    enr = cfg.get("module2", {}).get("enrichment") or {}
    rel = enr.get("geneset_dir") or cfg.get("paths", {}).get("genesets") or "data/genesets"
    root = module_root(cfg)
    p = Path(rel)
    return p if p.is_absolute() else root / p


def resolve_gmt_paths(cfg: dict[str, Any], gene_sets: list[str] | None = None) -> list[Path]:
    """Map library names to local ``*.gmt`` files under data/genesets/."""
    enr = cfg.get("module2", {}).get("enrichment") or {}
    names = list(gene_sets if gene_sets is not None else (enr.get("gene_sets") or []))
    # Default teaching set (no KEGG)
    if not names:
        names = ["MSigDB_Hallmark_2020", "Reactome_2022", "GO_Biological_Process_2023"]
    gdir = _geneset_dir(cfg)
    paths: list[Path] = []
    missing: list[str] = []
    for name in names:
        stem = name if name.endswith(".gmt") else f"{name}.gmt"
        path = gdir / stem
        if not path.exists():
            missing.append(str(path))
        else:
            paths.append(path)
    if missing:
        raise FileNotFoundError(
            "Missing local GMT file(s) required for offline enrichment:\n  - "
            + "\n  - ".join(missing)
            + "\nPlace Hallmark / Reactome_2022 / GO_Biological_Process_2023 under data/genesets/."
        )
    return paths


def _is_retryable_enrichr_error(exc: BaseException) -> bool:
    msg = str(exc).lower()
    name = exc.__class__.__name__.lower()
    needles = ("429", "rate", "timeout", "timed out", "connection", "temporarily", "503", "502")
    return any(n in msg for n in needles) or "api" in name


def _enrichr_once(gene_list: list[str], gene_sets: list[str], organism: str, cutoff: float):
    import gseapy as gp

    return gp.enrichr(
        gene_list=gene_list,
        gene_sets=gene_sets,
        organism=organism,
        outdir=None,
        cutoff=cutoff,
        no_plot=True,
        verbose=False,
    )


def _enrich_local_once(
    gene_list: list[str],
    gmt_paths: list[Path],
    background: list[str],
    cutoff: float,
):
    import gseapy as gp

    return gp.enrich(
        gene_list=gene_list,
        gene_sets=[str(p) for p in gmt_paths],
        background=background,
        outdir=None,
        cutoff=cutoff,
        verbose=False,
    )


def run_enrichment_for_groups(
    filtered_markers: pd.DataFrame,
    cfg: dict[str, Any],
    background_genes: list[str] | None = None,
) -> pd.DataFrame:
    """
    Over-representation enrichment per group.

    Default ``enrichment.mode=local_gmt``: gseapy.enrich against local GMTs with
    ``background=`` equal to genes measured in the object (no network).

    ``enrichment.mode=enrichr_online``: legacy Enrichr HTTP path (off by default).
    """
    enr_cfg = cfg["module2"]["enrichment"]
    if not enr_cfg.get("enabled", True):
        return pd.DataFrame()
    if filtered_markers.empty:
        return pd.DataFrame()

    mode = str(enr_cfg.get("mode", "local_gmt")).lower()
    gene_sets = list(
        enr_cfg.get("gene_sets")
        or ["MSigDB_Hallmark_2020", "Reactome_2022", "GO_Biological_Process_2023"]
    )
    # Hard exclude KEGG if someone leaves it in config
    gene_sets = [g for g in gene_sets if "kegg" not in g.lower()]
    organism = str(enr_cfg.get("organism", "human"))
    cutoff = float(enr_cfg.get("cutoff_padj", 0.05))

    gmt_paths: list[Path] = []
    if mode == "local_gmt":
        gmt_paths = resolve_gmt_paths(cfg, gene_sets)
        if not background_genes:
            raise ValueError("local_gmt mode requires background_genes (measured gene symbols)")
        background = sorted({str(g).strip() for g in background_genes if str(g).strip()})
    else:
        background = []

    retry_max = int(enr_cfg.get("retry_max", 4))
    backoff = float(enr_cfg.get("retry_backoff_sec", 5))
    pause = float(enr_cfg.get("pause_between_groups_sec", 2)) if mode != "local_gmt" else 0.0

    frames: list[pd.DataFrame] = []
    groups = list(filtered_markers.groupby("group", observed=True))
    for i, (group, sub) in enumerate(groups):
        genes = (
            sub["gene_symbol"]
            .astype(str)
            .str.strip()
            .replace("", pd.NA)
            .dropna()
            .drop_duplicates()
            .tolist()
        )
        if len(genes) < 5:
            frames.append(
                pd.DataFrame(
                    [
                        {
                            "group": group,
                            "Gene_set": ",".join(gene_sets),
                            "Term": "ENRICHMENT_SKIPPED: fewer than 5 input genes",
                            "Adjusted P-value": 1.0,
                            "Genes": ";".join(genes),
                            "n_input_genes": len(genes),
                            "mode": mode,
                        }
                    ]
                )
            )
            continue

        last_exc: BaseException | None = None
        res = None
        attempts = 1 if mode == "local_gmt" else retry_max
        for attempt in range(1, attempts + 1):
            try:
                if mode == "local_gmt":
                    enr = _enrich_local_once(genes, gmt_paths, background, cutoff)
                else:
                    enr = _enrichr_once(genes, gene_sets, organism, cutoff)
                res = enr.results.copy() if enr is not None else None
                last_exc = None
                break
            except Exception as exc:  # noqa: BLE001
                last_exc = exc
                if mode == "local_gmt" or attempt >= attempts or not _is_retryable_enrichr_error(exc):
                    break
                sleep_s = backoff * (2 ** (attempt - 1))
                print(
                    f"Enrichr retry {attempt}/{attempts} for group={group} "
                    f"after {exc.__class__.__name__}: sleeping {sleep_s:.1f}s"
                )
                time.sleep(sleep_s)

        if last_exc is not None:
            frames.append(
                pd.DataFrame(
                    [
                        {
                            "group": group,
                            "Gene_set": ",".join(gene_sets),
                            "Term": f"ENRICHMENT_FAILED: {last_exc.__class__.__name__}",
                            "Adjusted P-value": 1.0,
                            "Genes": ";".join(genes[:20]),
                            "n_input_genes": len(genes),
                            "error": str(last_exc),
                            "mode": mode,
                        }
                    ]
                )
            )
        elif res is not None and not res.empty:
            res = res.copy()
            res.insert(0, "group", group)
            res["n_input_genes"] = len(genes)
            res["mode"] = mode
            if mode == "local_gmt":
                res["n_background_genes"] = len(background)
            frames.append(res)

        if pause > 0 and i < len(groups) - 1:
            time.sleep(pause)

    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    padj_col = "Adjusted P-value" if "Adjusted P-value" in out.columns else None
    if padj_col is None and "FDR q-val" in out.columns:
        padj_col = "FDR q-val"
    if padj_col:
        out = out.sort_values(["group", padj_col], ascending=[True, True])
    return out


def run_prerank_gsea(
    ranked: pd.Series,
    cfg: dict[str, Any],
) -> pd.DataFrame:
    """
    Preranked GSEA against the same local GMTs (NES for Hallmark / Reactome / GO BP).

    ``ranked`` index = gene symbols, values = ranking metric (e.g. Wilcoxon score).
    """
    import gseapy as gp

    enr_cfg = cfg["module2"]["enrichment"]
    pr = enr_cfg.get("prerank") or {}
    if not pr.get("enabled", True):
        return pd.DataFrame()
    gene_sets = [
        g
        for g in (enr_cfg.get("gene_sets") or [])
        if "kegg" not in str(g).lower()
    ] or ["MSigDB_Hallmark_2020", "Reactome_2022", "GO_Biological_Process_2023"]
    gmt_paths = resolve_gmt_paths(cfg, gene_sets)
    rnk = ranked.dropna().astype(float)
    rnk = rnk[~rnk.index.duplicated(keep="first")]
    if rnk.empty:
        return pd.DataFrame()
    if int((rnk < 0).sum()) == 0:
        raise ValueError(
            "prerank ranked list has no negative scores; "
            "a one-sided Wilcoxon top-n slice is not a valid prerank input. "
            "Build ranks from a full-universe rank_genes_groups pass."
        )

    pre = gp.prerank(
        rnk=rnk,
        gene_sets=[str(p) for p in gmt_paths],
        min_size=int(pr.get("min_size", 15)),
        max_size=int(pr.get("max_size", 500)),
        permutation_num=int(pr.get("permutation_num", 100)),
        outdir=None,
        seed=int(cfg["module2"].get("random_seed", 0)),
        verbose=False,
    )
    res = pre.res2d.copy() if pre is not None else pd.DataFrame()
    if not res.empty:
        res["mode"] = "prerank_local_gmt"
    return res


def enrichment_failed_groups(enrichment_df: pd.DataFrame) -> list[str]:
    if enrichment_df is None or enrichment_df.empty or "Term" not in enrichment_df.columns:
        return []
    failed = enrichment_df[
        enrichment_df["Term"].astype(str).str.startswith("ENRICHMENT_FAILED")
        | enrichment_df["Term"].astype(str).str.startswith("ENRICHMENT_SKIPPED")
    ]
    if failed.empty or "group" not in failed.columns:
        return []
    return sorted(failed["group"].astype(str).unique().tolist())


def top_enrichment_terms(enrichment_df: pd.DataFrame, n_per_group: int = 10) -> pd.DataFrame:
    if enrichment_df.empty:
        return enrichment_df
    df = enrichment_df.copy()
    if "Term" in df.columns:
        df = df[
            ~df["Term"].astype(str).str.startswith("ENRICHMENT_FAILED")
            & ~df["Term"].astype(str).str.startswith("ENRICHMENT_SKIPPED")
        ]
    padj_col = "Adjusted P-value" if "Adjusted P-value" in df.columns else None
    if padj_col is None and "FDR q-val" in df.columns:
        padj_col = "FDR q-val"
    if padj_col is None:
        return df.groupby("group", observed=True).head(n_per_group)
    df = df[df[padj_col].astype(float) <= 0.05]
    return (
        df.sort_values(["group", padj_col], ascending=[True, True])
        .groupby("group", observed=True)
        .head(n_per_group)
        .reset_index(drop=True)
    )
