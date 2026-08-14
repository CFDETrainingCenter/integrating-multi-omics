"""Scan Donor_* trees into inventory tables."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from scripts.common.paths import donor_root, list_donor_keys, portable_path, resolve
from scripts.common.provenance import file_role, parse_block_metadata


def _iter_files(root: Path) -> Iterable[Path]:
    if not root.exists():
        return []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != ".DS_Store":
            yield path


def _inventory_donor(cfg: dict[str, Any], donor_key: str) -> pd.DataFrame:
    donor = cfg[donor_key]
    root = donor_root(cfg, donor["label"])
    rows: list[dict[str, Any]] = []
    for path in _iter_files(root):
        rel = path.relative_to(root)
        parts = rel.parts
        modality = parts[0] if parts else ""
        block_id = parts[1] if len(parts) > 1 else ""
        rows.append(
            {
                "donor_label": donor["label"],
                "donor_id": donor["donor_id"],
                "modality": modality,
                "block_or_dataset_id": block_id,
                "file_name": path.name,
                "relative_path": str(Path(donor["label"]) / rel),
                "absolute_path": portable_path(cfg, path),
                "file_format": path.suffix.lstrip(".") or "unknown",
                "size_bytes": path.stat().st_size,
                "size_mb": round(path.stat().st_size / 1e6, 1),
                "file_role": file_role(path.name),
            }
        )
    return pd.DataFrame(rows)






def inventory_all_donors(cfg: dict[str, Any]) -> pd.DataFrame:
    frames = [_inventory_donor(cfg, key) for key in list_donor_keys(cfg)]
    frames = [f for f in frames if not f.empty]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def summarize_inventory(inventory: pd.DataFrame) -> pd.DataFrame:
    if inventory.empty:
        return pd.DataFrame()
    return (
        inventory.groupby(["donor_label", "modality", "block_or_dataset_id"], dropna=False)
        .agg(
            n_files=("file_name", "count"),
            total_mb=("size_mb", "sum"),
            files=("file_name", lambda s: ", ".join(sorted(s))),
        )
        .reset_index()
        .sort_values(["donor_label", "modality", "block_or_dataset_id"])
    )


def snrna_block_metadata_table(cfg: dict[str, Any], donor_key: str = "donor_1") -> pd.DataFrame:
    donor = cfg[donor_key]
    root = donor_root(cfg, donor["label"])
    snrna = root / "snRNAseq"
    rows = []
    if not snrna.exists():
        return pd.DataFrame()
    for block_dir in sorted(p for p in snrna.iterdir() if p.is_dir()):
        tsvs = list(block_dir.glob("*.tsv"))
        if not tsvs:
            continue
        meta = parse_block_metadata(tsvs[0])
        present = {
            "has_raw_expr": (block_dir / "raw_expr.h5ad").exists(),
            "has_expr": (block_dir / "expr.h5ad").exists(),
            "has_secondary_analysis": (block_dir / "secondary_analysis.h5ad").exists(),
            "has_scvelo": (block_dir / "scvelo_annotated.h5ad").exists(),
        }
        rows.append(
            {
                "donor_label": donor["label"],
                "block_id": block_dir.name,
                **{
                    k: meta[k]
                    for k in (
                        "donor_id",
                        "age_value",
                        "sex",
                        "race",
                        "ethnicity",
                        "rnaseq_assay_method",
                        "assay_description",
                    )
                },
                **present,
                "module1_candidate": present["has_raw_expr"],
                "annotation_candidate": present["has_secondary_analysis"],
                "complete_core_set": present["has_raw_expr"] and present["has_secondary_analysis"],
                "is_teaching_block": block_dir.name
                == str(donor.get("teaching_block") or cfg.get("module1", {}).get("block_id") or ""),
            }
        )
    # Mark teaching blocks from module2.inputs when present
    teaching = {
        (x["donor_label"], x["block_id"])
        for x in (cfg.get("module2", {}).get("inputs") or [])
    }
    if teaching and rows:
        for r in rows:
            r["is_teaching_block"] = (r["donor_label"], r["block_id"]) in teaching
    return pd.DataFrame(rows)


def snrna_block_metadata_all(cfg: dict[str, Any]) -> pd.DataFrame:
    frames = [snrna_block_metadata_table(cfg, key) for key in list_donor_keys(cfg)]
    frames = [f for f in frames if not f.empty]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def write_inventory_outputs(cfg: dict[str, Any]) -> dict[str, Path]:
    from scripts.common.paths import ensure_output_dirs

    ensure_output_dirs(cfg)
    inv = inventory_all_donors(cfg)
    summary = summarize_inventory(inv)
    blocks = snrna_block_metadata_all(cfg)

    tables = resolve(cfg, "outputs_tables")
    paths = {
        "inventory": tables / "module1_file_inventory.tsv",
        "summary": tables / "module1_inventory_summary.tsv",
        "blocks": tables / "module1_snrna_blocks.tsv",
    }
    inv.to_csv(paths["inventory"], sep="\t", index=False)
    summary.to_csv(paths["summary"], sep="\t", index=False)
    blocks.to_csv(paths["blocks"], sep="\t", index=False)
    # Teaching-block subset for the learner verification view
    if not blocks.empty and "is_teaching_block" in blocks.columns:
        teach = blocks[blocks["is_teaching_block"]].copy()
        teach_path = tables / "module1_snrna_teaching_blocks.tsv"
        teach.to_csv(teach_path, sep="\t", index=False)
        paths["teaching_blocks"] = teach_path
    return paths
