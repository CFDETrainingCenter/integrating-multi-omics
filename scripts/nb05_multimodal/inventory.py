"""Inventory configured SNARE Muon multiome assets for Module 3."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from scripts.common.paths import module_root
from scripts.nb05_multimodal.load import multiome_paths


def _rel(cfg: dict[str, Any], path: Path) -> str:
    root = module_root(cfg)
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path)


def inventory_module3_assets(cfg: dict[str, Any]) -> pd.DataFrame:
    m3 = cfg["module3"]
    paths = multiome_paths(cfg)
    portal = m3.get("portal_url", "https://portal.hubmapconsortium.org")
    load_h5mu = bool((m3.get("params") or {}).get("load_h5mu", False))
    donor = m3["donor_label"]
    specs = [
        (
            "mofa_hdf5",
            "MOFA",
            "Fitted MOFA factors (rna + atac_cbg views)",
            "core",
            "default interpretation path (no full MuData download)",
        ),
        (
            "label_cache",
            "cache",
            "Committed Azimuth / WNN label cache",
            "core",
            "in-repo derived table; learners do not download the h5mu",
        ),
        (
            "secondary_h5mu",
            "SNARE_Muon",
            "Paired multiome MuData (rna + atac_cbg + atac_cbb)",
            "optional" if not load_h5mu else "core",
            "optional extension; not a learner download on the default path",
        ),
        (
            "pdf_combined",
            "figure",
            "Portal combined Leiden PDF",
            "recommended",
            "Reference figure from HuBMAP pipeline",
        ),
        (
            "pdf_rna",
            "figure",
            "Portal RNA Leiden PDF",
            "recommended",
            "Reference figure from HuBMAP pipeline",
        ),
        (
            "pdf_atac",
            "figure",
            "Portal ATAC Leiden PDF",
            "recommended",
            "Reference figure from HuBMAP pipeline",
        ),
        (
            "mudata_raw_h5mu",
            "SNARE_Muon",
            "Less-processed MuData (optional)",
            "optional",
            "Not required for Module 3 core path",
        ),
    ]
    rows = []
    matched = m3.get("matched_snrna") or {}
    for key, modality, role, core, notes in specs:
        path = paths.get(key)
        if path is None:
            continue
        path = Path(path)
        if key == "mudata_raw_h5mu" and not path.exists():
            continue
        rows.append(
            {
                "donor_label": donor,
                "donor_id": m3["donor_id"],
                "dataset_id": m3["dataset_id"],
                "portal_url": portal,
                "assay": m3.get("assay", "SNARE-seq2 [Salmon + ArchR + Muon]"),
                "modality": modality,
                "file_name": path.name,
                "local_path": _rel(cfg, path),
                "exists": path.exists(),
                "size_mb": round(path.stat().st_size / 1e6, 1) if path.exists() else None,
                "module_role": role,
                "core_or_optional": core,
                "notes": notes,
                "matched_snrna_primary_id": matched.get("primary_id", ""),
                "provenance_notes": (
                    f"HuBMAP dataset {m3['dataset_id']} from donor {m3['donor_id']}; "
                    f"Modules 1-2/4 use 10x snRNA-seq -- Module 3 is the SNARE multiome "
                    f"for {donor}"
                ),
            }
        )
    return pd.DataFrame(rows)


def inventory_module5_assets(cfg: dict[str, Any]) -> pd.DataFrame:
    return inventory_module3_assets(cfg)
