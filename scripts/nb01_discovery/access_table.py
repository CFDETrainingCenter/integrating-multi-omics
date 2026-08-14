"""Build Module 1 dataset access / provenance table (multi-donor + GTEx)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from scripts.common.paths import donor_root, ensure_output_dirs, module_root, resolve
from scripts.common.provenance import parse_block_metadata
from scripts.nb01_discovery.inventory import inventory_all_donors


PROVENANCE_COLUMNS = [
    "resource_name",
    "cfde_role",
    "native_portal",
    "donor_id",
    "primary_dataset_id",
    "processed_dataset_note",
    "file_name",
    "file_format",
    "file_size_mb",
    "data_type",
    "organ_or_tissue",
    "assay",
    "subject_or_donor_scope",
    "access_status",
    "core_or_optional",
    "module_role",
    "provenance_notes",
    "redistribution_status",
    "local_relative_path",
]


def _role_for_snrna_file(name: str) -> tuple[str, str]:
    if name.endswith(".tsv"):
        return "donor/sample metadata for provenance", "core"
    if name == "raw_expr.h5ad":
        return "Module 1 QC starting matrix", "core"
    if name == "secondary_analysis.h5ad":
        return "optional Azimuth labels for interpretation", "recommended"
    if name == "expr.h5ad":
        return "alternate processed matrix; not primary QC input", "optional"
    if name == "scvelo_annotated.h5ad":
        return "optional RNA velocity (not Module 3 core)", "optional"
    return "supporting", "optional"


def _snrna_rows(cfg: dict[str, Any], inventory: pd.DataFrame) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    snrna_inv = inventory[inventory["modality"] == "snRNAseq"].copy()
    for donor_label, sub in snrna_inv.groupby("donor_label"):
        donor_key = f"donor_{str(donor_label).split('_')[-1]}" if "_" in str(donor_label) else "donor_1"
        donor = cfg.get(donor_key, cfg["donor_1"])
        root = donor_root(cfg, donor_label)
        for block_id in sorted(sub["block_or_dataset_id"].dropna().unique()):
            block_dir = root / "snRNAseq" / block_id
            tsvs = list(block_dir.glob("*.tsv"))
            meta = parse_block_metadata(tsvs[0]) if tsvs else {}
            block_files = sub[sub["block_or_dataset_id"] == block_id]
            for _, f in block_files.iterrows():
                module_role, core = _role_for_snrna_file(f["file_name"])
                rows.append(
                    {
                        "resource_name": "HuBMAP",
                        "cfde_role": "DCC-native analysis-ready files (discovered via CFDE ecosystem)",
                        "native_portal": "https://portal.hubmapconsortium.org",
                        "donor_id": meta.get("donor_id", donor["donor_id"]),
                        "primary_dataset_id": block_id,
                        "processed_dataset_note": "Files sourced from RNAseq [Salmon] child of this primary",
                        "file_name": f["file_name"],
                        "file_format": f["file_format"],
                        "file_size_mb": f["size_mb"],
                        "data_type": "snRNA-seq",
                        "organ_or_tissue": donor.get("organ", "Lung"),
                        "assay": meta.get("rnaseq_assay_method") or "snRNAseq (10x Genomics v3)",
                        "subject_or_donor_scope": (
                            f"{donor['label']} | age={meta.get('age_value', donor.get('age'))} | "
                            f"sex={meta.get('sex', donor.get('sex'))} | race={meta.get('race', donor.get('race'))}"
                        ),
                        "access_status": "downloaded locally (public processed products)",
                        "core_or_optional": core,
                        "module_role": module_role,
                        "provenance_notes": (
                            f"Donor {donor['donor_id']}; primary block {block_id}; "
                            "folder named by primary ID, contains Salmon processed products"
                        ),
                        "redistribution_status": "not redistributed; learner downloads from source portal",
                        "local_relative_path": f["relative_path"],
                    }
                )
    return rows


def _atac_rows(cfg: dict[str, Any], inventory: pd.DataFrame) -> list[dict[str, Any]]:
    """Optional Donor_1 SnapATAC files (historical / unpaired; not the Module 5 core path)."""
    rows: list[dict[str, Any]] = []
    atac_inv = inventory[inventory["modality"] == "snATACseq"].copy()
    donor = cfg["donor_1"]
    for _, f in atac_inv.iterrows():
        name = f["file_name"]
        if name in {"cell_by_bin.h5ad", "cell_by_gene.h5ad"}:
            module_role, core = (
                "optional unpaired ATAC (not Module 3 core; Module 3 uses the configured SNARE Muon)",
                "optional",
            )
        else:
            module_role, core = "supporting ATAC file", "optional"

        rows.append(
            {
                "resource_name": "HuBMAP",
                "cfde_role": "DCC-native analysis-ready files (discovered via CFDE ecosystem)",
                "native_portal": "https://portal.hubmapconsortium.org",
                "donor_id": donor["donor_id"],
                "primary_dataset_id": f["block_or_dataset_id"],
                "processed_dataset_note": "SnapATAC processed aggregate from SNARE-seq2 ATAC libraries",
                "file_name": name,
                "file_format": f["file_format"],
                "file_size_mb": f["size_mb"],
                "data_type": "snATAC-seq (SNARE-seq2)",
                "organ_or_tissue": donor.get("organ", "Lung"),
                "assay": "snATACseq (SNARE-seq2) [SnapATAC]",
                "subject_or_donor_scope": (
                    f"{donor['label']} | age={donor.get('age')} | sex={donor.get('sex')} | "
                    f"race={donor.get('race')}"
                ),
                "access_status": "downloaded locally" if f.get("size_mb") else "may be absent on disk",
                "core_or_optional": core,
                "module_role": module_role,
                "provenance_notes": (
                    "Same donor as Donor_1 10x snRNA-seq; unpaired multimodal "
                    "(donor-matched, not same-cell multiome). Module 3 core path is the configured SNARE Muon."
                ),
                "redistribution_status": "not redistributed; learner downloads from source portal",
                "local_relative_path": f["relative_path"],
            }
        )
    return rows


def _snare_multiome_rows(cfg: dict[str, Any], inventory: pd.DataFrame) -> list[dict[str, Any]]:
    """Configured SNARE-seq2 Muon multiome products used by Module 3."""
    rows: list[dict[str, Any]] = []
    m3 = cfg.get("module3", {})
    donor_label = m3.get("donor_label", "")
    dataset_id = m3.get("dataset_id", "")
    donor = {}
    for key in ("donor_1", "donor_2", "donor_3", "donor_4"):
        block = cfg.get(key) or {}
        if block.get("label") == donor_label:
            donor = block
            break
    snare = inventory[inventory["modality"].isin(["SNAREseq", "SNARE-seq", "SNARE"])].copy()
    if dataset_id and not snare.empty and "block_or_dataset_id" in snare.columns:
        snare = snare[snare["block_or_dataset_id"].astype(str) == str(dataset_id)]
    if snare.empty:
        from scripts.nb05_multimodal.load import multiome_paths

        paths = multiome_paths(cfg)
        for key, role, core, dtype in (
            ("mofa_hdf5", "Module 3 core MOFA factors", "core", "MOFA HDF5"),
            ("secondary_h5mu", "Module 3 optional MuData (not a learner download)", "optional", "multiome MuData"),
            ("pdf_combined", "Module 3 portal Leiden PDF (combined)", "recommended", "PDF"),
            ("mudata_raw_h5mu", "Module 3 optional raw MuData", "optional", "multiome MuData"),
        ):
            path = paths.get(key)
            if path is None:
                continue
            if not path.exists() and key in {"mudata_raw_h5mu", "secondary_h5mu"}:
                continue
            try:
                rel = path.relative_to(module_root(cfg) / "data" / "source" / "HUBMAP")
            except ValueError:
                rel = path
            rows.append(
                {
                    "resource_name": "HuBMAP",
                    "cfde_role": "DCC-native analysis-ready files (discovered via CFDE ecosystem)",
                    "native_portal": m3.get("portal_url", "https://portal.hubmapconsortium.org"),
                    "donor_id": m3.get("donor_id", donor.get("donor_id", "")),
                    "primary_dataset_id": m3.get("dataset_id", ""),
                    "processed_dataset_note": "SNARE-seq2 [Salmon + ArchR + Muon] true multiome product",
                    "file_name": path.name,
                    "file_format": path.suffix.lstrip(".") or "unknown",
                    "file_size_mb": round(path.stat().st_size / 1e6, 1) if path.exists() else None,
                    "data_type": dtype,
                    "organ_or_tissue": donor.get("organ", "Lung"),
                    "assay": m3.get("assay", "SNARE-seq2 [Salmon + ArchR + Muon]"),
                    "subject_or_donor_scope": (
                        f"{m3.get('donor_label', donor_label)} | age={donor.get('age')} | "
                        f"sex={donor.get('sex')} | race={donor.get('race')}"
                    ),
                    "access_status": "downloaded locally" if path.exists() else "missing",
                    "core_or_optional": core,
                    "module_role": role,
                    "provenance_notes": (
                        f"Dataset {m3.get('dataset_id')}; true shared-nucleus multiome. "
                        "Modules 1-2 and 4 use 10x snRNA-seq; Module 3 uses this SNARE Muon product."
                    ),
                    "redistribution_status": "not redistributed; learner downloads from source portal",
                    "local_relative_path": str(rel),
                }
            )
        return rows

    for _, f in snare.iterrows():
        name = f["file_name"]
        if name == "secondary_analysis.h5mu":
            module_role, core = "Module 3 optional MuData (not a learner download)", "optional"
        elif name == "multiome_mofa.hdf5":
            module_role, core = "Module 3 core MOFA factors", "core"
        elif name.startswith("leiden_cluster") and name.endswith(".pdf"):
            module_role, core = "Module 3 portal Leiden PDF", "recommended"
        elif name == "mudata_raw.h5mu":
            module_role, core = "Module 3 optional raw MuData", "optional"
        else:
            module_role, core = "Module 3 supporting multiome file", "optional"
        rows.append(
            {
                "resource_name": "HuBMAP",
                "cfde_role": "DCC-native analysis-ready files (discovered via CFDE ecosystem)",
                "native_portal": m3.get("portal_url", "https://portal.hubmapconsortium.org"),
                "donor_id": m3.get("donor_id", f.get("donor_id", donor.get("donor_id", ""))),
                "primary_dataset_id": f.get("block_or_dataset_id", m3.get("dataset_id", "")),
                "processed_dataset_note": "SNARE-seq2 [Salmon + ArchR + Muon] true multiome product",
                "file_name": name,
                "file_format": f["file_format"],
                "file_size_mb": f["size_mb"],
                "data_type": "SNARE-seq2 multiome",
                "organ_or_tissue": donor.get("organ", "Lung"),
                "assay": m3.get("assay", "SNARE-seq2 [Salmon + ArchR + Muon]"),
                "subject_or_donor_scope": (
                    f"{m3.get('donor_label', f.get('donor_label', donor_label))} | "
                    f"age={donor.get('age')} | sex={donor.get('sex')} | race={donor.get('race')}"
                ),
                "access_status": "downloaded locally",
                "core_or_optional": core,
                "module_role": module_role,
                "provenance_notes": (
                    f"Dataset {m3.get('dataset_id', f.get('block_or_dataset_id'))}; "
                    "true shared-nucleus multiome. Modules 1-2 and 4 use 10x snRNA-seq; "
                    "Module 3 uses this SNARE Muon product."
                ),
                "redistribution_status": "not redistributed; learner downloads from source portal",
                "local_relative_path": f["relative_path"],
            }
        )
    return rows

def _gtex_rows(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    gtex_dir = module_root(cfg) / "data" / "source" / "gtex"
    rows = []
    # (relative path under gtex/, dataset id, role, core/optional, dtype)
    mapping = [
        (
            "gene_tpm_v10_lung.gct.gz",
            "GTEx_v10_adult_lung",
            "Module 3 primary tissue-level lung TPM reference (best HuBMAP gene-model match)",
            "core",
            "TPM",
        ),
        (
            "gene_reads_v10_lung.gct.gz",
            "GTEx_v10_adult_lung",
            "Module 3 optional read-count support matrix",
            "recommended",
            "read counts",
        ),
        (
            "v11/gene_tpm_adult_gtex_v11_lung.gct.gz",
            "GTEx_v11_adult_lung",
            "Optional version-comparison TPM (annotation-drift lesson)",
            "optional",
            "TPM",
        ),
        (
            "v11/gene_reads_adult_gtex_v11_lung.gct.gz",
            "GTEx_v11_adult_lung",
            "Optional version-comparison read counts",
            "optional",
            "read counts",
        ),
    ]
    for rel, dataset_id, role, core, dtype in mapping:
        path = gtex_dir / rel
        if not path.exists():
            continue
        name = Path(rel).name
        rows.append(
            {
                "resource_name": "GTEx",
                "cfde_role": "DCC-native tissue-level reference (not single-cell equivalent)",
                "native_portal": "https://gtexportal.org",
                "donor_id": "multi-subject adult GTEx lung cohort",
                "primary_dataset_id": dataset_id,
                "processed_dataset_note": "Adult lung bulk RNA-seq matrix (GCT)",
                "file_name": name,
                "file_format": "gct.gz",
                "file_size_mb": round(path.stat().st_size / 1e6, 1),
                "data_type": dtype,
                "organ_or_tissue": "Lung",
                "assay": "bulk RNA-seq",
                "subject_or_donor_scope": "adult GTEx lung samples",
                "access_status": "downloaded locally",
                "core_or_optional": core,
                "module_role": role,
                "provenance_notes": (
                    "Use as tissue-level reference context for HuBMAP snRNA-seq; "
                    "do not treat as single-nucleus equivalent. "
                    "v10 preferred for Ensembl gene-model overlap with HuBMAP."
                ),
                "redistribution_status": "not redistributed; learner downloads from source portal",
                "local_relative_path": f"gtex/{rel}",
            }
        )
    return rows


def build_access_provenance_table(cfg: dict[str, Any]) -> pd.DataFrame:
    inventory = inventory_all_donors(cfg)
    rows = (
        _snrna_rows(cfg, inventory)
        + _atac_rows(cfg, inventory)
        + _snare_multiome_rows(cfg, inventory)
        + _gtex_rows(cfg)
    )
    df = pd.DataFrame(rows)
    for col in PROVENANCE_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    return df[PROVENANCE_COLUMNS]


def write_access_provenance_table(cfg: dict[str, Any]) -> Path:
    ensure_output_dirs(cfg)
    df = build_access_provenance_table(cfg)
    out = resolve(cfg, "outputs_tables") / "module1_dataset_access_provenance.tsv"
    df.to_csv(out, sep="\t", index=False)
    return out
