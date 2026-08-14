"""Module 1 teaching-set verification and download-size measurement."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from scripts.common.paths import ensure_output_dirs, module_root, resolve, snrna_block_dir, portable_path


CORE_SNRNA_FILES = ("raw_expr.h5ad",)
OPTIONAL_SNRNA_FILES = ("expr.h5ad", "scvelo_annotated.h5ad", "secondary_analysis.h5ad")
METADATA_GLOB = "*.tsv"


def teaching_inputs(cfg: dict[str, Any]) -> list[dict[str, str]]:
    """Prefer module2.inputs (integration cohort); fall back to module1 single block."""
    inputs = cfg.get("module2", {}).get("inputs") or []
    if inputs:
        return [dict(x) for x in inputs]
    return [
        {
            "donor_label": cfg["module1"]["donor_label"],
            "block_id": cfg["module1"]["block_id"],
        }
    ]


def verify_teaching_files(cfg: dict[str, Any]) -> pd.DataFrame:
    """
    Check each expected teaching file exists and record measured size + sha256.

    Core snRNA file is raw_expr.h5ad only. Azimuth labels come from the cached
    TSV or an optional remote range-read of secondary_analysis (not required).
    """
    from scripts.common.runtime import sha256_optional

    rows: list[dict[str, Any]] = []
    root = module_root(cfg)

    for item in teaching_inputs(cfg):
        donor = item["donor_label"]
        primary = item.get("primary_id") or item.get("block_id")
        processed = item.get("processed_id") or item.get("block_id")
        # Local layout may still use primary folder names; resolution prefers processed.
        block_dir = snrna_block_dir(cfg, processed, donor)
        if not block_dir.exists() and primary:
            block_dir = snrna_block_dir(cfg, primary, donor)
        # The block metadata file is named for the primary block ID. Derive it rather than
        # globbing: a download manifest must name a file the learner does not have yet, and
        # a glob can only ever find files that are already on disk.
        tsvs = sorted(block_dir.glob(METADATA_GLOB)) if block_dir.exists() else []
        metadata_name = tsvs[0].name if tsvs else (f"{primary}.tsv" if primary else None)
        expected = list(CORE_SNRNA_FILES) + ([metadata_name] if metadata_name else [])
        for name in expected:
            path = block_dir / name
            rows.append(
                {
                    "donor_label": donor,
                    "primary_id": primary,
                    "processed_id": processed,
                    "block_id": processed,
                    "file_name": name,
                    "required": True,
                    "exists": path.exists(),
                    "size_mb": round(path.stat().st_size / 1e6, 2) if path.exists() else None,
                    "sha256": sha256_optional(path),
                    "local_path": portable_path(cfg, path),
                    "role": (
                        "Module 1 QC matrix"
                        if name == "raw_expr.h5ad"
                        else "block metadata"
                    ),
                }
            )
        for name in OPTIONAL_SNRNA_FILES:
            path = block_dir / name
            rows.append(
                {
                    "donor_label": donor,
                    "primary_id": primary,
                    "processed_id": processed,
                    "block_id": processed,
                    "file_name": name,
                    "required": False,
                    "exists": path.exists(),
                    "size_mb": round(path.stat().st_size / 1e6, 2) if path.exists() else None,
                    "sha256": sha256_optional(path),
                    "local_path": portable_path(cfg, path),
                    "role": (
                        "optional Azimuth source (use cached labels / remote range-read)"
                        if name == "secondary_analysis.h5ad"
                        else "optional (not required for core path)"
                    ),
                }
            )

    # Module 3 default: MOFA hdf5 only
    m3 = cfg.get("module3") or {}
    mofa = root / m3.get("multiome_dir", "") / (m3.get("files") or {}).get("mofa_hdf5", "multiome_mofa.hdf5")
    rows.append(
        {
            "donor_label": m3.get("donor_label", "Donor_2"),
            "primary_id": m3.get("dataset_id", "HBM828.GPVG.252"),
            "processed_id": m3.get("dataset_id", "HBM828.GPVG.252"),
            "block_id": m3.get("dataset_id", "HBM828.GPVG.252"),
            "file_name": mofa.name,
            "required": True,
            "exists": mofa.exists(),
            "size_mb": round(mofa.stat().st_size / 1e6, 2) if mofa.exists() else None,
            "sha256": sha256_optional(mofa),
            "local_path": portable_path(cfg, mofa),
            "role": "Module 3 MOFA interpretation (default path)",
        }
    )
    h5mu = root / m3.get("multiome_dir", "") / (m3.get("files") or {}).get(
        "secondary_h5mu", "secondary_analysis.h5mu"
    )
    rows.append(
        {
            "donor_label": m3.get("donor_label", "Donor_2"),
            "primary_id": m3.get("dataset_id", "HBM828.GPVG.252"),
            "processed_id": m3.get("dataset_id", "HBM828.GPVG.252"),
            "block_id": m3.get("dataset_id", "HBM828.GPVG.252"),
            "file_name": h5mu.name,
            "required": bool((m3.get("params") or {}).get("load_h5mu", False)),
            "exists": h5mu.exists(),
            "size_mb": round(h5mu.stat().st_size / 1e6, 2) if h5mu.exists() else None,
            "sha256": sha256_optional(h5mu),
            "local_path": portable_path(cfg, h5mu),
            "role": "Module 3 optional full MuData (off by default)",
        }
    )

    # Module 4 source files. Every path the config resolves under data/source/ belongs in
    # the download plan: a learner who fetches only the Module 1 to 3 files cannot run
    # Module 4. Keys are read from module4.data so this list cannot drift from the config.
    m4_data = (cfg.get("module4") or {}).get("data") or {}
    module4_files = [
        ("gtex_tpm", "GTEx", "v10_lung", "Module 4 GTEx lung TPM (levels and HuBMAP bridge)"),
        ("gtex_reads", "GTEx", "v10_lung", "Module 4 GTEx lung counts (aging DE)"),
        ("gtex_sample_attrs", "GTEx", "v10_annotations", "Module 4 GTEx sample metadata"),
        ("gtex_subject_pheno", "GTEx", "v10_annotations", "Module 4 GTEx subject metadata (age brackets)"),
        ("geo_gse150910_counts", "GEO", "GSE150910", "Module 4 IPF vs control counts"),
        ("geo_gse150910_series_matrix", "GEO", "GSE150910", "Module 4 GEO series and sample metadata"),
        ("osd248_counts", "OSDR", "OSD-248", "Module 4 RSEM counts (rRNArm) for the flight contrast"),
        ("osd248_genelab_de", "OSDR", "OSD-248", "Module 4 GeneLab DE table (independent cross-check)"),
        ("osd248_genelab_contrasts", "OSDR", "OSD-248", "Module 4 GeneLab contrast definitions"),
        ("osd248_isa_sample", "OSDR", "OSD-248", "Module 4 ISA sample metadata (arm assignment)"),
        ("osd248_isa_assay", "OSDR", "OSD-248", "Module 4 ISA assay metadata (library prep, sequencer)"),
    ]
    for key, label, accession, role in module4_files:
        rel = m4_data.get(key)
        if not rel:
            continue
        path = root / rel
        rows.append(
            {
                "donor_label": label,
                "primary_id": accession,
                "processed_id": accession,
                "block_id": accession,
                "file_name": path.name,
                "required": True,
                "exists": path.exists(),
                "size_mb": round(path.stat().st_size / 1e6, 2) if path.exists() else None,
                "sha256": sha256_optional(path),
                "local_path": str(rel),
                "role": role,
            }
        )

    return pd.DataFrame(rows)


def learner_download_manifest(cfg: dict[str, Any]) -> pd.DataFrame:
    """Required files only, measured sizes -- the learner download list."""
    df = verify_teaching_files(cfg)
    req = df[df["required"]].copy()
    req["included_in_learner_total"] = req["exists"] & req["required"]
    return req


def write_verification_outputs(cfg: dict[str, Any]) -> dict[str, Path]:
    ensure_output_dirs(cfg)
    tables = resolve(cfg, "outputs_tables")
    verify = verify_teaching_files(cfg)
    manifest = learner_download_manifest(cfg)

    missing = verify[(verify["required"]) & (~verify["exists"])]
    # size_mb is decimal megabytes (bytes/1e6). Report both decimal GB and binary GiB.
    total_mb = float(manifest.loc[manifest["exists"], "size_mb"].sum()) if len(manifest) else 0.0
    total_bytes = total_mb * 1e6
    summary = pd.DataFrame(
        [
            {
                "measured_utc": datetime.now(timezone.utc).isoformat(),
                "n_required": int(verify["required"].sum()),
                "n_required_present": int(((verify["required"]) & (verify["exists"])).sum()),
                "n_required_missing": int(len(missing)),
                "learner_total_mb": round(total_mb, 2),
                "learner_total_gb": round(total_mb / 1000.0, 3),
                "learner_total_gib": round(total_bytes / (1024.0**3), 3),
                "missing_files": "; ".join(
                    f"{r.donor_label}/{r.block_id}/{r.file_name}" for r in missing.itertuples()
                )
                or "",
                "note": (
                    "Total is measured on this machine and covers every file the config resolves "
                    "under data/source/ for Modules 1 to 4 (HuBMAP raw_expr + block metadata, "
                    "Module 3 MOFA object, GTEx, GEO and OSDR). "
                    "size_mb is decimal MB (bytes/1e6); learner_total_gb = MB/1000; "
                    "learner_total_gib = bytes/1024^3. "
                    "secondary_analysis.h5ad is optional; Azimuth labels use cached TSV "
                    "or remote range-read. Optional expr/scvelo/h5mu excluded. "
                    "SmartAPI live discovery deferred."
                ),
            }
        ]
    )

    paths = {
        "verification": tables / "module1_file_verification.tsv",
        "download_manifest": tables / "module1_learner_download_manifest.tsv",
        "download_summary": tables / "module1_learner_download_summary.tsv",
    }
    verify.to_csv(paths["verification"], sep="\t", index=False)
    manifest.to_csv(paths["download_manifest"], sep="\t", index=False)
    summary.to_csv(paths["download_summary"], sep="\t", index=False)
    return paths
