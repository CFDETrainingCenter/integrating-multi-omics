"""Loaders for GSE150910, OSD-248, GTEx, HuBMAP pseudobulk (Module 4)."""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from scripts.common.paths import module_root
from scripts.nb03_integration.gtex import load_gtex_lung_tpm, strip_ensembl_version


def _p(cfg: dict[str, Any], rel: str) -> Path:
    p = Path(rel)
    return p if p.is_absolute() else module_root(cfg) / p


# ---------------------------------------------------------------------------
# GSE150910
# ---------------------------------------------------------------------------


def load_gse150910_meta(cfg: dict[str, Any]) -> pd.DataFrame:
    """Parse GEO series matrix characteristics into one row per sample title."""
    path = _p(cfg, cfg["module4"]["data"]["geo_gse150910_series_matrix"])
    with gzip.open(path, "rt") as f:
        lines = f.readlines()
    acc = titles = None
    char_lines: list[list[str]] = []
    for line in lines:
        if line.startswith("!Sample_geo_accession"):
            acc = re.findall(r'"([^"]+)"', line)
        elif line.startswith("!Sample_title"):
            titles = re.findall(r'"([^"]+)"', line)
        elif line.startswith("!Sample_characteristics_ch1"):
            char_lines.append(re.findall(r'"([^"]+)"', line))
    if not acc or not titles:
        raise ValueError(f"Could not parse sample accession/title from {path}")
    n = len(acc)
    sample_chars = [{} for _ in range(n)]
    for cl in char_lines:
        for i, v in enumerate(cl):
            if ": " in v:
                k, val = v.split(": ", 1)
                sample_chars[i][k.strip().lower()] = val.strip()
    meta = pd.DataFrame(sample_chars)
    meta["gsm"] = acc
    meta["title"] = titles
    meta["diagnosis"] = meta["diagnosis"].astype(str).str.lower()
    return meta


def load_gse150910_counts(cfg: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Return (counts genesxsamples, sample_meta with condition).

    Teaching contrast: IPF vs control (CHP excluded).
    """
    counts_path = _p(cfg, cfg["module4"]["data"]["geo_gse150910_counts"])
    counts = pd.read_csv(counts_path, index_col=0)
    meta = load_gse150910_meta(cfg).set_index("title")
    include = {d.lower() for d in cfg["module4"]["geo"].get("include_diagnoses", ["ipf", "control"])}
    keep_titles = [t for t in counts.columns if t in meta.index and meta.loc[t, "diagnosis"] in include]
    counts = counts[keep_titles].copy()
    sample_meta = meta.loc[keep_titles].copy()
    sample_meta["condition"] = sample_meta["diagnosis"].map(
        lambda d: "IPF" if d == "ipf" else ("control" if d == "control" else d)
    )
    sample_meta["resource"] = "GSE150910"
    sample_meta["contrast_id"] = "H_DISEASE"
    return counts, sample_meta


# ---------------------------------------------------------------------------
# OSD-248
# ---------------------------------------------------------------------------


def _read_isa_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str)


def load_osd248_isa(cfg: dict[str, Any]) -> pd.DataFrame:
    """Join ISA sample + assay (QA Score) on Sample Name."""
    data = cfg["module4"]["data"]
    sample = _read_isa_table(_p(cfg, data["osd248_isa_sample"]))
    assay = _read_isa_table(_p(cfg, data["osd248_isa_assay"]))
    qa_col = cfg["module4"]["osd248"].get("qa_score_column", "Parameter Value[QA Score]")
    keep_assay = ["Sample Name"]
    if qa_col in assay.columns:
        keep_assay.append(qa_col)
    assay_s = assay[keep_assay].drop_duplicates(subset=["Sample Name"])
    merged = sample.merge(assay_s, on="Sample Name", how="left")
    if qa_col in merged.columns:
        merged["qa_score"] = pd.to_numeric(merged[qa_col], errors="coerce")
    else:
        merged["qa_score"] = np.nan
    return merged


def osd248_teaching_samples(cfg: dict[str, Any]) -> pd.DataFrame:
    """
    ISS-T Flight ~60 On ISS Carcass vs Ground Control ~60 On Earth Carcass.

    Filter on Factor Value columns -- not free-text group names.
    """
    isa = load_osd248_isa(cfg)
    o = cfg["module4"]["osd248"]
    dur = str(o.get("duration", "~60"))
    # Duration may appear as "~60" or "~60 day"
    dur_ok = isa["Factor Value[Duration]"].astype(str).str.replace(" day", "", regex=False).str.strip() == dur

    flight = (
        (isa["Factor Value[Spaceflight]"].astype(str) == o["spaceflight_flight"])
        & dur_ok
        & (isa["Factor Value[Euthanasia Location]"].astype(str) == o["euthanasia_location_flight"])
        & (isa["Factor Value[Dissection Condition]"].astype(str) == o["dissection_condition"])
    )
    ground = (
        (isa["Factor Value[Spaceflight]"].astype(str) == o["spaceflight_ground"])
        & dur_ok
        & (isa["Factor Value[Euthanasia Location]"].astype(str) == o["euthanasia_location_ground"])
        & (isa["Factor Value[Dissection Condition]"].astype(str) == o["dissection_condition"])
    )
    out = isa.loc[flight | ground].copy()
    out["condition"] = np.where(flight.loc[out.index], "flight", "ground")
    out["resource"] = "OSD-248"
    out["contrast_id"] = "M_FLIGHT"
    # Sample Name in ISA matches count column names
    out["sample_id"] = out["Sample Name"].astype(str)
    return out.reset_index(drop=True)


def load_osd248_counts(cfg: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Unnormalized counts for teaching contrast samples."""
    path = _p(cfg, cfg["module4"]["data"]["osd248_counts"])
    counts = pd.read_csv(path, index_col=0)
    samples = osd248_teaching_samples(cfg)
    cols = [c for c in samples["sample_id"] if c in counts.columns]
    missing = sorted(set(samples["sample_id"]) - set(cols))
    if missing:
        raise KeyError(f"OSD-248 count columns missing for samples: {missing[:5]}...")
    counts = counts[cols].copy()
    meta = samples.set_index("sample_id").loc[cols]
    return counts, meta


# ---------------------------------------------------------------------------
# HuBMAP pseudobulk + GTEx bridge
# ---------------------------------------------------------------------------


def load_hubmap_pseudobulk(cfg: dict[str, Any]) -> tuple[pd.Series, pd.Series, pd.Series | None]:
    """
    Composition-weighted (linear) and max-across-labels vectors from Module 2.

    Returns (weighted, max_across, fractions_optional).
    """
    data = cfg["module4"]["data"]
    w = pd.read_csv(_p(cfg, data["hubmap_pseudobulk"]), sep="\t", index_col=0).iloc[:, 0]
    w.name = "composition_weighted_linear"
    m = pd.read_csv(_p(cfg, data["hubmap_max_across"]), sep="\t", index_col=0).iloc[:, 0]
    m.name = "max_across_labels_linear"
    fracs = None
    fr_path = data.get("hubmap_fractions")
    if fr_path and _p(cfg, fr_path).exists():
        fracs = pd.read_csv(_p(cfg, fr_path), sep="\t", index_col=0).iloc[:, 0]
        fracs.name = "fraction"
    return w, m, fracs


def load_gtex_means_for_bridge(cfg: dict[str, Any]) -> pd.DataFrame:
    """GTEx lung mean TPM indexed by versionless Ensembl, with symbols."""
    gtex = load_gtex_lung_tpm(cfg)
    sample_cols = [c for c in gtex.columns if c not in {"gene_id", "gene_symbol"}]
    means = gtex[sample_cols].astype(float).mean(axis=1)
    out = pd.DataFrame(
        {
            "gene_symbol": gtex["gene_symbol"].astype(str),
            "gtex_lung_mean_tpm": means,
        }
    )
    out.index.name = "gene_id_no_version"
    return out


# ---------------------------------------------------------------------------
# Sample counts inventory
# ---------------------------------------------------------------------------


def sample_counts_table(cfg: dict[str, Any]) -> pd.DataFrame:
    """Actual n per arm from loaded data (honesty requirement)."""
    rows = []
    # GSE150910
    _, geo_meta = load_gse150910_counts(cfg)
    for cond, n in geo_meta["condition"].value_counts().items():
        rows.append(
            {
                "resource": "GSE150910",
                "contrast_id": "H_DISEASE",
                "arm": cond,
                "n": int(n),
                "notes": "CHP excluded from teaching contrast",
            }
        )
    # OSD-248 flight vs GC
    osd = osd248_teaching_samples(cfg)
    for cond, n in osd["condition"].value_counts().items():
        rows.append(
            {
                "resource": "OSD-248",
                "contrast_id": "M_FLIGHT",
                "arm": cond,
                "n": int(n),
                "notes": "ISS-T ~60 Carcass; duration-matched GC (not hindlimb unloading)",
            }
        )
    low = osd[
        osd["qa_score"]
        < float(cfg["module4"]["osd248"].get("rin_sensitivity_threshold", 6.0))
    ]
    rows.append(
        {
            "resource": "OSD-248",
            "contrast_id": "M_FLIGHT",
            "arm": "flight_qa_below_threshold",
            "n": int(len(low)),
            "notes": ";".join(low["sample_id"].tolist()) if len(low) else "none",
        }
    )
    # GTEx age contrast design (predeclared brackets)
    try:
        from scripts.nb06_crossres.contrasts import gtex_age_bracket_inventory

        inv = gtex_age_bracket_inventory(cfg)
        gtex_cfg = cfg["module4"].get("gtex") or {}
        young_b = gtex_cfg.get("young_age_brackets") or ["20-29", "30-39"]
        old_b = gtex_cfg.get("old_age_brackets") or ["60-69", "70-79"]
        for label, brackets in (("younger_brackets", young_b), ("older_brackets", old_b)):
            sub = inv[inv["age_bracket"].isin(brackets)]
            rows.append(
                {
                    "resource": "GTEx",
                    "contrast_id": "H_AGING_GTEX",
                    "arm": label,
                    "n": int(sub["n_subjects"].sum()),
                    "notes": (
                        f"brackets={list(brackets)}; samples={int(sub['n_samples'].sum())}; "
                        f"DE uses one sample/subject capped at "
                        f"{int(gtex_cfg.get('max_subjects_per_arm', 40))}/arm"
                    ),
                }
            )
    except Exception as exc:  # noqa: BLE001
        rows.append(
            {
                "resource": "GTEx",
                "contrast_id": "H_AGING_GTEX",
                "arm": "inventory_failed",
                "n": 0,
                "notes": str(exc),
            }
        )
    # HuBMAP age illustration (2 vs 2)
    try:
        from scripts.nb06_crossres.hubmap_donor_counts import donor_age_meta

        hm = donor_age_meta(cfg)
        for cond, sub in hm.groupby("condition"):
            rows.append(
                {
                    "resource": "HuBMAP",
                    "contrast_id": "H_AGING_HUBMAP",
                    "arm": cond,
                    "n": int(len(sub)),
                    "notes": (
                        "illustration only; donors="
                        + ",".join(sub["donor_label"].astype(str))
                        + f"; ages="
                        + ",".join(sub["age_value"].astype(str))
                    ),
                }
            )
    except Exception as exc:  # noqa: BLE001
        rows.append(
            {
                "resource": "HuBMAP",
                "contrast_id": "H_AGING_HUBMAP",
                "arm": "meta_failed",
                "n": 0,
                "notes": str(exc),
            }
        )
    # Level bridge inventory
    w, m, _ = load_hubmap_pseudobulk(cfg)
    rows.append(
        {
            "resource": "HuBMAP",
            "contrast_id": "bridge",
            "arm": "composition_weighted_genes",
            "n": int(len(w)),
            "notes": "linear-scale pseudobulk from Module 2 (levels only)",
        }
    )
    gtex = load_gtex_means_for_bridge(cfg)
    rows.append(
        {
            "resource": "GTEx",
            "contrast_id": "bridge",
            "arm": "lung_mean_tpm_genes",
            "n": int(len(gtex)),
            "notes": "mean TPM for levels/bridge only; age DE uses gene_reads",
        }
    )
    return pd.DataFrame(rows)


def native_levels_panel(cfg: dict[str, Any]) -> pd.DataFrame:
    """
    Gene panel native values side-by-side -- deliberately incomparable units.

    This is the demonstration table, not the analysis.
    """
    from scripts.nb06_crossres.orthologs import load_ortholog_table, map_human_to_mouse

    panel = [g.upper() for g in cfg["module4"]["params"]["gene_panel"]]
    weighted, max_vec, _ = load_hubmap_pseudobulk(cfg)
    gtex = load_gtex_means_for_bridge(cfg)
    geo_counts, geo_meta = load_gse150910_counts(cfg)
    osd_counts, osd_meta = load_osd248_counts(cfg)
    ortho = load_ortholog_table(cfg)
    mapped = map_human_to_mouse(panel, ortho).set_index("human_symbol")

    rows = []
    for sym in panel:
        row = {
            "gene_symbol_human": sym,
            "gene_symbol_mouse": (
                str(mapped.loc[sym, "mouse_symbol"]) if sym in mapped.index else ""
            ),
            "hubmap_composition_weighted_linear": np.nan,
            "hubmap_max_across_labels_linear": np.nan,
            "gtex_lung_mean_tpm": np.nan,
            "gse150910_mean_counts_control": np.nan,
            "gse150910_mean_counts_ipf": np.nan,
            "osd248_mean_counts_ground": np.nan,
            "osd248_mean_counts_flight": np.nan,
            "units_note": (
                "TPM vs linear CP10K-pseudobulk vs raw counts; human vs mouse -- not comparable"
            ),
        }
        hit = gtex[gtex["gene_symbol"].astype(str).str.upper() == sym]
        if len(hit):
            row["gtex_lung_mean_tpm"] = float(hit["gtex_lung_mean_tpm"].iloc[0])
            ens = str(hit.index[0])
            if ens in weighted.index:
                row["hubmap_composition_weighted_linear"] = float(weighted.loc[ens])
            if ens in max_vec.index:
                row["hubmap_max_across_labels_linear"] = float(max_vec.loc[ens])
        if sym in geo_counts.index:
            ctrl_cols = geo_meta.index[geo_meta["condition"] == "control"]
            ipf_cols = geo_meta.index[geo_meta["condition"] == "IPF"]
            row["gse150910_mean_counts_control"] = float(geo_counts.loc[sym, ctrl_cols].mean())
            row["gse150910_mean_counts_ipf"] = float(geo_counts.loc[sym, ipf_cols].mean())
        # OSD by mouse ensembl via ortholog table
        if sym in mapped.index and "mouse_ensembl" in ortho.columns:
            mens = ortho.loc[ortho["human_symbol"] == sym, "mouse_ensembl"].astype(str)
            mens = mens.str.replace(r"\.\d+$", "", regex=True)
            osd_idx = osd_counts.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
            for me in mens:
                hits = np.where(osd_idx == me)[0]
                if len(hits):
                    gid = osd_counts.index[int(hits[0])]
                    gnd = osd_meta.index[osd_meta["condition"] == "ground"]
                    flt = osd_meta.index[osd_meta["condition"] == "flight"]
                    row["osd248_mean_counts_ground"] = float(osd_counts.loc[gid, gnd].mean())
                    row["osd248_mean_counts_flight"] = float(osd_counts.loc[gid, flt].mean())
                    break
        rows.append(row)
    return pd.DataFrame(rows)
