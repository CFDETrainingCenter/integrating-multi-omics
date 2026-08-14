"""Per-resource differential expression (pydeseq2) -- logFC vectors only cross boundaries.

Module 4 rebuild (2026-08-10): four contrasts share one pydeseq2 path on count-like
inputs (GTEx gene_reads, GEO counts, OSD-248 RSEM rounded, HuBMAP summed/rounded).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


def _is_ercc_spike_id(gene_id: object) -> bool:
    """True for ERCC spike-in IDs (ERCC-...), not human ERCC* repair genes."""
    return str(gene_id).upper().startswith("ERCC-")


def _pydeseq2_contrast(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    *,
    condition_col: str = "condition",
    reference: str,
    contrast_level: str,
    contrast_id: str,
    min_count_sum: int = 10,
    size_factor_gene_ids: list[str] | None = None,
) -> pd.DataFrame:
    """
    Run pydeseq2 Wald test for ``contrast_level`` vs ``reference``.

    ``counts`` is genes x samples; ``meta`` indexed by sample id.

    Polarity is never left to alphabetical factor order:
    - condition is a categorical with ``reference`` as the first / Treatment baseline
    - DeseqStats contrast is always ``[factor, tested_level, ref_level]``

    ``size_factor_gene_ids`` (optional): gene IDs to use for size-factor estimation
    (pydeseq2 ``control_genes``). Spike-ins that are confounded with arm should be
    omitted here while remaining in the count matrix so their Wald stats are still
    returned for calibration plots.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    if reference == contrast_level:
        raise ValueError(f"{contrast_id}: reference and contrast_level are both {reference!r}")

    counts = counts.loc[:, meta.index]
    counts_i = np.floor(counts.fillna(0).clip(lower=0)).astype(int)
    keep = counts_i.sum(axis=1) >= int(min_count_sum)
    counts_i = counts_i.loc[keep]
    counts_t = counts_i.T
    counts_t.index = counts_t.index.astype(str)

    metadata = meta.copy()
    metadata.index = metadata.index.astype(str)
    levels_present = [str(x) for x in metadata[condition_col].astype(str).unique()]
    if reference not in levels_present or contrast_level not in levels_present:
        raise KeyError(
            f"{contrast_id}: need levels {reference!r} and {contrast_level!r} in "
            f"{condition_col}; found {sorted(levels_present)}"
        )
    # Explicit order: reference first. Remaining levels keep stable sort for determinism.
    other = sorted(x for x in levels_present if x not in {reference, contrast_level})
    cat_levels = [reference, contrast_level, *other]
    metadata[condition_col] = pd.Categorical(
        metadata[condition_col].astype(str),
        categories=cat_levels,
        ordered=False,
    )

    # Treatment() pins the design-matrix baseline; contrast= pins the reported LFC sign.
    design = f"~C({condition_col}, Treatment('{reference}'))"
    inference = DefaultInference(n_cpus=1)
    control_genes = None
    if size_factor_gene_ids is not None:
        keep_sf = [g for g in counts_t.columns if str(g) in set(map(str, size_factor_gene_ids))]
        if not keep_sf:
            raise ValueError(f"{contrast_id}: size_factor_gene_ids matched zero genes after filters")
        control_genes = keep_sf
    dds = DeseqDataSet(
        counts=counts_t,
        metadata=metadata,
        design=design,
        refit_cooks=True,
        inference=inference,
        control_genes=control_genes,
    )
    dds.deseq2()
    contrast = [condition_col, contrast_level, reference]
    stat = DeseqStats(
        dds,
        contrast=contrast,
        inference=inference,
    )
    stat.summary()
    res = stat.results_df.copy()
    res = res.rename(
        columns={
            "log2FoldChange": "logFC",
            "padj": "padj",
            "pvalue": "pvalue",
            "baseMean": "baseMean",
        }
    )
    res.index.name = "gene_symbol"
    res = res.reset_index()
    res["contrast_id"] = contrast_id
    res["reference"] = reference
    res["contrast_level"] = contrast_level
    res["design_formula"] = design
    res["contrast_vector"] = f"{condition_col}: {contrast_level} vs {reference}"
    res.attrs["polarity"] = {
        "contrast_id": contrast_id,
        "condition_col": condition_col,
        "reference": reference,
        "contrast_level": contrast_level,
        "categorical_levels": cat_levels,
        "design_formula": design,
        "deseq_stats_contrast": contrast,
        "alphabetical_first_level": sorted(levels_present)[0],
        "alphabetical_would_differ": sorted(levels_present)[0] != reference,
        "size_factors_exclude_ercc": bool(
            control_genes is not None
            and any(_is_ercc_spike_id(g) for g in counts_t.columns)
            and not any(_is_ercc_spike_id(g) for g in control_genes)
        ),
        "n_size_factor_genes": int(len(control_genes)) if control_genes is not None else int(counts_t.shape[1]),
        "n_genes_in_dds": int(counts_t.shape[1]),
    }
    return res


def contrast_h_disease(cfg: dict[str, Any]) -> pd.DataFrame:
    """GSE150910 IPF vs control (condition only -- published analysis used covariates)."""
    from scripts.nb06_crossres.load import load_gse150910_counts

    counts, meta = load_gse150910_counts(cfg)
    min_sum = int((cfg["module4"].get("thresholds") or {}).get("de_min_count_sum", 10))
    return _pydeseq2_contrast(
        counts,
        meta,
        condition_col="condition",
        reference="control",
        contrast_level="IPF",
        contrast_id="H_DISEASE",
        min_count_sum=min_sum,
    )


def contrast_m_flight(
    cfg: dict[str, Any],
    *,
    drop_low_qa: bool = False,
    exclude_ercc_from_size_factors: bool | None = None,
) -> pd.DataFrame:
    """OSD-248 ISS-T flight vs ground control (renamed from M_UNLOAD).

    Flight=ExFold Mix1 and GC=Mix2, so ERCC abundance is confounded with arm.
    By default, ``ERCC-`` rows stay in the DE table for calibration but are omitted
    from size-factor estimation (``module4.osd248.exclude_ercc_from_size_factors``).
    """
    from scripts.nb06_crossres.load import load_osd248_counts
    from scripts.nb06_crossres.orthologs import load_ortholog_table

    counts, meta = load_osd248_counts(cfg)
    thr = float(cfg["module4"]["osd248"].get("rin_sensitivity_threshold", 6.0))
    if drop_low_qa and "qa_score" in meta.columns:
        meta = meta[~(meta["qa_score"] < thr)].copy()
        counts = counts.loc[:, meta.index]
    min_sum = int((cfg["module4"].get("thresholds") or {}).get("de_min_count_sum", 10))
    if exclude_ercc_from_size_factors is None:
        exclude_ercc_from_size_factors = bool(
            (cfg["module4"].get("osd248") or {}).get("exclude_ercc_from_size_factors", True)
        )
    size_factor_ids = None
    if exclude_ercc_from_size_factors:
        size_factor_ids = [g for g in counts.index.astype(str) if not _is_ercc_spike_id(g)]
    res = _pydeseq2_contrast(
        counts,
        meta,
        condition_col="condition",
        reference="ground",
        contrast_level="flight",
        contrast_id="M_FLIGHT" + ("_qa_filtered" if drop_low_qa else ""),
        min_count_sum=min_sum,
        size_factor_gene_ids=size_factor_ids,
    )
    polarity = dict(getattr(res, "attrs", {}).get("polarity") or {})
    res = res.rename(columns={"gene_symbol": "mouse_ensembl"})
    ortho = load_ortholog_table(cfg)
    res["mouse_ensembl_base"] = res["mouse_ensembl"].astype(str).str.replace(
        r"\.\d+$", "", regex=True
    )
    ortho = ortho.copy()
    if "mouse_ensembl" in ortho.columns:
        ortho["mouse_ensembl_base"] = ortho["mouse_ensembl"].astype(str).str.replace(
            r"\.\d+$", "", regex=True
        )
        res = res.merge(
            ortho[
                ["mouse_ensembl_base", "mouse_symbol", "human_symbol"]
            ].drop_duplicates("mouse_ensembl_base"),
            on="mouse_ensembl_base",
            how="left",
        )
    else:
        res["mouse_symbol"] = pd.NA
        res["human_symbol"] = pd.NA
    res["gene_symbol"] = (
        res["human_symbol"].fillna(res["mouse_symbol"]).fillna(res["mouse_ensembl"])
    )
    if polarity:
        res.attrs["polarity"] = polarity
    return res


def _gtex_lung_age_table(cfg: dict[str, Any]) -> pd.DataFrame:
    """Lung samples joined to AGE brackets (one row per SAMPID)."""
    from scripts.common.paths import module_root

    m4 = cfg["module4"]
    data = m4.get("data") or {}
    pheno_path = module_root(cfg) / (
        data.get("gtex_subject_pheno")
        or "data/source/GTEX/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
    )
    attrs_path = module_root(cfg) / (
        data.get("gtex_sample_attrs")
        or "data/source/GTEX/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
    )
    pheno = pd.read_csv(pheno_path, sep="\t", dtype=str)
    attrs = pd.read_csv(attrs_path, sep="\t", dtype=str, low_memory=False)
    lung = attrs[attrs["SMTSD"].astype(str).str.contains("Lung", case=False, na=False)].copy()
    # Sample IDs look like GTEX-1117F-0526-SM-5EGHJ; subject is first two fields
    lung["SUBJID"] = lung["SAMPID"].astype(str).str.split("-").str[:2].str.join("-")
    lung = lung.merge(pheno[["SUBJID", "AGE"]], on="SUBJID", how="inner")
    return lung


def gtex_age_bracket_inventory(cfg: dict[str, Any]) -> pd.DataFrame:
    """Count lung samples and subjects per AGE bracket (declare design before DE)."""
    lung = _gtex_lung_age_table(cfg)
    rows = []
    for age, sub in lung.groupby("AGE", sort=True):
        rows.append(
            {
                "age_bracket": age,
                "n_samples": int(len(sub)),
                "n_subjects": int(sub["SUBJID"].nunique()),
            }
        )
    return pd.DataFrame(rows)


def contrast_h_aging_gtex(cfg: dict[str, Any]) -> pd.DataFrame:
    """
    GTEx lung older vs younger age brackets on gene_reads (pydeseq2).

    Brackets are predeclared in config (default 20-39 vs 60-79). Age is entangled
    with cause of death and ischemic time -- state that confound when reporting.
    """
    from pathlib import Path

    from scripts.common.paths import module_root
    from scripts.nb03_integration.gtex import read_gct_gz

    m4 = cfg["module4"]
    gtex_cfg = m4.get("gtex") or {}
    data = m4.get("data") or {}
    reads_path = module_root(cfg) / (
        data.get("gtex_reads") or "data/source/GTEX/gene_reads_v10_lung.gct.gz"
    )
    if not Path(reads_path).exists():
        raise FileNotFoundError(f"Missing GTEx gene_reads: {reads_path}")

    young_brackets = [
        str(x) for x in gtex_cfg.get("young_age_brackets") or ["20-29", "30-39"]
    ]
    old_brackets = [
        str(x) for x in gtex_cfg.get("old_age_brackets") or ["60-69", "70-79"]
    ]
    lung = _gtex_lung_age_table(cfg)
    young = lung[lung["AGE"].isin(young_brackets)].copy()
    old = lung[lung["AGE"].isin(old_brackets)].copy()

    gct = read_gct_gz(reads_path)
    sample_cols = [c for c in gct.columns if c not in {"gene_id", "gene_symbol"}]

    def _pick_ids(arm: pd.DataFrame) -> list[str]:
        # Prefer one sample per subject to avoid pseudoreplication
        if bool(gtex_cfg.get("one_sample_per_subject", True)):
            arm = arm.drop_duplicates(subset=["SUBJID"], keep="first")
        ids = [s for s in arm["SAMPID"].astype(str) if s in sample_cols]
        max_n = int(gtex_cfg.get("max_subjects_per_arm", 40))
        seed = int(gtex_cfg.get("random_seed", 0))
        if len(ids) > max_n:
            rng = np.random.default_rng(seed)
            ids = list(rng.choice(ids, size=max_n, replace=False))
        return ids

    y_ids = _pick_ids(young)
    o_ids = _pick_ids(old)
    if len(y_ids) < 5 or len(o_ids) < 5:
        raise RuntimeError(
            f"H_AGING_GTEX arms too small after sampling: young={len(y_ids)} old={len(o_ids)}"
        )

    mat = gct[y_ids + o_ids].copy()
    symbols = gct["gene_symbol"].astype(str).str.upper()
    mat.index = symbols.values
    mat = mat[~mat.index.duplicated(keep="first")]

    subjid_map = lung.drop_duplicates("SAMPID").set_index("SAMPID")["SUBJID"].astype(str)
    age_map = lung.drop_duplicates("SAMPID").set_index("SAMPID")["AGE"].astype(str)
    all_ids = y_ids + o_ids
    meta = pd.DataFrame(
        {
            "condition": ["younger"] * len(y_ids) + ["older"] * len(o_ids),
            "SUBJID": [subjid_map.get(s, "-".join(str(s).split("-")[:2])) for s in all_ids],
            "AGE": [age_map.get(s, "") for s in all_ids],
        },
        index=all_ids,
    )
    min_sum = int((m4.get("thresholds") or {}).get("de_min_count_sum", 10))
    res = _pydeseq2_contrast(
        mat,
        meta,
        condition_col="condition",
        reference="younger",
        contrast_level="older",
        contrast_id="H_AGING_GTEX",
        min_count_sum=min_sum,
    )
    res.attrs["design"] = {
        "young_age_brackets": young_brackets,
        "old_age_brackets": old_brackets,
        "n_younger": len(y_ids),
        "n_older": len(o_ids),
        "one_sample_per_subject": bool(gtex_cfg.get("one_sample_per_subject", True)),
        "max_subjects_per_arm": int(gtex_cfg.get("max_subjects_per_arm", 40)),
        "random_seed": int(gtex_cfg.get("random_seed", 0)),
        "note": (
            "Age entangled with cause of death and ischemic time; "
            "brackets predeclared before DE"
        ),
    }
    return res


def contrast_h_aging(cfg: dict[str, Any]) -> pd.DataFrame | None:
    """Backward-compatible name -> ``contrast_h_aging_gtex``."""
    return contrast_h_aging_gtex(cfg)


def contrast_h_aging_hubmap(cfg: dict[str, Any]) -> pd.DataFrame:
    """
    HuBMAP donor summed counts: older (Donor_3+7) vs younger (Donor_1+3).

    Illustration only (2 vs 2). Do not loosen padj thresholds. Label everywhere
    as design illustration, never as a powered result.
    """
    from scripts.nb06_crossres.hubmap_donor_counts import (
        donor_age_meta,
        summed_counts_by_donor,
    )

    pb = summed_counts_by_donor(cfg)
    meta_rows = donor_age_meta(cfg).set_index("donor_label")
    donor_cols = [c for c in meta_rows.index if c in pb.columns]
    if len(donor_cols) < 4:
        raise RuntimeError(
            f"H_AGING_HUBMAP needs four donors in count matrix; found {donor_cols}"
        )

    counts = pb[donor_cols].copy()
    # Index by gene symbol when available (human NES space)
    if "gene_symbol" in pb.columns:
        symbols = pb["gene_symbol"].astype(str).str.upper()
        counts.index = symbols.values
        counts = counts[~counts.index.duplicated(keep="first")]
    else:
        counts.index = counts.index.astype(str).str.replace(r"\.\d+$", "", regex=True)

    meta = meta_rows.loc[donor_cols, ["condition", "age_value", "contrast_id"]].copy()
    min_sum = int((cfg["module4"].get("thresholds") or {}).get("de_min_count_sum", 10))
    res = _pydeseq2_contrast(
        counts,
        meta,
        condition_col="condition",
        reference="younger",
        contrast_level="older",
        contrast_id="H_AGING_HUBMAP",
        min_count_sum=min_sum,
    )
    res.attrs["design"] = {
        "older_donors": list(
            meta_rows.index[meta_rows["condition"] == "older"]
        ),
        "younger_donors": list(
            meta_rows.index[meta_rows["condition"] == "younger"]
        ),
        "ages": meta_rows["age_value"].to_dict(),
        "illustration_only": True,
        "note": (
            "2 vs 2 design illustration. HuBMAP older ages sit inside GTEx middle "
            "brackets; power is an order of magnitude below H_AGING_GTEX. "
            "Counts layer summed then rounded (not integer UMIs)."
        ),
    }
    return res


def logfc_series(de: pd.DataFrame, gene_col: str = "gene_symbol") -> pd.Series:
    s = de.set_index(gene_col)["logFC"].astype(float)
    s = s[~s.index.duplicated(keep="first")]
    return s
