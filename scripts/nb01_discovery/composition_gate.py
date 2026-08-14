"""Azimuth label composition summary for the teaching cohort (Phase 1).

Reports epithelial airway-lineage membership aligned with the Module 2 four-bucket
map. This is a composition summary, not a certificate that a fifth DE group exists.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

from scripts.common.paths import ensure_output_dirs, resolve, snrna_block_dir
from scripts.nb01_discovery.verify import teaching_inputs

# Airway-lineage epithelial labels from the four-bucket Azimuth map (P0-6).
# Includes Tuft (epithelial), not a separate airway_secretory DE group.
AIRWAY_LINEAGE_EPITHELIAL: frozenset[str] = frozenset(
    {
        "Basal resting",
        "Suprabasal",
        "Multiciliated (non-nasal)",
        "Ionocyte",
        "Tuft",
        "Transitional Club-AT2",
        "Club (non-nasal)",
        "Club (nasal)",
        "Goblet (bronchial)",
        "Goblet (nasal)",
        "SMG duct",
        "SMG serous (bronchial)",
        "SMG mucous",
        "Deuterosomal",
    }
)

# Backward-compatible alias (substring keys previously disagreed on Tuft).
AIRWAY_KEYS = tuple(sorted(AIRWAY_LINEAGE_EPITHELIAL))


def _azimuth_series(path: Path) -> pd.Series:
    """Read azimuth_label from secondary_analysis without loading full X."""
    a = ad.read_h5ad(path, backed="r")
    try:
        if "azimuth_label" in a.obs.columns:
            return a.obs["azimuth_label"].astype(str).copy()
        for col in ("predicted.ann_finest_level", "cell_type", "annotation"):
            if col in a.obs.columns:
                return a.obs[col].astype(str).copy()
        raise KeyError(f"No azimuth-like label column in {path}")
    finally:
        a.file.close()


CROSSTAB_COLUMNS = ["donor_label", "block_id", "azimuth_label", "n_nuclei"]


def _crosstab_from_cache(cfg: dict[str, Any]) -> pd.DataFrame:
    """Per-donor label counts from the committed Azimuth cache.

    secondary_analysis.h5ad is excluded from the minimum learner download, so the
    cache is the only label source on that path. It was built from the same obs
    column of the same files, so the counts match the local-file route exactly.
    """
    from scripts.common.paths import module_root

    rel = ((cfg.get("module2") or {}).get("labels") or {}).get("cache")
    if not rel:
        return pd.DataFrame(columns=CROSSTAB_COLUMNS)
    path = module_root(cfg) / rel
    if not path.exists():
        return pd.DataFrame(columns=CROSSTAB_COLUMNS)

    cache = pd.read_csv(path, sep="\t")
    wanted = {
        str(item["donor_label"]): str(item.get("block_id") or item.get("primary_id") or "")
        for item in teaching_inputs(cfg)
    }
    cache = cache[cache["donor_label"].astype(str).isin(wanted)]
    if cache.empty:
        return pd.DataFrame(columns=CROSSTAB_COLUMNS)
    out = (
        cache.groupby(["donor_label", "azimuth_label"], as_index=False)
        .size()
        .rename(columns={"size": "n_nuclei"})
    )
    out["block_id"] = out["donor_label"].astype(str).map(wanted)
    # Match the local-file route byte for byte: value_counts orders by descending
    # count within a donor, and donors follow teaching_inputs order. A fallback
    # that emits a differently ordered file makes the artifact route-dependent.
    return _order_crosstab(out, cfg)


def _order_crosstab(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    """Deterministic row order, applied to both routes.

    value_counts breaks count ties using pandas internals, which the cache route
    cannot reproduce and which is not stable across versions. Sorting explicitly
    makes the artifact identical whichever source supplied the labels.
    """
    if df.empty:
        return pd.DataFrame(columns=CROSSTAB_COLUMNS)
    donor_order = {str(item["donor_label"]): i for i, item in enumerate(teaching_inputs(cfg))}
    out = df.copy()
    out["_donor_rank"] = out["donor_label"].astype(str).map(donor_order).fillna(len(donor_order))
    out = out.sort_values(
        ["_donor_rank", "n_nuclei", "azimuth_label"], ascending=[True, False, True]
    )
    return out[CROSSTAB_COLUMNS].reset_index(drop=True)


def donor_azimuth_crosstab(cfg: dict[str, Any]) -> pd.DataFrame:
    frames = []
    for item in teaching_inputs(cfg):
        block_dir = snrna_block_dir(cfg, item.get("block_id"), item["donor_label"])
        sec = block_dir / "secondary_analysis.h5ad"
        if not sec.exists():
            continue
        labels = _azimuth_series(sec)
        vc = labels.value_counts().rename("n_nuclei").reset_index()
        vc = vc.rename(columns={vc.columns[0]: "azimuth_label"})
        vc.insert(0, "donor_label", item["donor_label"])
        vc.insert(1, "block_id", item.get("block_id") or item.get("primary_id"))
        frames.append(vc)
    if not frames:
        return _crosstab_from_cache(cfg)
    return _order_crosstab(pd.concat(frames, ignore_index=True), cfg)


def _is_airway_lineage(label: object) -> bool:
    return str(label) in AIRWAY_LINEAGE_EPITHELIAL


def epithelial_airway_lineage_summary(crosstab: pd.DataFrame) -> pd.DataFrame:
    """Per-donor counts of epithelial airway-lineage labels (composition only)."""
    if crosstab.empty:
        # A column-less frame writes a zero-byte TSV that pandas cannot read back
        # (EmptyDataError). Emptiness is a finding; it still needs a header.
        return pd.DataFrame(
            columns=[
                "donor_label",
                "block_id",
                "airway_lineage_nuclei",
                "labels_present",
                "note",
            ]
        )
    mask = crosstab["azimuth_label"].astype(str).map(_is_airway_lineage)
    sub = crosstab.loc[mask].copy()
    if sub.empty:
        donors = crosstab[["donor_label", "block_id"]].drop_duplicates()
        donors["airway_lineage_nuclei"] = 0
        donors["labels_present"] = ""
        donors["note"] = "composition summary only; not a DE group certificate"
        return donors
    g = (
        sub.groupby(["donor_label", "block_id"], as_index=False)
        .agg(
            airway_lineage_nuclei=("n_nuclei", "sum"),
            labels_present=("azimuth_label", lambda s: "; ".join(sorted(set(s)))),
        )
    )
    g["note"] = "composition summary only; not a DE group certificate"
    return g


# Backward-compatible name used by older call sites
def airway_secretory_gate(crosstab: pd.DataFrame) -> pd.DataFrame:
    return epithelial_airway_lineage_summary(crosstab)


def write_composition_gate(cfg: dict[str, Any]) -> dict[str, Path]:
    ensure_output_dirs(cfg)
    tables = resolve(cfg, "outputs_tables")
    ct = donor_azimuth_crosstab(cfg)
    summary = epithelial_airway_lineage_summary(ct)
    airway = ct[ct["azimuth_label"].astype(str).map(_is_airway_lineage)]
    wide = (
        airway.pivot_table(
            index=["donor_label", "block_id"],
            columns="azimuth_label",
            values="n_nuclei",
            aggfunc="sum",
            fill_value=0,
        )
        .reset_index()
        if not airway.empty
        else pd.DataFrame()
    )

    paths = {
        "azimuth_long": tables / "module1_donor_azimuth_composition.tsv",
        "airway_lineage": tables / "module1_epithelial_airway_lineage_composition.tsv",
        "airway_wide": tables / "module1_epithelial_airway_lineage_by_donor.tsv",
    }
    ct.to_csv(paths["azimuth_long"], sep="\t", index=False)
    summary.to_csv(paths["airway_lineage"], sep="\t", index=False)
    if wide.empty:
        note = pd.DataFrame(
            {
                "note": [
                    "no epithelial airway-lineage labels in teaching cohort "
                    "(four-bucket epithelial membership)"
                ]
            }
        )
        note.to_csv(paths["airway_wide"], sep="\t", index=False)
    else:
        wide.to_csv(paths["airway_wide"], sep="\t", index=False)

    # Remove legacy fifth-group gate artifacts if present (P0-6 / Module1 reply).
    for legacy in (
        tables / "module1_airway_secretory_gate.tsv",
        tables / "module1_airway_secretory_by_donor.tsv",
    ):
        if legacy.exists():
            legacy.unlink()
    return paths
