"""Extract ISS-T flight vs ground GeneLab DE columns into a small committed TSV.

Reads the ~111 MB rRNArm differential-expression CSV plus the contrasts CSV,
keeps ENSEMBL/SYMBOL identity columns and the ISS-T ~60d Carcass flight-vs-GC
stats, and writes ``outputs/tables/module4_osd248_genelab_isst_de_subset.tsv``.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd

from scripts.common.paths import load_config, module_root

# Teaching contrast: Space Flight ~60 On ISS Carcass vs Ground Control ~60 On Earth Carcass
ISST_CONTRAST_SUBSTR = (
    "(Space Flight & ~60 day & On ISS & Carcass)v"
    "(Ground Control & ~60 day & On Earth & Carcass)"
)

ID_COLS = ("ENSEMBL", "SYMBOL", "GENENAME", "ENTREZID")
STAT_PREFIXES = ("Log2fc_", "Stat_", "P.value_", "Adj.p.value_")
GROUP_MEAN_PREFIXES = ("Group.Mean_", "Group.Stdev_")
GROUP_FLIGHT = "Space Flight & ~60 day & On ISS & Carcass"
GROUP_GC = "Ground Control & ~60 day & On Earth & Carcass"


def _resolve_data_path(cfg: dict[str, Any], key: str, default_rel: str) -> Path:
    data = (cfg.get("module4") or {}).get("data") or {}
    rel = data.get(key) or default_rel
    p = Path(rel)
    if p.is_absolute():
        return p
    return module_root(cfg) / p


def isst_contrast_columns(columns: list[str]) -> list[str]:
    """Select GeneLab DE columns for the ISS-T flight vs GC teaching contrast."""
    keep: list[str] = []
    for c in columns:
        if c in ID_COLS:
            keep.append(c)
            continue
        if ISST_CONTRAST_SUBSTR in c and any(c.startswith(p) for p in STAT_PREFIXES):
            keep.append(c)
            continue
        if c in {
            f"Group.Mean_({GROUP_FLIGHT})",
            f"Group.Stdev_({GROUP_FLIGHT})",
            f"Group.Mean_({GROUP_GC})",
            f"Group.Stdev_({GROUP_GC})",
        }:
            keep.append(c)
    return keep


def extract_genelab_isst_de_subset(
    cfg: dict[str, Any] | None = None,
    *,
    de_csv: Path | str | None = None,
    contrasts_csv: Path | str | None = None,
    out_tsv: Path | str | None = None,
) -> Path:
    """Extract ISS-T flight vs ground columns and write a small TSV."""
    cfg = cfg or load_config()
    de_path = Path(de_csv) if de_csv else _resolve_data_path(
        cfg,
        "osd248_genelab_de",
        "data/source/OSDR/OSD-248/GLDS-248_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv",
    )
    contrasts_path = Path(contrasts_csv) if contrasts_csv else _resolve_data_path(
        cfg,
        "osd248_genelab_contrasts",
        "data/source/OSDR/OSD-248/GLDS-248_rna_seq_contrasts_GLbulkRNAseq.csv",
    )
    out_path = Path(out_tsv) if out_tsv else _resolve_data_path(
        cfg,
        "osd248_genelab_isst_de_subset",
        "outputs/tables/module4_osd248_genelab_isst_de_subset.tsv",
    )

    if contrasts_path.exists():
        # Soft check: confirm the teaching contrast appears among published contrasts
        contrast_header = pd.read_csv(contrasts_path, nrows=0).columns.tolist()
        if not any(ISST_CONTRAST_SUBSTR in str(c) for c in contrast_header):
            # contrasts file uses the bare contrast string as a column name (no Log2fc_ prefix)
            if not any(
                "Space Flight & ~60 day & On ISS & Carcass" in str(c)
                and "Ground Control & ~60 day & On Earth & Carcass" in str(c)
                for c in contrast_header
            ):
                raise KeyError(
                    f"ISS-T teaching contrast not found in {contrasts_path}"
                )

    header = pd.read_csv(de_path, nrows=0).columns.tolist()
    cols = isst_contrast_columns(header)
    if not any(c.startswith("Log2fc_") for c in cols):
        raise KeyError(
            f"No Log2fc column for ISS-T teaching contrast in {de_path}"
        )

    df = pd.read_csv(de_path, usecols=cols, low_memory=False)
    # Friendly short names for the teaching contrast stats
    rename = {}
    for c in cols:
        if ISST_CONTRAST_SUBSTR in c:
            if c.startswith("Log2fc_"):
                rename[c] = "log2fc_isst_flight_vs_gc"
            elif c.startswith("Stat_"):
                rename[c] = "stat_isst_flight_vs_gc"
            elif c.startswith("P.value_"):
                rename[c] = "pvalue_isst_flight_vs_gc"
            elif c.startswith("Adj.p.value_"):
                rename[c] = "padj_isst_flight_vs_gc"
        elif c == f"Group.Mean_({GROUP_FLIGHT})":
            rename[c] = "group_mean_isst_flight"
        elif c == f"Group.Stdev_({GROUP_FLIGHT})":
            rename[c] = "group_stdev_isst_flight"
        elif c == f"Group.Mean_({GROUP_GC})":
            rename[c] = "group_mean_isst_gc"
        elif c == f"Group.Stdev_({GROUP_GC})":
            rename[c] = "group_stdev_isst_gc"
    df = df.rename(columns=rename)
    df.insert(0, "contrast_id", "M_FLIGHT_ISST_flight_vs_gc")
    df.attrs["source_de"] = str(de_path)
    df.attrs["contrast"] = ISST_CONTRAST_SUBSTR

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    return out_path


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--de-csv", default=None)
    parser.add_argument("--contrasts-csv", default=None)
    parser.add_argument("--out-tsv", default=None)
    args = parser.parse_args(argv)
    cfg = load_config()
    out = extract_genelab_isst_de_subset(
        cfg,
        de_csv=args.de_csv,
        contrasts_csv=args.contrasts_csv,
        out_tsv=args.out_tsv,
    )
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
