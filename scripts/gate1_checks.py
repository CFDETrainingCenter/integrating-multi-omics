"""Gate 1 assertions for Stage 1 (no full module re-run).

Reports:
1. Four-bucket Azimuth coverage (empty unmapped list)
2. New cell_class sizes vs expected
3. Prerank guard rejects one-sided ranked lists
4. is_ercc_symbol distinguishes spike-ins from ERCC repair genes
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main() -> int:
    from scripts.nb03_integration.integrate import (
        _AZIMUTH_TO_COARSE,
        assert_azimuth_mapping_coverage,
        cell_class_membership_frame,
    )
    from scripts.nb04_de.enrich import run_prerank_gsea
    from scripts.nb06_crossres.compare import is_ercc_symbol
    from scripts.common.paths import load_config

    cfg = load_config()
    frac_path = ROOT / "outputs" / "tables" / "module2_celltype_fractions.tsv"
    # Prefer absolute counts from donor crosstab when present; else scale fractions.
    xtab_path = ROOT / "outputs" / "tables" / "module2_donor_by_label_crosstab.tsv"
    expected = {
        "epithelial": 6928,
        "endothelial": 6205,
        "stromal": 3190,
        "immune": 2556,
        "unlabeled": 5,
    }

    if xtab_path.exists():
        xt = pd.read_csv(xtab_path, sep="\t")
        # Wide: donor_label rows x azimuth_label columns
        if "donor_label" in xt.columns:
            counts = xt.drop(columns=["donor_label"]).sum(axis=0)
        elif "azimuth_label" in xt.columns and "n" in xt.columns:
            counts = xt.groupby("azimuth_label", observed=True)["n"].sum()
        elif "azimuth_label" in xt.columns:
            value_cols = [c for c in xt.columns if c != "azimuth_label"]
            counts = xt.set_index("azimuth_label")[value_cols].sum(axis=1)
        else:
            label_col = xt.columns[0]
            counts = xt.set_index(label_col).sum(axis=1)
        labels = counts.index.astype(str)
        membership = cell_class_membership_frame(
            pd.Series(counts.index.astype(str).repeat(counts.astype(int).to_numpy()))
        )
    else:
        fr = pd.read_csv(frac_path, sep="\t")
        label_col = "azimuth_label" if "azimuth_label" in fr.columns else fr.columns[0]
        frac_col = "fraction" if "fraction" in fr.columns else fr.columns[1]
        n_total = 18884
        counts = (fr.set_index(label_col)[frac_col] * n_total).round().astype(int)
        labels = counts.index.astype(str)
        membership = cell_class_membership_frame(
            pd.Series(counts.index.astype(str).repeat(counts.to_numpy()))
        )

    unmapped = assert_azimuth_mapping_coverage(labels)
    print("=== GATE 1.1 coverage ===")
    print("labels present:", len(list(labels)))
    for lab in sorted(labels):
        if str(lab) in {"unlabeled"}:
            compartment = "unlabeled"
        else:
            compartment = _AZIMUTH_TO_COARSE.get(str(lab), f"UNMAPPED:{lab}")
        n = int(counts.get(lab, counts.get(str(lab), 0)))
        print(f"  {lab!r:40s} -> {compartment:16s} n={n}")
    print("unmapped:", unmapped)

    sizes = (
        membership.groupby("cell_class", observed=True)["n_nuclei"]
        .sum()
        .reindex(["epithelial", "endothelial", "stromal", "immune", "unlabeled"])
        .fillna(0)
        .astype(int)
    )
    print("\n=== GATE 1.2 group sizes ===")
    print(sizes.to_string())
    print("expected:", expected)
    ok_sizes = True
    for k, v in expected.items():
        got = int(sizes.get(k, 0))
        if got != v:
            # Allow +/-1 from fraction rounding when crosstab absent
            if abs(got - v) > 1:
                ok_sizes = False
                print(f"SIZE MISMATCH {k}: got {got} expected {v}")

    print("\n=== GATE 1.3 prerank guard ===")
    one_sided = pd.Series(
        {"GENEA": 3.0, "GENEB": 2.0, "GENEC": 1.0, "GENED": 0.5},
        name="wilcoxon_score_demo",
    )
    guard_fired = False
    try:
        run_prerank_gsea(one_sided, cfg)
    except ValueError as exc:
        guard_fired = "no negative scores" in str(exc).lower() or "one-sided" in str(exc).lower()
        print("raised:", exc)
    print("guard_fired:", guard_fired)

    print("\n=== GATE 1.4 ERCC ===")
    ercc_human = is_ercc_symbol("ERCC1")
    ercc_spike = is_ercc_symbol("ERCC-00130")
    print(f'is_ercc_symbol("ERCC1") = {ercc_human} (want False)')
    print(f'is_ercc_symbol("ERCC-00130") = {ercc_spike} (want True)')

    ok = (
        unmapped == []
        and ok_sizes
        and guard_fired
        and (ercc_human is False)
        and (ercc_spike is True)
    )
    print("\n=== GATE 1 VERDICT ===")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
