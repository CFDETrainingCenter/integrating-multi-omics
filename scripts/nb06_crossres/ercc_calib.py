"""ERCC ExFold calibration for OSD-248 ISS-T flight vs ground (Module 4 / P6-6).

Flight = Mix 1, ground control = Mix 2 (ISA Parameter Value[Spike-in Mix Number]).
Expected log2(Mix1/Mix2) by Thermo ExFold subgroup: A=+2, B=0, C=-0.58, D=-1.
ERCC is perfectly confounded with arm in this contrast -- calibration, not biology.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from scripts.common.paths import module_root
from scripts.common.plotting import save_figure
from scripts.nb06_crossres.compare import is_ercc_symbol

# Teaching contrast: ISS-T flight vs GC (see FIX_LIST P6-6).
FLIGHT_MIX = 1
GROUND_MIX = 2
EXPECTED_BY_SUBGROUP = {"A": 2.0, "B": 0.0, "C": -0.58, "D": -1.0}


def ercc_expected_table(cfg: dict[str, Any]) -> pd.DataFrame:
    """Load Thermo Fisher ExFold expected Mix1/Mix2 fold changes (92 ERCC IDs)."""
    data = (cfg.get("module4") or {}).get("data") or {}
    rel = data.get("ercc_exfold_expected") or "data/reference/ercc_exfold_expected_fc.txt"
    path = Path(rel)
    if not path.is_absolute():
        path = module_root(cfg) / path
    if not path.exists():
        raise FileNotFoundError(f"Missing ERCC ExFold reference table: {path}")
    df = pd.read_csv(path, sep="\t")
    # Normalize Thermo column names
    rename = {
        "ERCC ID": "ercc_id",
        "subgroup": "subgroup",
        "log2(Mix 1/Mix 2)": "expected_log2fc_mix1_over_mix2",
        "expected fold-change ratio": "expected_fc_ratio",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    df["ercc_id"] = df["ercc_id"].astype(str).str.strip()
    df["subgroup"] = df["subgroup"].astype(str).str.strip().str.upper()
    if "expected_log2fc_mix1_over_mix2" not in df.columns:
        df["expected_log2fc_mix1_over_mix2"] = df["subgroup"].map(EXPECTED_BY_SUBGROUP)
    df["expected_log2fc_mix1_over_mix2"] = df["expected_log2fc_mix1_over_mix2"].astype(float)
    df["flight_mix"] = FLIGHT_MIX
    df["ground_mix"] = GROUND_MIX
    df["note"] = (
        "ISS-T teaching contrast: flight=Mix1, GC=Mix2; observed log2FC is flight/GC = Mix1/Mix2"
    )
    return df[["ercc_id", "subgroup", "expected_log2fc_mix1_over_mix2", "flight_mix", "ground_mix", "note"]]


def ercc_read_fraction_by_sample(cfg: dict[str, Any]) -> pd.DataFrame:
    """ERCC spike-in read fraction per OSD-248 teaching sample, with arm label.

    Systematically higher ERCC fraction in flight vs ground is the size-factor
    confound mechanism when Mix1/Mix2 differ and spike-ins enter normalization.
    """
    from scripts.nb06_crossres.load import load_osd248_counts

    counts, meta = load_osd248_counts(cfg)
    counts_i = np.floor(counts.fillna(0).clip(lower=0)).astype(float)
    is_ercc = counts_i.index.astype(str).map(lambda g: str(g).upper().startswith("ERCC-"))
    ercc_sum = counts_i.loc[is_ercc].sum(axis=0)
    total = counts_i.sum(axis=0)
    frac = (ercc_sum / total.replace(0, np.nan)).astype(float)
    out = pd.DataFrame(
        {
            "sample_id": counts_i.columns.astype(str),
            "condition": meta.loc[counts_i.columns, "condition"].astype(str).values,
            "ercc_counts": ercc_sum.values,
            "total_counts": total.values,
            "ercc_fraction": frac.values,
            "n_ercc_genes_nonzero": (counts_i.loc[is_ercc] > 0).sum(axis=0).astype(int).values,
        }
    )
    return out.sort_values(["condition", "sample_id"]).reset_index(drop=True)


def ercc_read_fraction_arm_summary(by_sample: pd.DataFrame) -> dict[str, Any]:
    """Compact arm-level summary for methods / run_params."""
    if by_sample is None or by_sample.empty:
        return {}
    summary: dict[str, Any] = {"n_samples": int(len(by_sample))}
    for arm, sub in by_sample.groupby("condition"):
        summary[f"{arm}_n"] = int(len(sub))
        summary[f"{arm}_mean_ercc_fraction"] = float(sub["ercc_fraction"].mean())
        summary[f"{arm}_median_ercc_fraction"] = float(sub["ercc_fraction"].median())
        summary[f"{arm}_mean_ercc_counts"] = float(sub["ercc_counts"].mean())
    if {"flight", "ground"}.issubset(set(by_sample["condition"])):
        f = float(by_sample.loc[by_sample["condition"] == "flight", "ercc_fraction"].mean())
        g = float(by_sample.loc[by_sample["condition"] == "ground", "ercc_fraction"].mean())
        summary["flight_minus_ground_mean_fraction"] = f - g
        summary["flight_over_ground_mean_fraction"] = (f / g) if g > 0 else float("nan")
    return summary


def ercc_calibration_frame(de_m: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    """Join M_FLIGHT DE ERCC rows to ExFold expected log2FC."""
    expected = ercc_expected_table(cfg)
    if de_m is None or de_m.empty or "gene_symbol" not in de_m.columns:
        return pd.DataFrame()
    ercc = de_m.loc[de_m["gene_symbol"].map(is_ercc_symbol)].copy()
    if ercc.empty:
        return pd.DataFrame()
    ercc["ercc_id"] = ercc["gene_symbol"].astype(str).str.upper()
    # Some tables may store bare IDs without case consistency
    out = expected.merge(
        ercc,
        on="ercc_id",
        how="left",
        suffixes=("", "_de"),
    )
    out["observed_log2fc"] = pd.to_numeric(out.get("logFC"), errors="coerce")
    out["observed_padj"] = pd.to_numeric(out.get("padj"), errors="coerce")
    out["residual"] = out["observed_log2fc"] - out["expected_log2fc_mix1_over_mix2"]
    out["in_de_table"] = out["observed_log2fc"].notna()
    return out


def ercc_calibration_summary(calib: pd.DataFrame) -> dict[str, Any]:
    """Compact calibration stats for run_params / methods.

    Headline metrics are OLS slope, intercept, and R2 (observed ~ expected).
    Spearman is retained as a secondary rank check.
    """
    if calib is None or calib.empty:
        return {"n_ercc_expected": 0, "n_ercc_observed": 0}
    obs = calib.loc[calib["in_de_table"]].copy()
    stats: dict[str, Any] = {
        "n_ercc_expected": int(len(calib)),
        "n_ercc_observed": int(len(obs)),
        "flight_mix": FLIGHT_MIX,
        "ground_mix": GROUND_MIX,
        "ols_slope": float("nan"),
        "ols_intercept": float("nan"),
        "ols_r2": float("nan"),
        "spearman_obs_vs_exp": float("nan"),
        "pearson_obs_vs_exp": float("nan"),
        "mean_abs_residual": float("nan"),
    }
    if len(obs) >= 3:
        from scipy.stats import pearsonr, spearmanr

        x = obs["expected_log2fc_mix1_over_mix2"].astype(float).to_numpy()
        y = obs["observed_log2fc"].astype(float).to_numpy()
        # OLS: y = a + b x
        b, a = np.polyfit(x, y, 1)
        yhat = a + b * x
        ss_res = float(np.sum((y - yhat) ** 2))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        stats["ols_slope"] = float(b)
        stats["ols_intercept"] = float(a)
        stats["ols_r2"] = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
        stats["spearman_obs_vs_exp"] = float(spearmanr(x, y).correlation)
        stats["pearson_obs_vs_exp"] = float(pearsonr(x, y)[0])
        stats["mean_abs_residual"] = float(np.nanmean(np.abs(obs["residual"])))
    by = (
        obs.groupby("subgroup", observed=True)
        .agg(
            n=("ercc_id", "count"),
            expected=("expected_log2fc_mix1_over_mix2", "first"),
            mean_observed=("observed_log2fc", "mean"),
            mean_residual=("residual", "mean"),
        )
        .reset_index()
    )
    stats["by_subgroup"] = by.to_dict(orient="records")
    return stats


def plot_ercc_calibration(
    calib: pd.DataFrame,
    path: Path | str,
    *,
    summary: dict[str, Any] | None = None,
    dpi: int = 150,
) -> Path | None:
    """Observed vs expected log2FC for ERCC ExFold subgroups (separate from biological volcano)."""
    import matplotlib.pyplot as plt

    if calib is None or calib.empty:
        return None
    obs = calib.loc[calib["in_de_table"]].copy()
    if obs.empty:
        return None
    colors = {"A": "#4C78A8", "B": "#F58518", "C": "#54A24B", "D": "#E45756"}
    fig, ax = plt.subplots(figsize=(5.5, 5.2))
    for sub, sub_df in obs.groupby("subgroup", observed=True):
        ax.scatter(
            sub_df["expected_log2fc_mix1_over_mix2"],
            sub_df["observed_log2fc"],
            s=28,
            alpha=0.85,
            c=colors.get(str(sub), "#666666"),
            label=f"Subgroup {sub} (exp={EXPECTED_BY_SUBGROUP.get(str(sub), '?')})",
            edgecolors="k",
            linewidths=0.3,
        )
    lims = [-1.5, 2.5]
    ax.plot(lims, lims, ls="--", color="#888888", lw=1.0, label="y = x (perfect recovery)")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Expected log2FC (Mix1 / Mix2)")
    ax.set_ylabel("Observed log2FC (flight / GC)")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    rho = (summary or {}).get("spearman_obs_vs_exp")
    slope = (summary or {}).get("ols_slope")
    intercept = (summary or {}).get("ols_intercept")
    r2 = (summary or {}).get("ols_r2")
    n_obs = (summary or {}).get("n_ercc_observed", len(obs))
    if slope == slope and r2 == r2:
        title = (
            "ERCC ExFold calibration (ISS-T: flight=Mix1, GC=Mix2)\n"
            f"n={n_obs}/92; OLS slope={float(slope):.3f}, intercept={float(intercept):.3f}, "
            f"R^2={float(r2):.3f}"
        )
    elif rho == rho:
        title = (
            "ERCC ExFold calibration (ISS-T: flight=Mix1, GC=Mix2)\n"
            f"n={n_obs}/92; Spearman rho~{float(rho):.3f}"
        )
    else:
        title = f"ERCC ExFold calibration (ISS-T: flight=Mix1, GC=Mix2); n={n_obs}/92"
    ax.set_title(title, fontsize=10)
    ax.text(
        0.02,
        0.02,
        "Spike-in mix is confounded with arm by design;\n"
        "ERCC cannot normalize this contrast. Separate from biology volcano.",
        transform=ax.transAxes,
        fontsize=7.5,
        va="bottom",
        color="#333333",
    )
    fig.tight_layout()
    return save_figure(path, dpi=dpi, close=True)
