"""Shared-mechanism concordance for Module 4 (replaces global rho as headline).

Three biology pairs plus the three-way intersection. First-run pathways ship as
tier 3 (concordance only) with leading-edge overlap columns; tier 1/2 assignment
waits until the Jaccard distribution is inspected. Suspect ribosome/translation/
OxPhos terms are flagged, never excluded.
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from scripts.common.paths import module_root, resolve
from scripts.nb06_crossres.orthologs import load_ortholog_table

BIOLOGY_PAIRS: list[tuple[str, str, str]] = [
    ("H_DISEASE", "H_AGING_GTEX", "disease_aging"),
    ("H_DISEASE", "M_FLIGHT", "disease_flight"),
    ("H_AGING_GTEX", "M_FLIGHT", "aging_flight"),
]

CROSS_SPECIES_CONTRASTS = {"M_FLIGHT"}

SUSPECT_PATTERNS = [
    re.compile(p, re.I)
    for p in [
        r"\bribosom",
        r"\btranslation\b",
        r"\btranslational\b",
        r"oxidative phosphorylation",
        r"\boxphos\b",
        r"electron transport",
        r"respiratory chain",
        r"\beif[0-9]",
        r"peptide chain elongation",
        r"nonsense.mediated decay",
        r"\brrna\b",
        r"cytosolic ribosom",
        r"mitochondrial translation",
        r"selenoamino acid metabolism",
    ]
]


def _strip_term(term: str) -> str:
    t = str(term)
    return re.sub(r"^.*\.gmt__", "", t)


def _parse_lead(raw: object) -> set[str]:
    if raw is None or (isinstance(raw, float) and np.isnan(raw)):
        return set()
    s = str(raw).strip()
    if not s or s.lower() in {"nan", "none"}:
        return set()
    return {g.strip().upper() for g in s.split(";") if g.strip()}


def _is_suspect(term: str) -> bool:
    return any(p.search(term) for p in SUSPECT_PATTERNS)


def load_nes_table(cfg: dict[str, Any], contrast_id: str) -> pd.DataFrame:
    tables = resolve(cfg, "outputs_tables")
    path = tables / f"module4_nes_{contrast_id}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Missing NES table for {contrast_id}: {path}")
    df = pd.read_csv(path, sep="\t")
    df = df.copy()
    df["Term"] = df["Term"].astype(str).map(_strip_term)
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    fdr_col = "FDR q-val" if "FDR q-val" in df.columns else "FDR"
    df["FDR"] = pd.to_numeric(df[fdr_col], errors="coerce")
    df["lead_genes"] = df["Lead_genes"].map(_parse_lead) if "Lead_genes" in df.columns else [set()] * len(df)
    df["direction"] = np.sign(df["NES"]).astype(int)
    df = df.dropna(subset=["NES", "FDR"]).drop_duplicates(subset=["Term"], keep="first")
    return df


def mechanism_params(cfg: dict[str, Any]) -> dict[str, Any]:
    m = ((cfg.get("module4") or {}).get("mechanisms") or {})
    return {
        "fdr_q": float(m.get("fdr_q", 0.25)),
        "nes_abs_min": float(m.get("nes_abs_min", 1.5)),
        "fdr_q_strict": float(m.get("fdr_q_strict", 0.05)),
        "cluster_jaccard_min": float(m.get("cluster_jaccard_min", 0.5)),
        "n_permutations": int(m.get("n_permutations", 2000)),
        "random_seed": int(m.get("random_seed", 0)),
        "top_intersection_genes": int(m.get("top_intersection_genes", 40)),
    }


def pass_mask(df: pd.DataFrame, *, fdr_q: float, nes_abs_min: float) -> pd.Series:
    return (df["FDR"] < fdr_q) & (df["NES"].abs() >= nes_abs_min)


def ortholog_human_universe(cfg: dict[str, Any]) -> set[str]:
    ortho = load_ortholog_table(cfg)
    return set(ortho["human_symbol"].astype(str).str.upper())


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return float("nan")
    u = len(a | b)
    return float(len(a & b) / u) if u else float("nan")


def leading_edge_overlap(
    lead_a: set[str],
    lead_b: set[str],
    *,
    restrict_to: set[str] | None,
    cross_species: bool,
) -> dict[str, Any]:
    """Symmetric LE overlap. For cross-species pairs, restrict both sides first."""
    a_raw, b_raw = set(lead_a), set(lead_b)
    if cross_species and restrict_to is not None:
        a_use = a_raw & restrict_to
        b_use = b_raw & restrict_to
    else:
        a_use, b_use = a_raw, b_raw
    inter = a_use & b_use
    return {
        "le_a_n_raw": int(len(a_raw)),
        "le_b_n_raw": int(len(b_raw)),
        "le_a_n_restricted": int(len(a_use)),
        "le_b_n_restricted": int(len(b_use)),
        "le_intersection_n": int(len(inter)),
        "le_jaccard": _jaccard(a_use, b_use),
        "le_intersection_genes": ";".join(sorted(inter)),
        "le_restricted_to_ortholog_mappable": bool(cross_species),
    }


def analytic_null(
    n_shared: int,
    k_a: int,
    k_b: int,
    p_a_up: float,
    p_b_up: float,
) -> dict[str, float]:
    if n_shared <= 0:
        return {
            "E_both_pass": float("nan"),
            "P_agree": float("nan"),
            "E_concordant": float("nan"),
        }
    e_both = (k_a * k_b) / n_shared
    p_agree = p_a_up * p_b_up + (1.0 - p_a_up) * (1.0 - p_b_up)
    return {
        "E_both_pass": float(e_both),
        "P_agree": float(p_agree),
        "E_concordant": float(e_both * p_agree),
    }


def permute_concordant_count(
    pass_a: np.ndarray,
    sign_a: np.ndarray,
    pass_b: np.ndarray,
    sign_b: np.ndarray,
    *,
    n_perm: int,
    seed: int,
) -> dict[str, float]:
    """
    Shuffle pathway labels on contrast B (pass + sign together), keep marginals.

    Returns mean/median expected concordant and a one-sided p (obs >= null).
    Caller supplies observed count separately for p-value.
    """
    rng = np.random.default_rng(seed)
    n = len(pass_a)
    if n == 0:
        return {"perm_mean": float("nan"), "perm_median": float("nan"), "perm_values": []}
    # Bundle B's labels and permute assignment onto fixed pathway order of A
    b_pack = np.stack([pass_b.astype(int), sign_b.astype(int)], axis=1)
    counts = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        shuffled = b_pack[rng.permutation(n)]
        both = pass_a & (shuffled[:, 0].astype(bool))
        agree = sign_a[both] == shuffled[both, 1]
        counts[i] = float(agree.sum())
    return {
        "perm_mean": float(np.mean(counts)),
        "perm_median": float(np.median(counts)),
        "perm_sd": float(np.std(counts, ddof=1)) if n_perm > 1 else float("nan"),
        "perm_values": counts,
    }


def analytic_null_three_way(
    n_shared: int,
    k: dict[str, int],
    p_up: dict[str, float],
    ids: tuple[str, str, str],
) -> dict[str, float]:
    """
    Expected all-three-pass with all-agree signs.

    Eight sign patterns; two are all-agree (all +, all -).
    P(all agree) = p1*p2*p3 + (1-p1)*(1-p2)*(1-p3).
    Independence assumption; permutation is the quoteable null.
    """
    a, b, c = ids
    if n_shared <= 0:
        return {"E_all_pass": float("nan"), "P_all_agree": float("nan"), "E_concordant": float("nan")}
    e_all = (k[a] * k[b] * k[c]) / (n_shared ** 2)
    p_agree = p_up[a] * p_up[b] * p_up[c] + (1 - p_up[a]) * (1 - p_up[b]) * (1 - p_up[c])
    return {
        "E_all_pass": float(e_all),
        "P_all_agree": float(p_agree),
        "E_concordant": float(e_all * p_agree),
    }


def permute_three_way_concordant(
    frames: dict[str, pd.DataFrame],
    terms: list[str],
    ids: tuple[str, str, str],
    *,
    fdr_q: float,
    nes_abs_min: float,
    n_perm: int,
    seed: int,
) -> dict[str, float]:
    """Permute pathway labels independently for contrasts B and C."""
    rng = np.random.default_rng(seed)
    a, b, c = ids
    n = len(terms)
    packs = {}
    for cid in ids:
        sub = frames[cid].set_index("Term").reindex(terms)
        passed = pass_mask(sub.fillna({"FDR": 1.0, "NES": 0.0}), fdr_q=fdr_q, nes_abs_min=nes_abs_min).to_numpy()
        signs = np.sign(pd.to_numeric(sub["NES"], errors="coerce").fillna(0).to_numpy()).astype(int)
        packs[cid] = np.stack([passed.astype(int), signs], axis=1)

    fixed = packs[a]
    counts = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        pb = packs[b][rng.permutation(n)]
        pc = packs[c][rng.permutation(n)]
        both = fixed[:, 0].astype(bool) & pb[:, 0].astype(bool) & pc[:, 0].astype(bool)
        if not both.any():
            counts[i] = 0.0
            continue
        sa, sb, sc = fixed[both, 1], pb[both, 1], pc[both, 1]
        counts[i] = float(((sa == sb) & (sa == sc)).sum())
    return {
        "perm_mean": float(np.mean(counts)),
        "perm_median": float(np.median(counts)),
        "perm_sd": float(np.std(counts, ddof=1)) if n_perm > 1 else float("nan"),
        "perm_values": counts,
    }


def _pair_concordant_rows(
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    *,
    contrast_a: str,
    contrast_b: str,
    pair_id: str,
    fdr_q: float,
    nes_abs_min: float,
    fdr_q_strict: float,
    mappable: set[str],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    merged = df_a.merge(df_b, on="Term", how="inner", suffixes=("_a", "_b"))
    n_shared = int(len(merged))
    pass_a = pass_mask(
        merged.rename(columns={"FDR_a": "FDR", "NES_a": "NES"}),
        fdr_q=fdr_q,
        nes_abs_min=nes_abs_min,
    )
    pass_b = pass_mask(
        merged.rename(columns={"FDR_b": "FDR", "NES_b": "NES"}),
        fdr_q=fdr_q,
        nes_abs_min=nes_abs_min,
    )
    # rebuild cleanly
    pass_a = (merged["FDR_a"] < fdr_q) & (merged["NES_a"].abs() >= nes_abs_min)
    pass_b = (merged["FDR_b"] < fdr_q) & (merged["NES_b"].abs() >= nes_abs_min)
    k_a, k_b = int(pass_a.sum()), int(pass_b.sum())
    if k_a > 0:
        p_a_up = float((merged.loc[pass_a, "direction_a"] > 0).mean())
    else:
        p_a_up = float("nan")
    if k_b > 0:
        p_b_up = float((merged.loc[pass_b, "direction_b"] > 0).mean())
    else:
        p_b_up = float("nan")

    both = pass_a & pass_b
    agree = merged["direction_a"] == merged["direction_b"]
    conc = merged.loc[both & agree].copy()

    cross = contrast_a in CROSS_SPECIES_CONTRASTS or contrast_b in CROSS_SPECIES_CONTRASTS
    rows = []
    for _, r in conc.iterrows():
        ov = leading_edge_overlap(
            r["lead_genes_a"],
            r["lead_genes_b"],
            restrict_to=mappable if cross else None,
            cross_species=cross,
        )
        strict = (
            (float(r["FDR_a"]) < fdr_q_strict)
            and (float(r["FDR_b"]) < fdr_q_strict)
            and (abs(float(r["NES_a"])) >= nes_abs_min)
            and (abs(float(r["NES_b"])) >= nes_abs_min)
        )
        rows.append(
            {
                "pair_id": pair_id,
                "contrast_a": contrast_a,
                "contrast_b": contrast_b,
                "Term": r["Term"],
                "NES_a": float(r["NES_a"]),
                "NES_b": float(r["NES_b"]),
                "FDR_a": float(r["FDR_a"]),
                "FDR_b": float(r["FDR_b"]),
                "direction": "up" if int(r["direction_a"]) > 0 else "down",
                "tier": 3,
                "passes_fdr_0.05_both": bool(strict),
                "suspect_flag": _is_suspect(str(r["Term"])),
                **ov,
            }
        )
    conc_df = pd.DataFrame(rows)
    obs = int(len(conc_df))
    meta = {
        "pair_id": pair_id,
        "contrast_a": contrast_a,
        "contrast_b": contrast_b,
        "N_shared_pathways": n_shared,
        "k_a": k_a,
        "k_b": k_b,
        "p_a_up": p_a_up,
        "p_b_up": p_b_up,
        "observed_concordant": obs,
        "observed_concordant_unflagged": int((~conc_df["suspect_flag"]).sum()) if obs else 0,
        "observed_concordant_fdr05": int(conc_df["passes_fdr_0.05_both"].sum()) if obs else 0,
        "pass_a": pass_a.to_numpy(),
        "pass_b": pass_b.to_numpy(),
        "sign_a": merged["direction_a"].to_numpy().astype(int),
        "sign_b": merged["direction_b"].to_numpy().astype(int),
        "terms": merged["Term"].tolist(),
    }
    meta.update(analytic_null(n_shared, k_a, k_b, p_a_up if p_a_up == p_a_up else 0.5, p_b_up if p_b_up == p_b_up else 0.5))
    return conc_df, meta


def cluster_mechanisms(
    conc: pd.DataFrame,
    *,
    jaccard_min: float,
    top_genes: int,
) -> pd.DataFrame:
    """Connected components on LE Jaccard >= threshold within each pair_id."""
    if conc is None or conc.empty:
        return pd.DataFrame()

    out_rows = []
    for pair_id, sub in conc.groupby("pair_id", sort=False):
        sub = sub.reset_index(drop=True)
        n = len(sub)
        leads = []
        for g in sub["le_intersection_genes"].fillna(""):
            leads.append({x for x in str(g).split(";") if x} if g else set())
        # Prefer restricted LE sets reconstructed from intersection + sizes is lossy;
        # cluster on intersection genes union individually stored? Use intersection for
        # within-pair similarity of *shared* biology; better: recompute from raw columns.
        # We stored only intersection string. Rebuild pairwise Jaccard from that is wrong.
        # Fix: store restricted LE sets during concordant build  -  for clustering use
        # intersection genes as nodes' features is weak. Recompute from le_intersection
        # is insufficient. Use the intersection gene sets for clustering of pathways
        # that share the same LE overlap genes  -  actually doc says Jaccard on leading-edge
        # genes (the restricted LEs), not on intersection.
        #
        # We need restricted LE sets. Re-parse from columns: we have le_a_n_restricted but
        # not the gene lists. Add them in concordant rows going forward.
        # For now if columns le_a_genes_restricted / le_b exist use them.
        if "le_a_genes_restricted" in sub.columns and "le_b_genes_restricted" in sub.columns:
            sets_a = [set(str(x).split(";")) - {""} for x in sub["le_a_genes_restricted"].fillna("")]
            sets_b = [set(str(x).split(";")) - {""} for x in sub["le_b_genes_restricted"].fillna("")]
            node_sets = [sa | sb for sa, sb in zip(sets_a, sets_b)]
        else:
            node_sets = leads  # fallback: intersection genes only

        parent = list(range(n))

        def find(i: int) -> int:
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        def union(i: int, j: int) -> None:
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[rj] = ri

        for i in range(n):
            for j in range(i + 1, n):
                if _jaccard(node_sets[i], node_sets[j]) >= jaccard_min:
                    union(i, j)

        clusters: dict[int, list[int]] = defaultdict(list)
        for i in range(n):
            clusters[find(i)].append(i)

        for cid, idxs in enumerate(sorted(clusters.values(), key=lambda ix: -len(ix))):
            members = sub.iloc[idxs]
            # Representative: largest |NES| mean
            score = members["NES_a"].abs() + members["NES_b"].abs()
            rep = members.loc[score.idxmax(), "Term"]
            union_genes = set()
            for g in members["le_intersection_genes"].fillna(""):
                union_genes |= {x for x in str(g).split(";") if x}
            flagged = bool(members["suspect_flag"].any())
            all_flagged = bool(members["suspect_flag"].all())
            out_rows.append(
                {
                    "pair_id": pair_id,
                    "cluster_id": f"{pair_id}__c{cid+1:02d}",
                    "n_member_terms": int(len(members)),
                    "representative_term": rep,
                    "member_terms": "; ".join(members["Term"].astype(str).tolist()),
                    "direction": members["direction"].iloc[0],
                    "suspect_flag_any": flagged,
                    "suspect_flag_all": all_flagged,
                    "shared_leading_edge_n": int(len(union_genes)),
                    "shared_leading_edge_genes": ";".join(sorted(union_genes)[:top_genes]),
                    "tier": 3,
                }
            )
    return pd.DataFrame(out_rows)


def run_mechanism_analysis(cfg: dict[str, Any]) -> dict[str, pd.DataFrame]:
    """Build concordant / null / mechanism tables from existing NES outputs."""
    params = mechanism_params(cfg)
    mappable = ortholog_human_universe(cfg)
    frames = {
        cid: load_nes_table(cfg, cid)
        for cid in ["H_DISEASE", "H_AGING_GTEX", "M_FLIGHT"]
    }

    conc_parts: list[pd.DataFrame] = []
    null_rows: list[dict[str, Any]] = []

    for ca, cb, pair_id in BIOLOGY_PAIRS:
        # Enrich overlap columns with restricted gene lists for clustering
        df_a = frames[ca]
        df_b = frames[cb]
        # Temporarily monkey-patch overlap to also return gene lists  -  inline here
        merged = df_a.merge(df_b, on="Term", how="inner", suffixes=("_a", "_b"))
        pass_a = (merged["FDR_a"] < params["fdr_q"]) & (merged["NES_a"].abs() >= params["nes_abs_min"])
        pass_b = (merged["FDR_b"] < params["fdr_q"]) & (merged["NES_b"].abs() >= params["nes_abs_min"])
        k_a, k_b = int(pass_a.sum()), int(pass_b.sum())
        p_a_up = float((merged.loc[pass_a, "direction_a"] > 0).mean()) if k_a else float("nan")
        p_b_up = float((merged.loc[pass_b, "direction_b"] > 0).mean()) if k_b else float("nan")
        both = pass_a & pass_b
        agree = merged["direction_a"] == merged["direction_b"]
        conc = merged.loc[both & agree].copy()
        cross = ca in CROSS_SPECIES_CONTRASTS or cb in CROSS_SPECIES_CONTRASTS
        rows = []
        for _, r in conc.iterrows():
            a_raw = set(r["lead_genes_a"])
            b_raw = set(r["lead_genes_b"])
            if cross:
                a_use, b_use = a_raw & mappable, b_raw & mappable
            else:
                a_use, b_use = a_raw, b_raw
            inter = a_use & b_use
            strict = (
                float(r["FDR_a"]) < params["fdr_q_strict"]
                and float(r["FDR_b"]) < params["fdr_q_strict"]
            )
            rows.append(
                {
                    "pair_id": pair_id,
                    "contrast_a": ca,
                    "contrast_b": cb,
                    "Term": r["Term"],
                    "NES_a": float(r["NES_a"]),
                    "NES_b": float(r["NES_b"]),
                    "FDR_a": float(r["FDR_a"]),
                    "FDR_b": float(r["FDR_b"]),
                    "direction": "up" if int(r["direction_a"]) > 0 else "down",
                    "tier": 3,
                    "passes_fdr_0.05_both": bool(strict),
                    "suspect_flag": _is_suspect(str(r["Term"])),
                    "le_a_n_raw": len(a_raw),
                    "le_b_n_raw": len(b_raw),
                    "le_a_n_restricted": len(a_use),
                    "le_b_n_restricted": len(b_use),
                    "le_intersection_n": len(inter),
                    "le_jaccard": _jaccard(a_use, b_use),
                    "le_intersection_genes": ";".join(sorted(inter)),
                    "le_a_genes_restricted": ";".join(sorted(a_use)),
                    "le_b_genes_restricted": ";".join(sorted(b_use)),
                    "le_restricted_to_ortholog_mappable": bool(cross),
                    "note_m_flight_le_already_human": bool(
                        ca == "M_FLIGHT" or cb == "M_FLIGHT"
                    ),
                }
            )
        conc_df = pd.DataFrame(rows)
        conc_parts.append(conc_df)
        obs = int(len(conc_df))
        an = analytic_null(
            int(len(merged)),
            k_a,
            k_b,
            p_a_up if p_a_up == p_a_up else 0.5,
            p_b_up if p_b_up == p_b_up else 0.5,
        )
        perm = permute_concordant_count(
            pass_a.to_numpy(),
            merged["direction_a"].to_numpy().astype(int),
            pass_b.to_numpy(),
            merged["direction_b"].to_numpy().astype(int),
            n_perm=params["n_permutations"],
            seed=params["random_seed"],
        )
        pv = perm.pop("perm_values")
        perm_p = float(np.mean(pv >= obs)) if len(pv) else float("nan")
        null_rows.append(
            {
                "set_id": pair_id,
                "kind": "pairwise",
                "contrast_a": ca,
                "contrast_b": cb,
                "contrast_c": "",
                "N_shared_pathways": int(len(merged)),
                "k_a": k_a,
                "k_b": k_b,
                "k_c": "",
                "p_a_up": p_a_up,
                "p_b_up": p_b_up,
                "p_c_up": "",
                "observed_concordant": obs,
                "observed_concordant_unflagged": int((~conc_df["suspect_flag"]).sum())
                if obs
                else 0,
                "observed_concordant_fdr05": int(conc_df["passes_fdr_0.05_both"].sum())
                if obs
                else 0,
                **an,
                "analytic_ratio_obs_over_E": (
                    obs / an["E_concordant"] if an["E_concordant"] not in (0, float("nan")) and an["E_concordant"] == an["E_concordant"] and an["E_concordant"] > 0 else float("nan")
                ),
                "perm_mean": perm["perm_mean"],
                "perm_median": perm["perm_median"],
                "perm_sd": perm["perm_sd"],
                "perm_ratio_obs_over_mean": (
                    obs / perm["perm_mean"] if perm["perm_mean"] and perm["perm_mean"] == perm["perm_mean"] and perm["perm_mean"] > 0 else float("nan")
                ),
                "perm_p_ge_obs": perm_p,
                "n_permutations": params["n_permutations"],
                "fdr_q": params["fdr_q"],
                "nes_abs_min": params["nes_abs_min"],
            }
        )

    # --- three-way ---
    ids = ("H_DISEASE", "H_AGING_GTEX", "M_FLIGHT")
    m3 = frames[ids[0]][["Term", "NES", "FDR", "direction", "lead_genes"]].rename(
        columns={c: f"{c}_0" if c != "Term" else c for c in ["NES", "FDR", "direction", "lead_genes"]}
    )
    # cleaner merge
    t = frames[ids[0]][["Term"]].copy()
    for i, cid in enumerate(ids):
        part = frames[cid][["Term", "NES", "FDR", "direction", "lead_genes"]].rename(
            columns={
                "NES": f"NES_{i}",
                "FDR": f"FDR_{i}",
                "direction": f"direction_{i}",
                "lead_genes": f"lead_genes_{i}",
            }
        )
        t = t.merge(part, on="Term", how="inner")
    n_shared3 = int(len(t))
    passes = []
    ks = {}
    p_ups = {}
    for i, cid in enumerate(ids):
        p = (t[f"FDR_{i}"] < params["fdr_q"]) & (t[f"NES_{i}"].abs() >= params["nes_abs_min"])
        passes.append(p)
        ks[cid] = int(p.sum())
        p_ups[cid] = float((t.loc[p, f"direction_{i}"] > 0).mean()) if ks[cid] else float("nan")
    all_pass = passes[0] & passes[1] & passes[2]
    same_sign = (t["direction_0"] == t["direction_1"]) & (t["direction_0"] == t["direction_2"])
    three = t.loc[all_pass & same_sign].copy()
    three_rows = []
    for _, r in three.iterrows():
        # Restrict all three LEs to mappable (flight is cross-species)
        sets = []
        for i in range(3):
            raw = set(r[f"lead_genes_{i}"])
            sets.append(raw & mappable)
        inter = sets[0] & sets[1] & sets[2]
        three_rows.append(
            {
                "pair_id": "three_way",
                "contrast_a": ids[0],
                "contrast_b": ids[1],
                "contrast_c": ids[2],
                "Term": r["Term"],
                "NES_a": float(r["NES_0"]),
                "NES_b": float(r["NES_1"]),
                "NES_c": float(r["NES_2"]),
                "FDR_a": float(r["FDR_0"]),
                "FDR_b": float(r["FDR_1"]),
                "FDR_c": float(r["FDR_2"]),
                "direction": "up" if int(r["direction_0"]) > 0 else "down",
                "tier": 3,
                "passes_fdr_0.05_both": bool(
                    float(r["FDR_0"]) < params["fdr_q_strict"]
                    and float(r["FDR_1"]) < params["fdr_q_strict"]
                    and float(r["FDR_2"]) < params["fdr_q_strict"]
                ),
                "suspect_flag": _is_suspect(str(r["Term"])),
                "le_a_n_raw": len(r["lead_genes_0"]),
                "le_b_n_raw": len(r["lead_genes_1"]),
                "le_c_n_raw": len(r["lead_genes_2"]),
                "le_a_n_restricted": len(sets[0]),
                "le_b_n_restricted": len(sets[1]),
                "le_c_n_restricted": len(sets[2]),
                "le_intersection_n": len(inter),
                "le_jaccard": (
                    float(len(inter) / len(sets[0] | sets[1] | sets[2]))
                    if (sets[0] | sets[1] | sets[2])
                    else float("nan")
                ),
                "le_intersection_genes": ";".join(sorted(inter)),
                "le_a_genes_restricted": ";".join(sorted(sets[0])),
                "le_b_genes_restricted": ";".join(sorted(sets[1])),
                "le_c_genes_restricted": ";".join(sorted(sets[2])),
                "le_restricted_to_ortholog_mappable": True,
            }
        )
    three_df = pd.DataFrame(three_rows)
    conc_parts.append(three_df)
    obs3 = int(len(three_df))
    an3 = analytic_null_three_way(n_shared3, ks, {k: (v if v == v else 0.5) for k, v in p_ups.items()}, ids)
    # Build frames for permute_three_way - need Term-indexed with FDR/NES
    perm3 = permute_three_way_concordant(
        frames,
        t["Term"].tolist(),
        ids,
        fdr_q=params["fdr_q"],
        nes_abs_min=params["nes_abs_min"],
        n_perm=params["n_permutations"],
        seed=params["random_seed"] + 1,
    )
    pv3 = perm3.pop("perm_values")
    null_rows.append(
        {
            "set_id": "three_way",
            "kind": "three_way",
            "contrast_a": ids[0],
            "contrast_b": ids[1],
            "contrast_c": ids[2],
            "N_shared_pathways": n_shared3,
            "k_a": ks[ids[0]],
            "k_b": ks[ids[1]],
            "k_c": ks[ids[2]],
            "p_a_up": p_ups[ids[0]],
            "p_b_up": p_ups[ids[1]],
            "p_c_up": p_ups[ids[2]],
            "observed_concordant": obs3,
            "observed_concordant_unflagged": int((~three_df["suspect_flag"]).sum()) if obs3 else 0,
            "observed_concordant_fdr05": int(three_df["passes_fdr_0.05_both"].sum()) if obs3 else 0,
            "E_both_pass": an3["E_all_pass"],
            "P_agree": an3["P_all_agree"],
            "E_concordant": an3["E_concordant"],
            "analytic_ratio_obs_over_E": (
                obs3 / an3["E_concordant"]
                if an3["E_concordant"] and an3["E_concordant"] == an3["E_concordant"] and an3["E_concordant"] > 0
                else float("nan")
            ),
            "perm_mean": perm3["perm_mean"],
            "perm_median": perm3["perm_median"],
            "perm_sd": perm3["perm_sd"],
            "perm_ratio_obs_over_mean": (
                obs3 / perm3["perm_mean"]
                if perm3["perm_mean"] and perm3["perm_mean"] == perm3["perm_mean"] and perm3["perm_mean"] > 0
                else float("nan")
            ),
            "perm_p_ge_obs": float(np.mean(pv3 >= obs3)) if len(pv3) else float("nan"),
            "n_permutations": params["n_permutations"],
            "fdr_q": params["fdr_q"],
            "nes_abs_min": params["nes_abs_min"],
        }
    )

    concordant = pd.concat(conc_parts, ignore_index=True) if conc_parts else pd.DataFrame()
    # Cluster pairwise only (not three-way) for mechanism table; also cluster three-way separately
    mech_parts = []
    if not concordant.empty:
        pairwise = concordant[concordant["pair_id"] != "three_way"]
        if not pairwise.empty:
            mech_parts.append(
                cluster_mechanisms(
                    pairwise,
                    jaccard_min=params["cluster_jaccard_min"],
                    top_genes=params["top_intersection_genes"],
                )
            )
        if not three_df.empty:
            # Adapt three-way to cluster_mechanisms shape (needs NES_a/NES_b)
            tw = three_df.copy()
            tw["NES_b"] = tw["NES_b"]  # already
            mech_parts.append(
                cluster_mechanisms(
                    tw,
                    jaccard_min=params["cluster_jaccard_min"],
                    top_genes=params["top_intersection_genes"],
                )
            )
    mechanisms = pd.concat(mech_parts, ignore_index=True) if mech_parts else pd.DataFrame()
    null_df = pd.DataFrame(null_rows)

    # Summary sentence helpers
    summary_rows = []
    for _, r in null_df.iterrows():
        summary_rows.append(
            {
                "set_id": r["set_id"],
                "n_mechanisms_total": (
                    int(mechanisms.loc[mechanisms["pair_id"] == r["set_id"]].shape[0])
                    if not mechanisms.empty
                    else 0
                ),
                "n_mechanisms_unflagged": (
                    int(
                        (
                            ~mechanisms.loc[mechanisms["pair_id"] == r["set_id"], "suspect_flag_any"]
                        ).sum()
                    )
                    if not mechanisms.empty and r["set_id"] in set(mechanisms["pair_id"])
                    else 0
                ),
                "n_concordant_pathways": int(r["observed_concordant"]),
                "n_concordant_unflagged": int(r["observed_concordant_unflagged"]),
                "perm_mean": r["perm_mean"],
                "perm_p_ge_obs": r["perm_p_ge_obs"],
            }
        )
    summary = pd.DataFrame(summary_rows)

    return {
        "concordant_pathways": concordant,
        "concordance_null": null_df,
        "mechanisms": mechanisms,
        "mechanism_summary": summary,
    }


def save_mechanism_outputs(cfg: dict[str, Any], tables: dict[str, pd.DataFrame]) -> dict[str, Path]:
    out_dir = resolve(cfg, "outputs_tables")
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "concordant_pathways": out_dir / "module4_concordant_pathways.tsv",
        "concordance_null": out_dir / "module4_concordance_null.tsv",
        "mechanisms": out_dir / "module4_mechanisms.tsv",
        "mechanism_summary": out_dir / "module4_mechanism_summary.tsv",
    }
    for key, path in paths.items():
        df = tables.get(key, pd.DataFrame())
        # Drop huge set columns before write if present as Python sets (already strings)
        df.to_csv(path, sep="\t", index=False)
    return paths


def plot_nes_scatter_with_concordant(
    nes_joined: pd.DataFrame,
    concordant_terms: set[str],
    path: Path | str,
    *,
    stats: dict[str, Any] | None = None,
    pair_label: str | None = None,
    dpi: int = 150,
) -> Path | None:
    """NES scatter with concordant pathways highlighted."""
    from scripts.nb06_crossres.compare import plot_nes_scatter

    # Reuse base plot by temporarily marking highlights via Term list
    return plot_nes_scatter(
        nes_joined,
        path,
        stats=stats,
        highlight=sorted(concordant_terms)[:12],
        pair_label=pair_label,
        dpi=dpi,
    )
