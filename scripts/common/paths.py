"""Resolve project paths from config/paths.yaml."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


# CFDE_Module_Build/
MODULE_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = MODULE_ROOT / "config" / "paths.yaml"

# Dead decoy keys removed in Stage 1 (P5-1). Raising keeps learners on live knobs.
_FORBIDDEN_MODULE_KEYS: dict[str, frozenset[str]] = {
    "module1": frozenset({"thresholds", "params"}),
    "module2": frozenset({"thresholds", "params"}),
}


def _assert_no_decoy_keys(cfg: dict[str, Any]) -> None:
    for mod, forbidden in _FORBIDDEN_MODULE_KEYS.items():
        block = cfg.get(mod)
        if not isinstance(block, dict):
            continue
        bad = sorted(forbidden.intersection(block))
        if bad:
            raise ValueError(
                f"config/{mod} contains removed decoy key(s) {bad}; "
                f"use live knobs (module1.qc / module2.embedding+de+harmony) instead"
            )
    geo = (cfg.get("module4") or {}).get("geo") or {}
    if isinstance(geo, dict) and "exclude_diagnoses" in geo:
        raise ValueError(
            "config/module4.geo.exclude_diagnoses is unused; remove it "
            "(include_diagnoses is the live filter)"
        )


def load_config(config_path: Path | None = None) -> dict[str, Any]:
    path = Path(config_path) if config_path else CONFIG_PATH
    with path.open() as fh:
        cfg = yaml.safe_load(fh)
    cfg["_module_root"] = MODULE_ROOT
    cfg["_config_path"] = path
    _assert_no_decoy_keys(cfg)
    return cfg


def module_root(cfg: dict[str, Any] | None = None) -> Path:
    if cfg is None:
        cfg = load_config()
    return Path(cfg["_module_root"])


def portable_path(cfg: dict[str, Any] | None, path: Path | str) -> str:
    """Return a repo-relative POSIX path when ``path`` is under the module root.

    Absolute build-machine paths must not ship in learner-facing tables or
    ``*_run_params.json``. Falls back to the original string if the path is
    outside the package root.
    """
    root = module_root(cfg).resolve()
    p = Path(path)
    try:
        resolved = p.resolve()
    except OSError:
        resolved = p
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError:
        text = str(path).replace("\\", "/")
        for marker in ("/scripts/", "/outputs/", "/docs/", "/config/", "/data/"):
            if marker in text:
                return marker.lstrip("/") + text.split(marker, 1)[-1]
        return text


def resolve(cfg: dict[str, Any], key: str) -> Path:
    """Resolve a paths.<key> entry to an absolute Path."""
    rel = cfg["paths"][key]
    return module_root(cfg) / rel


def ensure_output_dirs(cfg: dict[str, Any]) -> None:
    for key in (
        "outputs_tables",
        "outputs_figures",
        "outputs_reports",
        "processed_hubmap",
        "processed_integrated",
    ):
        if key in cfg.get("paths", {}):
            resolve(cfg, key).mkdir(parents=True, exist_ok=True)








def list_donor_keys(cfg: dict[str, Any]) -> list[str]:
    """Return config keys like donor_1 ... donor_N that are present."""
    keys = [
        k
        for k in cfg
        if isinstance(k, str) and k.startswith("donor_") and isinstance(cfg[k], dict)
    ]

    def _num(k: str) -> int:
        try:
            return int(k.split("_", 1)[1])
        except (IndexError, ValueError):
            return 0

    return sorted(keys, key=_num)


def donor_root(cfg: dict[str, Any], donor_label: str | None = None) -> Path:
    """Resolve a donor folder by label (Donor_1 ... Donor_N)."""
    label = donor_label or cfg.get("module1", {}).get("donor_label", "Donor_1")
    normalized = str(label).strip().replace(" ", "_")
    if normalized.lower().startswith("donor_"):
        num = normalized.split("_", 1)[1]
        folder = f"Donor_{num}"
        key = f"donor_{num}_root"
    else:
        folder = str(label)
        key = "donor_1_root"
    if key in cfg.get("paths", {}):
        return resolve(cfg, key)
    return resolve(cfg, "hubmap_source") / folder


def _input_ids_for_donor(
    cfg: dict[str, Any],
    donor_label: str | None,
) -> dict[str, str]:
    """Return primary_id / processed_id / block_id for a donor from module2.inputs or module1."""
    label = donor_label or cfg.get("module1", {}).get("donor_label")
    for item in cfg.get("module2", {}).get("inputs") or []:
        if label and item.get("donor_label") == label:
            return {
                "primary_id": str(item.get("primary_id") or item.get("block_id") or ""),
                "processed_id": str(item.get("processed_id") or ""),
                "block_id": str(item.get("block_id") or item.get("primary_id") or ""),
            }
    m1 = cfg.get("module1") or {}
    return {
        "primary_id": str(m1.get("primary_id") or m1.get("block_id") or ""),
        "processed_id": str(m1.get("processed_id") or ""),
        "block_id": str(m1.get("block_id") or m1.get("primary_id") or ""),
    }


def snrna_block_dir(
    cfg: dict[str, Any],
    block_id: str | None = None,
    donor_label: str | None = None,
) -> Path:
    """Resolve snRNAseq block folder: try processed_id first, then primary_id / block_id."""
    ids = _input_ids_for_donor(cfg, donor_label)
    root = donor_root(cfg, donor_label)
    candidates: list[str] = []
    # Explicit block_id from caller still preferred when it names an existing folder,
    # but processed_id is tried first when resolving from config.
    if block_id:
        # When caller passes a primary/local key, still try matching processed first
        # if this block_id equals the configured primary_id for the donor.
        if ids.get("processed_id") and block_id in {
            ids.get("primary_id"),
            ids.get("block_id"),
            ids.get("processed_id"),
        }:
            candidates.append(ids["processed_id"])
        candidates.append(block_id)
    else:
        if ids.get("processed_id"):
            candidates.append(ids["processed_id"])
        if ids.get("primary_id"):
            candidates.append(ids["primary_id"])
        if ids.get("block_id"):
            candidates.append(ids["block_id"])
        m1_block = (cfg.get("module1") or {}).get("block_id")
        if m1_block:
            candidates.append(str(m1_block))

    seen: set[str] = set()
    ordered: list[str] = []
    for c in candidates:
        if c and c not in seen:
            seen.add(c)
            ordered.append(c)

    for bid in ordered:
        path = root / "snRNAseq" / bid
        if path.exists():
            return path
    # Fall back to first candidate even if missing (callers report FileNotFoundError)
    fallback = ordered[0] if ordered else (block_id or ids.get("block_id") or "UNKNOWN")
    return root / "snRNAseq" / fallback


def snatac_dir(cfg: dict[str, Any], donor_label: str = "Donor_1") -> Path:
    return donor_root(cfg, donor_label) / "snATACseq"
