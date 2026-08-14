"""One-time generator for the Module 3 Azimuth / WNN label cache.

Reads ``secondary_analysis.h5mu`` with h5py in mode ``r`` and never touches an
``X`` matrix. Output is committed so learners do not download the MuData file.

Run from the package root::

    python -m scripts.nb05_multimodal.cache
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd

CACHE_COLUMNS = [
    "barcode",
    "leiden_wnn",
    "azimuth_label",
    "azimuth_id",
    "prediction_score",
    "wnn_umap_1",
    "wnn_umap_2",
]


def _decode_arr(values) -> np.ndarray:
    out = []
    for v in values:
        if isinstance(v, bytes):
            out.append(v.decode("utf-8", errors="replace"))
        else:
            out.append(str(v))
    return np.asarray(out, dtype=object)


def _read_categorical_or_string(group: h5py.Group, key: str) -> np.ndarray:
    if key not in group:
        raise KeyError(f"Missing {group.name}/{key}")
    obj = group[key]
    if isinstance(obj, h5py.Dataset):
        vals = obj[:]
        if vals.dtype.kind in ("S", "O"):
            return _decode_arr(vals)
        return np.asarray(vals)
    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        cats = _decode_arr(obj["categories"][:])
        codes = np.asarray(obj["codes"][:])
        return np.array([cats[c] if c >= 0 else None for c in codes], dtype=object)
    raise TypeError(f"Cannot decode {group.name}/{key}")


def _sha256_file(path: Path, chunk: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def _mofa_sample_list(mofa_path: Path) -> np.ndarray:
    with h5py.File(mofa_path, "r") as f:
        return _decode_arr(f["samples/group1"][:]).astype(str)


def build_label_cache(h5mu_path: Path, out_tsv: Path, dataset_id: str, mofa_path: Path) -> Path:
    """Write the committed Module 3 label cache aligned to MOFA sample order."""
    h5mu_path = Path(h5mu_path)
    out_tsv = Path(out_tsv)
    mofa_path = Path(mofa_path)
    if not h5mu_path.exists():
        raise FileNotFoundError(f"MuData required to build the label cache: {h5mu_path}")
    if not mofa_path.exists():
        raise FileNotFoundError(f"MOFA file required to align the cache: {mofa_path}")

    mofa_samples = _mofa_sample_list(mofa_path)

    with h5py.File(h5mu_path, "r") as f:
        joint_bc = _decode_arr(f["obs/_index"][:]).astype(str)
        rna_bc = _decode_arr(f["mod/rna/obs/_index"][:]).astype(str)
        leiden_wnn = _read_categorical_or_string(f["obs"], "leiden_wnn")
        azimuth_label = _read_categorical_or_string(f["mod/rna/obs"], "azimuth_label")
        azimuth_id = _read_categorical_or_string(f["mod/rna/obs"], "azimuth_id")
        prediction_score = np.asarray(f["mod/rna/obs/prediction_score"][:], dtype=np.float64)
        umap = np.asarray(f["obsm/X_umap"][:], dtype=np.float64)

    joint_set = set(joint_bc.tolist())
    mofa_set = set(mofa_samples.tolist())
    if joint_set != mofa_set:
        only_mofa = sorted(mofa_set - joint_set)[:8]
        only_h5mu = sorted(joint_set - mofa_set)[:8]
        raise ValueError(
            "MOFA samples/group1 and h5mu obs/_index barcode sets differ. "
            f"only_mofa={only_mofa} only_h5mu={only_h5mu}"
        )

    rna_pos = {b: i for i, b in enumerate(rna_bc)}
    joint_pos = {b: i for i, b in enumerate(joint_bc)}
    missing_rna = [b for b in mofa_samples if b not in rna_pos]
    if missing_rna:
        raise ValueError(
            f"RNA obs is missing {len(missing_rna)} MOFA barcodes "
            f"(example {missing_rna[:5]})"
        )

    # Cache rows follow MOFA sample order so the notebook can join positionally.
    order_joint = np.array([joint_pos[b] for b in mofa_samples])
    order_rna = np.array([rna_pos[b] for b in mofa_samples])
    if not np.array_equal(joint_bc, mofa_samples):
        print(
            "cache.py: h5mu obs/_index order differs from MOFA samples; "
            "writing cache in MOFA sample order"
        )

    frame = pd.DataFrame(
        {
            "barcode": mofa_samples,
            "leiden_wnn": leiden_wnn[order_joint],
            "azimuth_label": azimuth_label[order_rna],
            "azimuth_id": azimuth_id[order_rna],
            "prediction_score": prediction_score[order_rna],
            "wnn_umap_1": umap[order_joint, 0],
            "wnn_umap_2": umap[order_joint, 1],
        }
    )
    if list(frame.columns) != CACHE_COLUMNS:
        raise RuntimeError(f"Cache columns {list(frame.columns)} != {CACHE_COLUMNS}")
    if len(frame) != 21978:
        print(f"cache.py: row count is {len(frame)} (probe expected 21978)")

    built_utc = datetime.now(timezone.utc).isoformat()
    print(f"cache.py: hashing {h5mu_path} (this is slow once)...")
    sha = _sha256_file(h5mu_path)
    header = (
        f"# provenance: dataset_id={dataset_id} "
        f"h5mu_sha256={sha} built_utc={built_utc}\n"
    )
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8") as fh:
        fh.write(header)
        frame.to_csv(fh, sep="\t", index=False)

    reports = out_tsv.parent.parent / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    def _rel(path: Path) -> str:
        try:
            return str(path.resolve().relative_to(out_tsv.parents[2].resolve()))
        except ValueError:
            return str(path)
    provenance = {
        "dataset_id": dataset_id,
        "h5mu_path": _rel(h5mu_path),
        "h5mu_sha256": sha,
        "mofa_path": _rel(mofa_path),
        "out_tsv": _rel(out_tsv),
        "n_rows": int(len(frame)),
        "columns": CACHE_COLUMNS,
        "built_utc": built_utc,
        "n_azimuth_labels": int(frame["azimuth_label"].nunique(dropna=True)),
        "n_leiden_wnn": int(frame["leiden_wnn"].nunique(dropna=True)),
    }
    prov_path = reports / "module3_label_cache_provenance.json"
    prov_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {out_tsv} ({len(frame)} rows) and {prov_path}")
    return out_tsv


def main(argv: list[str] | None = None) -> int:
    del argv
    root = Path(__file__).resolve().parents[2]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    from scripts.common.paths import load_config, module_root

    cfg: dict[str, Any] = load_config()
    m3 = cfg["module3"]
    root = module_root(cfg)
    h5mu = root / m3["multiome_dir"] / m3["files"]["secondary_h5mu"]
    mofa = root / m3["multiome_dir"] / m3["files"]["mofa_hdf5"]
    out = root / m3["files"]["label_cache"]
    build_label_cache(h5mu, out, dataset_id=str(m3["dataset_id"]), mofa_path=mofa)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
