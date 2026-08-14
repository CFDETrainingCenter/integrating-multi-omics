"""I/O helpers for HuBMAP metadata and AnnData objects."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


def read_hubmap_metadata_tsv(path: Path | str) -> pd.DataFrame:
    """Read HuBMAP Key/Value metadata TSV."""
    df = pd.read_csv(path, sep="\t")
    expected = {"HuBMAP ID", "Entity", "Key", "Value"}
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")
    return df


def metadata_kv(df: pd.DataFrame, entity: str | None = None) -> dict[str, str]:
    """Collapse Key/Value rows to a dict; optionally filter by Entity."""
    sub = df if entity is None else df[df["Entity"] == entity]
    out: dict[str, str] = {}
    for _, row in sub.iterrows():
        key = str(row["Key"])
        val = "" if pd.isna(row["Value"]) else str(row["Value"])
        out[key] = val
    return out


def donor_and_dataset_ids(df: pd.DataFrame) -> tuple[list[str], list[str]]:
    donor_ids = (
        df.loc[df["Entity"] == "Donor", "HuBMAP ID"].dropna().astype(str).unique().tolist()
    )
    dataset_ids = (
        df.loc[df["Entity"] != "Donor", "HuBMAP ID"].dropna().astype(str).unique().tolist()
    )
    return donor_ids, dataset_ids


def read_h5ad(path: Path | str, backed: str | None = None):
    """Lazy import scanpy/anndata to keep discovery notebooks light."""
    import scanpy as sc

    return sc.read_h5ad(path, backed=backed)


def write_h5ad(adata, path: Path | str, **kwargs: Any) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path, **kwargs)
    return path


def environment_versions() -> dict[str, str]:
    """Record versions for packages that affect analytical results."""
    from importlib.metadata import PackageNotFoundError, version as pkg_version

    versions: dict[str, str] = {}
    # display_name -> distribution name on PyPI / conda
    packages = {
        "pandas": "pandas",
        "numpy": "numpy",
        "scipy": "scipy",
        "matplotlib": "matplotlib",
        "anndata": "anndata",
        "scanpy": "scanpy",
        "yaml": "PyYAML",
        "harmonypy": "harmonypy",
        "leidenalg": "leidenalg",
        "igraph": "igraph",
        "umap": "umap-learn",
        "sklearn": "scikit-learn",
        "gseapy": "gseapy",
        "mudata": "mudata",
        "muon": "muon",
        "pydeseq2": "pydeseq2",
        "h5py": "h5py",
    }
    import_fallbacks = {
        "yaml": "yaml",
        "umap": "umap",
        "sklearn": "sklearn",
    }
    for display, dist in packages.items():
        try:
            versions[display] = pkg_version(dist)
            continue
        except PackageNotFoundError:
            pass
        except Exception as exc:  # noqa: BLE001
            versions[display] = f"unavailable ({exc.__class__.__name__})"
            continue
        mod_name = import_fallbacks.get(display, display)
        try:
            mod = __import__(mod_name)
            versions[display] = getattr(mod, "__version__", "unknown")
        except Exception as exc:  # noqa: BLE001
            versions[display] = f"unavailable ({exc.__class__.__name__})"
    return versions
