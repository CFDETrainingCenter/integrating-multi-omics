"""Fetch Azimuth / obs labels from a public HuBMAP secondary_analysis.h5ad via HTTP range.

Ported from _refinement/reference/remote_obs_probe.py. Downloads only the HDF5
pages h5py touches (no fsspec). Used by the QC / label-join path when
module2.labels.fetch_remote is enabled, or to rebuild the cached TSV.
"""

from __future__ import annotations

import io
import urllib.request
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd

UA = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/127.0 Safari/537.36"
)

HUBMAP_ASSET_URL = "https://assets.hubmapconsortium.org/{uuid}/secondary_analysis.h5ad"


class HTTPRangeFile(io.RawIOBase):
    """Minimal seekable read-only file over HTTP byte ranges."""

    def __init__(self, url: str, block: int = 1 << 20):
        self.url = url
        self.block = block
        self.pos = 0
        self.bytes_fetched = 0
        self.requests = 0
        self.size = self._head()

    def _head(self) -> int:
        req = urllib.request.Request(
            self.url,
            method="GET",
            headers={"User-Agent": UA, "Range": "bytes=0-0"},
        )
        with urllib.request.urlopen(req, timeout=60) as r:
            cr = r.headers.get("Content-Range", "")
        return int(cr.split("/")[-1])

    def _fetch(self, start: int, end: int) -> bytes:
        end = min(end, self.size - 1)
        req = urllib.request.Request(
            self.url,
            headers={"User-Agent": UA, "Range": f"bytes={start}-{end}"},
        )
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read()
        self.requests += 1
        self.bytes_fetched += len(data)
        return data

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def seek(self, offset: int, whence: int = 0) -> int:
        self.pos = (
            offset
            if whence == 0
            else self.pos + offset
            if whence == 1
            else self.size + offset
        )
        return self.pos

    def tell(self) -> int:
        return self.pos

    def read(self, n: int = -1) -> bytes:
        if n < 0:
            n = self.size - self.pos
        if n == 0 or self.pos >= self.size:
            return b""
        data = self._fetch(self.pos, self.pos + n - 1)
        self.pos += len(data)
        return data

    def readinto(self, b) -> int:
        d = self.read(len(b))
        b[: len(d)] = d
        return len(d)


def categories_for(h5: h5py.File, col: str) -> tuple[list[str] | None, np.ndarray]:
    """Resolve a categorical column's category strings (old and new AnnData layouts).

    Donor_1 uses obs/__categories/; newer donors use Group with categories/codes.
    """
    obs = h5["obs"]
    node = obs[col]
    if isinstance(node, h5py.Group):
        raw = node["categories"][:]
        codes = node["codes"][:]
    else:
        codes = node[:]
        raw = None
        for cand in (f"{col}_categories", col):
            if "uns" in h5 and cand in h5["uns"]:
                raw = h5["uns"][cand][:]
                break
        if raw is None and "__categories" in obs:
            raw = obs["__categories"][col][:]
    if raw is None:
        return None, codes
    cats = [c.decode() if isinstance(c, bytes) else str(c) for c in raw]
    return cats, codes


def secondary_analysis_url(uuid: str) -> str:
    return HUBMAP_ASSET_URL.format(uuid=str(uuid).strip())


def fetch_azimuth_labels(
    url: str,
    *,
    column: str = "azimuth_label",
) -> tuple[pd.Series, dict[str, Any]]:
    """Return barcode-indexed azimuth labels from a remote secondary_analysis.h5ad."""
    f = HTTPRangeFile(url)
    with h5py.File(f, "r") as h:
        cats, codes = categories_for(h, column)
        if cats is None:
            raise KeyError(f"categories for obs/{column} not found in {url}")
        bc_attr = h["obs"].attrs.get("_index", b"_index")
        bc_key = bc_attr.decode() if isinstance(bc_attr, bytes) else str(bc_attr)
        barcodes = h["obs"][bc_key][:]
        barcodes = [
            b.decode() if isinstance(b, bytes) else str(b) for b in barcodes
        ]
        labels = [cats[int(c)] if int(c) >= 0 else "" for c in codes]
    series = pd.Series(labels, index=pd.Index(barcodes, name="barcode"), name=column)
    meta = {
        "url": url,
        "column": column,
        "n_nuclei": int(len(series)),
        "n_categories": int(len(cats)),
        "bytes_fetched": int(f.bytes_fetched),
        "requests": int(f.requests),
        "remote_size_bytes": int(f.size),
    }
    return series, meta


def fetch_azimuth_labels_for_uuid(
    uuid: str,
    *,
    column: str = "azimuth_label",
) -> tuple[pd.Series, dict[str, Any]]:
    return fetch_azimuth_labels(secondary_analysis_url(uuid), column=column)


def fetch_cohort_azimuth_labels(cfg: dict[str, Any]) -> pd.DataFrame:
    """Fetch azimuth_label for each module2 input that has a uuid; long TSV shape."""
    rows: list[dict[str, Any]] = []
    metas: list[dict[str, Any]] = []
    for item in cfg.get("module2", {}).get("inputs") or []:
        uuid = item.get("uuid")
        if not uuid:
            continue
        series, meta = fetch_azimuth_labels_for_uuid(str(uuid))
        metas.append(
            {
                "donor_label": item.get("donor_label"),
                "processed_id": item.get("processed_id") or item.get("block_id"),
                "primary_id": item.get("primary_id"),
                **meta,
            }
        )
        for barcode, label in series.items():
            rows.append(
                {
                    "donor_label": item.get("donor_label"),
                    "processed_id": item.get("processed_id") or item.get("block_id"),
                    "primary_id": item.get("primary_id"),
                    "uuid": uuid,
                    "barcode": barcode,
                    "azimuth_label": label,
                }
            )
    out = pd.DataFrame(rows)
    out.attrs["fetch_meta"] = metas
    return out


def write_azimuth_label_cache(cfg: dict[str, Any], path: Path | str | None = None) -> Path:
    """Fetch remote labels and write module2.labels.cache TSV."""
    labels_cfg = (cfg.get("module2") or {}).get("labels") or {}
    rel = path or labels_cfg.get("cache") or "outputs/tables/module2_azimuth_labels.tsv"
    root = Path(cfg["_module_root"])
    out = Path(rel) if Path(rel).is_absolute() else root / rel
    out.parent.mkdir(parents=True, exist_ok=True)
    df = fetch_cohort_azimuth_labels(cfg)
    df.to_csv(out, sep="\t", index=False)
    return out
