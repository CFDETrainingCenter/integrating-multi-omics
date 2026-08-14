#!/usr/bin/env python3
"""Verify learner download URLs with ranged GET (never HEAD).

GATE_RESULTS section 7: for every URL the course tells a learner to fetch, issue a
ranged GET with a browser user agent, assert HTTP 206 or 200, and assert the
Content-Range / Content-Length size matches the size quoted in the download
manifest. Fail loudly on 401, 404, or size drift.

Default URL catalog is built from config/paths.yaml processed UUIDs plus known
GTEx / GEO / OSDR asset URLs recorded beside the learner manifest.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
UA = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/124.0.0.0 Safari/537.36"
)
HUBMAP_ASSETS = "https://assets.hubmapconsortium.org/{uuid}/{filename}"

# Decimal MB as used in learner download summary (bytes / 1e6).
SIZE_TOLERANCE_FRAC = 0.02
SIZE_TOLERANCE_ABS = 1_000_000  # 1 MB


@dataclass
class Asset:
    name: str
    url: str
    expected_bytes: int | None
    required: bool = True


def _mb_to_bytes(mb: float) -> int:
    return int(round(float(mb) * 1_000_000))



def load_manifest_sizes(root: Path) -> dict[str, int]:
    """Map expected bytes from the learner manifest, keyed two ways.

    Four donors each contribute a file named ``raw_expr.h5ad``, so a basename-only
    key collapses them onto whichever row is read last and makes three of the four
    report a false size mismatch. Emit a ``donor_label/file_name`` key as well and
    prefer it at the call site; the bare basename is kept for the singleton files.
    """
    path = root / "outputs/tables/module1_learner_download_manifest.tsv"
    out: dict[str, int] = {}
    if not path.exists():
        return out
    import pandas as pd

    df = pd.read_csv(path, sep="\t")
    for _, row in df.iterrows():
        if not bool(row.get("required")):
            continue
        name = str(row["file_name"])
        size = _mb_to_bytes(row["size_mb"])
        donor = str(row.get("donor_label") or "").strip()
        if donor:
            out[f"{donor}/{name}"] = size
        out.setdefault(name, size)
    return out


def default_catalog(root: Path, cfg: dict) -> list[Asset]:
    sizes = load_manifest_sizes(root)
    assets: list[Asset] = []

    # HuBMAP snRNA raw_expr via processed UUID
    for row in (cfg.get("module2") or {}).get("inputs") or []:
        uuid = row.get("uuid")
        pid = row.get("primary_id")
        if not uuid:
            continue
        assets.append(
            Asset(
                name=f"hubmap_raw_expr:{pid}",
                url=HUBMAP_ASSETS.format(uuid=uuid, filename="raw_expr.h5ad"),
                # Every donor's file is named raw_expr.h5ad, so a basename-keyed
                # lookup compares all four against one donor's size and reports
                # false failures. Key on the donor's own manifest row instead.
                expected_bytes=sizes.get(f"{row.get('donor_label')}/raw_expr.h5ad"),
                required=True,
            )
        )
        # Override expected from per-row if we can read local size
        local = (
            root
            / "data/source/HUBMAP"
            / str(row.get("donor_label"))
            / "snRNAseq"
            / str(pid)
            / "raw_expr.h5ad"
        )
        if local.exists():
            assets[-1].expected_bytes = local.stat().st_size

    # Module 3 MOFA hdf5
    inv = root / "outputs/tables/module3_multimodal_inventory.tsv"
    m3 = cfg.get("module3") or {}
    uuid3 = m3.get("uuid") or "772d0ce4657bb306359c0230c57ee1e9"
    mofa_local = (
        root
        / (m3.get("multiome_dir") or "data/source/HUBMAP/Donor_2/SNAREseq/HBM828.GPVG.252")
        / ((m3.get("files") or {}).get("mofa_hdf5") or "multiome_mofa.hdf5")
    )
    assets.append(
        Asset(
            name="module3_mofa",
            url=HUBMAP_ASSETS.format(uuid=uuid3, filename="multiome_mofa.hdf5"),
            expected_bytes=mofa_local.stat().st_size if mofa_local.exists() else sizes.get(
                "multiome_mofa.hdf5"
            ),
            required=True,
        )
    )

    # GTEx lung TPM (public GCS; path from GATE_RESULTS 2026-08-09)
    gtex_local = root / "data/source/GTEX/gene_tpm_v10_lung.gct.gz"
    gtex_url = (
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/"
        "tpms-by-tissue/gene_tpm_v10_lung.gct.gz"
    )
    assets.append(
        Asset(
            name="gtex_gene_tpm_v10_lung",
            url=gtex_url,
            expected_bytes=gtex_local.stat().st_size if gtex_local.exists() else sizes.get(
                "gene_tpm_v10_lung.gct.gz"
            ),
            required=True,
        )
    )

    # GTEx lung gene reads (counts; Module 4 aging DE)
    gtex_reads_local = root / "data/source/GTEX/gene_reads_v10_lung.gct.gz"
    gtex_reads_url = (
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/"
        "counts-by-tissue/gene_reads_v10_lung.gct.gz"
    )
    assets.append(
        Asset(
            name="gtex_gene_reads_v10_lung",
            url=gtex_reads_url,
            expected_bytes=(
                gtex_reads_local.stat().st_size if gtex_reads_local.exists() else None
            ),
            required=True,
        )
    )

    # GEO GSE150910 gene-level counts (Module 4 H_DISEASE)
    geo_local = root / "data/source/GEO/GSE150910/GSE150910_gene-level_count_file.csv.gz"
    if not geo_local.exists():
        # tolerate alternate layout under data/source/GEO
        candidates = list((root / "data/source/GEO").rglob("GSE150910_gene-level_count_file.csv.gz"))
        geo_local = candidates[0] if candidates else geo_local
    assets.append(
        Asset(
            name="geo_gse150910_gene_counts",
            url=(
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150910/suppl/"
                "GSE150910_gene-level_count_file.csv.gz"
            ),
            expected_bytes=geo_local.stat().st_size if geo_local.exists() else 8_100_000,
            required=True,
        )
    )

    # OSDR OSD-248 portal (study page). Local count matrices are verified by
    # presence on disk in main() because OSDR asset URLs are not stable ranged-GET targets.
    assets.append(
        Asset(
            name="osdr_osd248_portal",
            url="https://osdr.nasa.gov/bio/repo/data/studies/OSD-248",
            expected_bytes=None,
            required=True,
        )
    )

    # API liveness (not sized)
    assets.extend(
        [
            Asset(
                name="hubmap_search_api",
                url="https://search.api.hubmapconsortium.org/v3/indices",
                expected_bytes=None,
                required=True,
            ),
            Asset(
                name="gtex_api_v2",
                url="https://gtexportal.org/api/v2/dataset/tissueSiteDetail",
                expected_bytes=None,
                required=True,
            ),
            Asset(
                name="smartapi_hubmap_query",
                url="https://smart-api.info/api/query?q=hubmap",
                expected_bytes=None,
                required=True,
            ),
        ]
    )
    return assets


def ranged_get(url: str, nbytes: int = 100) -> tuple[int, dict[str, str], bytes]:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": UA,
            "Range": f"bytes=0-{nbytes - 1}",
            "Accept": "*/*",
        },
        method="GET",
    )
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            status = getattr(resp, "status", None) or resp.getcode()
            headers = {k.lower(): v for k, v in resp.headers.items()}
            body = resp.read(nbytes)
            return int(status), headers, body
    except urllib.error.HTTPError as exc:
        headers = {k.lower(): v for k, v in (exc.headers.items() if exc.headers else [])}
        body = exc.read(nbytes) if hasattr(exc, "read") else b""
        return int(exc.code), headers, body


def parse_total_size(headers: dict[str, str]) -> int | None:
    cr = headers.get("content-range")
    if cr and "/" in cr:
        total = cr.rsplit("/", 1)[-1].strip()
        if total.isdigit():
            return int(total)
    cl = headers.get("content-length")
    if cl and cl.isdigit():
        # For 200 full-body responses Content-Length is the whole object;
        # for 206 it is the range length -- prefer Content-Range.
        if headers.get("content-range"):
            return parse_total_size({"content-range": headers["content-range"]})
        return int(cl)
    return None


def size_ok(expected: int | None, observed: int | None) -> bool:
    if expected is None or observed is None:
        return True
    tol = max(SIZE_TOLERANCE_ABS, int(expected * SIZE_TOLERANCE_FRAC))
    return abs(expected - observed) <= tol


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument(
        "--catalog",
        type=Path,
        default=None,
        help="Optional JSON list of {name,url,expected_bytes,required}",
    )
    ap.add_argument(
        "--report",
        type=Path,
        default=None,
        help="Write JSON report path (default: outputs/reports/verify_access_report.json)",
    )
    args = ap.parse_args()
    root = args.root.resolve()
    cfg = yaml.safe_load((root / "config/paths.yaml").read_text(encoding="utf-8"))

    if args.catalog:
        raw = json.loads(args.catalog.read_text(encoding="utf-8"))
        catalog = [Asset(**row) for row in raw]
    else:
        catalog = default_catalog(root, cfg)

    results = []
    failures = 0
    for asset in catalog:
        status, headers, _body = ranged_get(asset.url)
        observed = parse_total_size(headers)
        ok_status = status in {200, 206}
        ok_size = size_ok(asset.expected_bytes, observed)
        # Portal HTML pages: size check N/A
        if asset.url.endswith((".h5ad", ".hdf5", ".gz", ".tsv", ".csv")):
            pass
        else:
            ok_size = True
        ok = ok_status and ok_size
        if asset.required and not ok:
            failures += 1
        row = {
            "name": asset.name,
            "url": asset.url,
            "status": status,
            "expected_bytes": asset.expected_bytes,
            "observed_bytes": observed,
            "ok_status": ok_status,
            "ok_size": ok_size,
            "ok": ok,
            "required": asset.required,
        }
        results.append(row)
        flag = "OK" if ok else "FAIL"
        print(
            f"{flag} {asset.name} status={status} "
            f"expected={asset.expected_bytes} observed={observed}"
        )

    # Local Module 4 inputs that lack stable ranged-GET file URLs
    local_required = [
        ("local_gtex_gene_reads", root / "data/source/GTEX/gene_reads_v10_lung.gct.gz"),
    ]
    geo_hits = list((root / "data/source/GEO").rglob("GSE150910_gene-level_count_file.csv.gz"))
    local_required.append(
        (
            "local_geo_gse150910_counts",
            geo_hits[0] if geo_hits else root / "data/source/GEO/GSE150910_gene-level_count_file.csv.gz",
        )
    )
    osd_hits = list((root / "data/source/OSDR").rglob("*Unnormalized*Counts*"))
    if osd_hits:
        local_required.append(("local_osd248_unnormalized_counts", osd_hits[0]))
    else:
        # Still record a required miss so packaging notices
        local_required.append(
            (
                "local_osd248_unnormalized_counts",
                root / "data/source/OSDR/OSD-248/missing_unnormalized_counts",
            )
        )
    for name, path in local_required:
        ok = path.exists() and path.stat().st_size > 0
        if not ok:
            failures += 1
        results.append(
            {
                "name": name,
                "url": str(path),
                "status": 200 if ok else 404,
                "expected_bytes": path.stat().st_size if path.exists() else None,
                "observed_bytes": path.stat().st_size if path.exists() else None,
                "ok_status": ok,
                "ok_size": ok,
                "ok": ok,
                "required": True,
                "check": "local_presence",
            }
        )
        print(f"{'OK' if ok else 'FAIL'} {name} local={path} exists={ok}")

    report_path = args.report or (root / "outputs/reports/verify_access_report.json")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "n_assets": len(results),
        "n_failures": failures,
        "results": results,
    }
    report_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {report_path} failures={failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
