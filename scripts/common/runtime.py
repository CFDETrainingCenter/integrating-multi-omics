"""Runtime and input-integrity helpers for reproducible runs."""

from __future__ import annotations

import hashlib
import os
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator


def sha256_file(path: Path | str, *, chunk: int = 1 << 20) -> str:
    """Hex sha256 of a file on disk."""
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def sha256_optional(path: Path | str | None) -> str | None:
    if path is None:
        return None
    p = Path(path)
    if not p.exists() or not p.is_file():
        return None
    return sha256_file(p)


def peak_rss_mb() -> float | None:
    """Best-effort peak resident set size in MiB (macOS/Linux)."""
    try:
        import resource

        usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # macOS reports bytes; Linux reports kilobytes
        if sys_platform_is_darwin():
            return round(float(usage) / (1024.0 * 1024.0), 2)
        return round(float(usage) / 1024.0, 2)
    except Exception:  # noqa: BLE001
        return None


def sys_platform_is_darwin() -> bool:
    return os.uname().sysname == "Darwin"


@contextmanager
def measure_run() -> Iterator[dict[str, Any]]:
    """Wall-clock and process peak RSS for a scoped block.

    ``peak_rss_mb`` is ``ru_maxrss`` for the process so far (monotone across a
    multi-donor loop), not a per-donor peak. Call ``finalize_timing(timing)``
    before any code that snapshots the dict (e.g. writing run_params).
    """
    out: dict[str, Any] = {
        "_compute_seconds": None,
        "peak_rss_mb": None,
        "peak_rss_mb_start": peak_rss_mb(),
        "peak_rss_note": (
            "process peak RSS so far (ru_maxrss); monotone across a loop, not per-step"
        ),
        "_t0": time.perf_counter(),
    }
    try:
        yield out
    finally:
        finalize_timing(out)


def finalize_timing(timing: dict[str, Any]) -> dict[str, Any]:
    """Stamp wall-clock and peak RSS onto a measure_run dict (safe to call twice)."""
    t0 = timing.get("_t0")
    if t0 is not None:
        timing["_compute_seconds"] = round(time.perf_counter() - float(t0), 3)
    timing["peak_rss_mb"] = peak_rss_mb()
    timing.setdefault(
        "peak_rss_note",
        "process peak RSS so far (ru_maxrss); monotone across a loop, not per-step",
    )
    return timing


def feature_space_counts(var_names) -> dict[str, int]:
    """Count exonic vs Salmon intron (`*-I`) features in a gene index."""
    names = [str(x) for x in var_names]
    n_intron = sum(1 for n in names if n.endswith("-I"))
    n_exonic = len(names) - n_intron
    return {
        "n_features_total": len(names),
        "n_features_exonic": n_exonic,
        "n_features_intron": n_intron,
    }


def matrix_memory_lines(adata) -> list[str]:
    """Describe dense vs sparse footprint of `.X` and `layers['counts']`."""
    try:
        import scipy.sparse as sp
    except Exception:  # noqa: BLE001
        return []

    lines: list[str] = []
    for name, M in (("X", getattr(adata, "X", None)), ("counts", (getattr(adata, "layers", {}) or {}).get("counts"))):
        if M is None:
            continue
        dense = not sp.issparse(M)
        if dense:
            nbytes = int(getattr(M, "nbytes", 0) or 0)
        else:
            nbytes = int(M.data.nbytes + M.indices.nbytes + M.indptr.nbytes)
        lines.append(
            f"{name}: {'DENSE' if dense else 'sparse'} {getattr(M, 'dtype', '?')} {nbytes / 1e9:.2f} GB"
        )
    return lines


def rss_checkpoint(label: str, adata=None, *, print_fn=print) -> dict[str, Any]:
    """Print peak RSS and optional matrix sparsity/density at a named step."""
    info: dict[str, Any] = {
        "label": label,
        "peak_rss_mb": peak_rss_mb(),
    }
    print_fn(f"[RSS {label}] peak_rss_mb={info['peak_rss_mb']}")
    if adata is not None:
        for line in matrix_memory_lines(adata):
            print_fn(f"[RSS {label}] {line}")
        info["matrix_memory"] = matrix_memory_lines(adata)
    return info
