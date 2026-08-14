#!/usr/bin/env python3
"""Orphan / stale-output detector for module runs.

Snapshot filenames and mtimes under outputs/ (and optional extra dirs) before a
run, compare after, and report:
  - expected artifacts whose mtime did not advance (stale / orphan risk)
  - new files not in the before snapshot (unexpected writes)
  - missing expected files

Usage:
  python scripts/detect_orphans.py snapshot --label pre_m4
  # ... run module ...
  python scripts/detect_orphans.py compare --label pre_m4 --expect outputs/tables/module4_*.tsv
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SNAP_DIR = ROOT / "outputs" / "reports" / "_orphans"


def _iter_files(roots: list[Path]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for base in roots:
        if not base.exists():
            continue
        for p in base.rglob("*"):
            if not p.is_file():
                continue
            rel = str(p.relative_to(ROOT))
            st = p.stat()
            out[rel] = {"mtime": st.st_mtime, "size": st.st_size}
    return out


def cmd_snapshot(args: argparse.Namespace) -> int:
    roots = [ROOT / r for r in args.roots]
    snap = {
        "label": args.label,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "roots": [str(r) for r in roots],
        "files": _iter_files(roots),
    }
    SNAP_DIR.mkdir(parents=True, exist_ok=True)
    path = SNAP_DIR / f"{args.label}.json"
    path.write_text(json.dumps(snap, indent=2) + "\n", encoding="utf-8")
    print(f"snapshot {path} n_files={len(snap['files'])}")
    return 0


def _expand_expect(patterns: list[str]) -> list[str]:
    found: list[str] = []
    for pat in patterns:
        matches = sorted(ROOT.glob(pat))
        if not matches and (ROOT / pat).exists():
            matches = [ROOT / pat]
        for m in matches:
            if m.is_file():
                found.append(str(m.relative_to(ROOT)))
    return sorted(set(found))


def cmd_compare(args: argparse.Namespace) -> int:
    before_path = SNAP_DIR / f"{args.label}.json"
    if not before_path.exists():
        raise SystemExit(f"missing snapshot: {before_path}")
    before = json.loads(before_path.read_text(encoding="utf-8"))
    roots = [ROOT / r for r in args.roots]
    after = _iter_files(roots)
    before_files = before.get("files") or {}

    expect = _expand_expect(args.expect or [])
    stale = []
    missing = []
    for rel in expect:
        if rel not in after:
            missing.append(rel)
            continue
        prev = before_files.get(rel)
        if prev is None:
            # first write this run
            continue
        if after[rel]["mtime"] <= prev["mtime"] + 1e-6:
            stale.append(rel)

    unexpected = sorted(
        rel
        for rel in after
        if rel not in before_files and (not expect or any(rel.startswith(e.split("*")[0]) for e in (args.expect or []) or ["outputs/"]))
    )
    # Narrow unexpected to outputs/ only
    unexpected = [u for u in unexpected if u.startswith("outputs/")]

    report = {
        "label": args.label,
        "n_expect": len(expect),
        "n_stale": len(stale),
        "n_missing": len(missing),
        "n_unexpected": len(unexpected),
        "stale": stale,
        "missing": missing,
        "unexpected_new": unexpected[:200],
    }
    out = SNAP_DIR / f"{args.label}_compare.json"
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in ("n_expect", "n_stale", "n_missing", "n_unexpected")}, indent=2))
    for s in stale[:30]:
        print("STALE", s)
    for m in missing[:30]:
        print("MISSING", m)
    return 1 if (stale or missing) else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("snapshot", help="Write mtime snapshot")
    sp.add_argument("--label", required=True)
    sp.add_argument(
        "--roots",
        nargs="+",
        default=["outputs/tables", "outputs/figures", "outputs/reports"],
    )
    sp.set_defaults(func=cmd_snapshot)

    cp = sub.add_parser("compare", help="Compare current tree to snapshot")
    cp.add_argument("--label", required=True)
    cp.add_argument(
        "--roots",
        nargs="+",
        default=["outputs/tables", "outputs/figures", "outputs/reports"],
    )
    cp.add_argument(
        "--expect",
        nargs="*",
        default=[],
        help="Glob(s) under package root that must advance mtime",
    )
    cp.set_defaults(func=cmd_compare)

    args = ap.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
