#!/usr/bin/env python3
"""Clear stored outputs from notebooks/ (ship refresh helper)."""
from __future__ import annotations

import json
from pathlib import Path


def clear_notebook(path: Path) -> int:
    nb = json.loads(path.read_text(encoding="utf-8"))
    n = 0
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        n += len(cell.get("outputs") or [])
        cell["outputs"] = []
        cell["execution_count"] = None
    # Jupyter rewrites kernelspec.name to the registered env on save; learners
    # only have "python3". Keep display_name for the human-facing env label.
    ks = nb.setdefault("metadata", {}).setdefault("kernelspec", {})
    ks["name"] = "python3"
    if not ks.get("display_name"):
        ks["display_name"] = "CFDE_lung_env"
    if not ks.get("language"):
        ks["language"] = "python"
    path.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    return n


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    total = 0
    for path in sorted((root / "notebooks").glob("*.ipynb")):
        n = clear_notebook(path)
        total += n
        print(f"{path.relative_to(root)}: cleared {n} output payloads")
    print(f"done; {total} payloads removed")


if __name__ == "__main__":
    main()
