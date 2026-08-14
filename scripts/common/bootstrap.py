"""Notebook/sys.path bootstrap so `import scripts...` works."""

from __future__ import annotations

import sys
from pathlib import Path


def add_module_root_to_path() -> Path:
    """Insert CFDE_Module_Build on sys.path and return that root."""
    here = Path(__file__).resolve()
    root = here.parents[2]  # .../CFDE_Module_Build
    root_s = str(root)
    if root_s not in sys.path:
        sys.path.insert(0, root_s)
    return root
