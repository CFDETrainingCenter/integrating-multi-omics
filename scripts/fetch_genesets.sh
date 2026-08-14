#!/usr/bin/env bash
# Fetch offline GMT libraries for Module 2 / Module 4 enrichment.
# Run from the package root (script resolves ROOT from its own location).
#
# Writes:
#   data/genesets/MSigDB_Hallmark_2020.gmt
#   data/genesets/Reactome_2022.gmt
#   data/genesets/GO_Biological_Process_2023.gmt
#   data/genesets/GENESETS_MANIFEST.json
#
# Requires: CFDE_lung_env (gseapy), network access.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

OUT="${ROOT}/data/genesets"
mkdir -p "${OUT}"

# Prefer project env if available
if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate CFDE_lung_env
fi

cd "${ROOT}"
python <<'PY'
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import gseapy as gp

out = Path("data/genesets")
out.mkdir(parents=True, exist_ok=True)

# Enrichr / gseapy library names -> local filenames (no KEGG)
libs = {
    "MSigDB_Hallmark_2020": "MSigDB_Hallmark_2020.gmt",
    "Reactome_2022": "Reactome_2022.gmt",
    "GO_Biological_Process_2023": "GO_Biological_Process_2023.gmt",
}

meta = {
    "downloaded_utc": datetime.now(timezone.utc).isoformat(),
    "source": "gseapy.get_library (Enrichr / maayanlab)",
    "organism": "Human",
    "libraries": {},
}

for lib, fname in libs.items():
    path = out / fname
    print(f"Fetching {lib} ...")
    gene_sets = gp.get_library(name=lib, organism="Human")
    lines = []
    for term, genes in gene_sets.items():
        # GMT: term, description, gene1, gene2, ...
        lines.append("\t".join([str(term), "na", *[str(g) for g in genes]]))
    path.write_text("\n".join(lines) + "\n")
    meta["libraries"][lib] = {
        "file": fname,
        "bytes": path.stat().st_size,
        "n_terms": len(gene_sets),
        "path": str(path.resolve()),
    }
    print(f"  wrote {path}  terms={len(gene_sets)}  bytes={path.stat().st_size}")

manifest = out / "GENESETS_MANIFEST.json"
manifest.write_text(json.dumps(meta, indent=2) + "\n")
print(f"Manifest: {manifest}")
print("Done.")
PY

echo
echo "Verify:"
ls -lh "${OUT}"/*.gmt "${OUT}"/GENESETS_MANIFEST.json
