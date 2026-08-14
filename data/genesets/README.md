# Gene set libraries (offline)

Required GMT files for Module 2 and Module 4 enrichment / preranked GSEA:

| File | Library |
|---|---|
| `MSigDB_Hallmark_2020.gmt` | MSigDB Hallmark 2020 |
| `Reactome_2022.gmt` | Reactome 2022 |
| `GO_Biological_Process_2023.gmt` | GO Biological Process 2023 |

**KEGG is excluded** from this course (config and loaders filter it out).

Place the three `.gmt` files in this directory (`data/genesets/`). Notebooks default to
`enrichment.mode: local_gmt` and do not call Enrichr at runtime.

Manifest: after files are present, record download date / source in `GENESETS_MANIFEST.json`
(optional; created when you fetch libraries).
