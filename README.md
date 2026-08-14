# CFDE lung transcriptomics training

A five-module course in lung transcriptomics built on Common Fund Data Ecosystem resources. You
discover datasets through CFDE surfaces, download them from the DCC portals, and carry one analysis
from raw nuclei through quality control, cross-donor integration, differential expression and
pathway enrichment, paired multiome interpretation, and a cross-ecosystem contrast comparison.

Primary omics files are **not** included and are not redistributed. The folders they belong in are
already named and empty; you download into them. Module 0 sets up the environment and Module 1
covers discovery, download and verification, with terminal commands for every file in its
Appendix A.

## Attribution

Copyright (c) 2026 Brian Billings.

**Title.** CFDE lung transcriptomics course (HuBMAP / GTEx bridge with GEO and OSDR arms).

**Licenses.** Course prose, documents, notebooks, and derived tables and figures:
[CC BY 4.0](LICENSE). Code under `scripts/`: [MIT](LICENSE-CODE). Third-party gene-set libraries
under `data/genesets/` remain under their upstream terms; see
[THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md).

## What is in this package

| Included | Not included |
|---|---|
| Six module documents plus the instructor answer keys, under `docs/` | HuBMAP / GTEx / GEO / OSDR primary files |
| Notebooks 01 to 04, outputs **cleared** | Alternate HuBMAP blocks and SNARE products other than Donor_2 `HBM828.GPVG.252` |
| Analysis package under `scripts/` and `config/paths.yaml` | Regenerated `data/processed/*.h5ad` |
| Empty named folders for every teaching HuBMAP ID | Full optional MuData (~14.7 GiB) |
| Committed figures and reports under `outputs/` | Tables the pipeline regenerates when it runs |
| Gene sets, ortholog and ERCC reference files | |

`outputs/tables/` holds only the files the pipeline reads rather than writes: the two Azimuth label
caches, the committed GeneLab differential-expression subset, the NES sensitivity table, the offline
discovery caches, and the download manifest. Everything else appears when you run the notebooks.

## Documents

| Module | Document | Time | Notebook |
|---|---|---|---|
| 0 | `docs/module0.docx` | 20 min reading | none, environment setup and code orientation |
| 1 | `docs/module1.docx` | 45 min | `notebooks/01_discovery_and_qc.ipynb` |
| 2 | `docs/module2.docx` | 45 min | `notebooks/02_integration_de_pathways.ipynb` |
| 3 | `docs/module3.docx` | 30 min | `notebooks/03_multiome.ipynb` |
| 4 | `docs/module4.docx` | 30 min | `notebooks/04_cross_ecosystem.ipynb` |
| 5 | `docs/module5.docx` | 30 min | none, communication and assessment |

Three hours of instruction across Modules 1 to 5, excluding Module 0.

`docs/answer_keys.docx` is an instructor document covering the Module 5 knowledge check and the
capstone rubric. It is not learner material.

## What you need

| Requirement | Value |
|---|---|
| Python | 3.11 |
| Environment | conda/mamba `CFDE_lung_env`, or a pip venv from `requirements.txt` |
| RAM | **16 GB recommended**; 8 GB minimum (measured peaks: Module 2 DE 7.6 GB, Module 3 7.1 GB, Module 4 5.5 GB) |
| Disk | 5 GB, covering the download plus the objects the notebooks write |
| Download | About 815 MB across 20 files, itemized in Module 0 Appendix B with retrieval commands in Module 1 Appendix A |
| Gene sets | Shipped under `data/genesets/` (Hallmark, Reactome, GO BP; no KEGG) |

## Quick start

```bash
mamba env create -f environment.yml
mamba activate CFDE_lung_env
python -m ipykernel install --user --name cfde_lung_env --display-name "Python (CFDE_lung_env)"
```

Read `docs/module0.docx` first: it covers the environment, the code layout, and the ten Python
patterns the analysis packages use. Then work Module 1, which walks discovery, download and
verification before any analysis runs.

```bash
python scripts/verify_access.py          # confirm every taught URL is reachable
jupyter lab notebooks/01_discovery_and_qc.ipynb
```

Pip alternative:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Cohort

**Donor_1, Donor_2, Donor_3, Donor_4**, four HuBMAP lung donors selected for comparative value
across age, sex and reported race and ethnicity.

| Donor | Age | Sex | snRNA block | Also provides |
|---|---|---|---|---|
| Donor_1 | 37 | Male | `HBM347.VDHS.379` | Module 1 QC scaffold |
| Donor_2 | 25 | Female | `HBM473.NKMR.872` | SNARE-seq2 multiome `HBM828.GPVG.252` (Module 3) |
| Donor_3 | 56.8 | Female | `HBM484.RGWP.797` | |
| Donor_4 | 52.96 | Male | `HBM379.QCJT.435` | |

Module 4 adds GTEx v10 lung, GEO series GSE150910, and NASA OSDR study OSD-248.

## Interpreting without running

Notebooks ship without stored cell outputs on purpose, so a clean re-run is forced after download.
The results are readable before you run anything:

1. `docs/module{0-5}.docx`, the module documents, with figures and tables embedded
2. `outputs/figures/`, every figure the documents cite
3. `outputs/reports/`, methods records and biological interpretation drafts

## Design notes

- No primary data redistribution: DCC portals for files, CFDE surfaces for discovery.
- GTEx enters at **Module 4**, not Module 2.
- Gene sets: Hallmark, Reactome and GO BP only (**no KEGG**).
- Module 4's headline is mechanism concordance against a permutation null; the global NES
  correlation is supporting context only.
- Every number in the documents is read from a committed table rather than transcribed.

## Optional smoke run

After the download is complete:

```bash
python scripts/run_modules_smoke.py --module all
```

## Folders you populate

`pdfs/` and `screencasts/` are placeholders for the rendered PDF documents and the module
screencasts.
