# Third-party notices

This course redistributes a small set of gene-set libraries under `data/genesets/`
(2.17 MB measured) and a derived human-mouse ortholog table under `data/reference/`.
Primary omics datasets (HuBMAP, GTEx, GEO, NASA OSDR) are **not** redistributed; learners
download them from the source portals under those portals' terms.

The notices below are required for redistribution of the gene-set libraries. See also
`data/genesets/GENESETS_MANIFEST.json`.

## Disclaimer of warranties

Required by CC BY 4.0 section 3(a)(1)(A)(iv) for the CC BY material redistributed here.

Unless otherwise separately undertaken by the Licensor, to the extent possible the Licensor offers
the Licensed Material as-is and as-available, and makes no representations or warranties of any kind
concerning the Licensed Material, whether express, implied, statutory or other. The full text of this
disclaimer, and the corresponding limitation of liability, is section 5 of the license, reproduced in
full in `LICENSE`. The same as-is, as-available terms apply to the analytical outputs of this course.

---

## MSigDB Hallmark 2020

- **Material:** `data/genesets/MSigDB_Hallmark_2020.gmt`
- **Copyright:** Broad Institute / Molecular Signatures Database (MSigDB)
- **License:** Creative Commons Attribution 4.0 International (CC BY 4.0)
- **License URL:** https://creativecommons.org/licenses/by/4.0/
- **Upstream:** https://www.gsea-msigdb.org/gsea/msigdb/

## Reactome 2022

- **Material:** `data/genesets/Reactome_2022.gmt`
- **Copyright:** Reactome
- **License:** CC0 1.0 Universal (CC0 1.0) Public Domain Dedication
- **License URL:** https://creativecommons.org/publicdomain/zero/1.0/
- **Upstream:** https://reactome.org/

## Gene Ontology Biological Process 2023 (Enrichr library export)

- **Material:** `data/genesets/GO_Biological_Process_2023.gmt`
- **Copyright:** Gene Ontology Consortium
- **License:** Creative Commons Attribution 4.0 International (CC BY 4.0)
- **License URL:** https://creativecommons.org/licenses/by/4.0/
- **Citation policy:** https://geneontology.org/docs/go-citation-policy/
- **Upstream:** https://geneontology.org/
- **Zenodo archive (concept DOI for GO releases):** https://doi.org/10.5281/zenodo.1205166

When publishing analyses that depend on GO contents, cite the GO release date and the
corresponding Zenodo DOI for that release (not only this Enrichr library label).

## Enrichr / Ma'ayan Lab (redistribution channel)

These GMT files were obtained via `gseapy.get_library` from Enrichr gene-set libraries
maintained by the Ma'ayan Lab.

- **Enrichr:** https://maayanlab.cloud/Enrichr/
- **gseapy:** https://github.com/zqfang/GSEApy

## Ensembl orthologs (derived table)

- **Material:** `data/reference/human_mouse_orthologs_one2one.tsv`
- **Source:** Ensembl BioMart one-to-one orthologs (human <-> mouse)
- **Terms:** https://www.ensembl.org/info/about/legal/index.html
- **Note:** This is a derived teaching table, not a full Ensembl redistribution.

## HuBMAP Module 3 label cache (derived table)

- **Material:** `outputs/tables/module3_label_cache_HBM828.GPVG.252.tsv`
- **Source dataset:** HuBMAP `HBM828.GPVG.252` (SNARE-seq2 Salmon + ArchR + Muon), CC BY 4.0
- **Portal:** https://portal.hubmapconsortium.org/browse/dataset/HBM828.GPVG.252
- **Note:** Derived Azimuth labels, WNN Leiden assignments, and WNN UMAP coordinates
  read from `secondary_analysis.h5mu`. The MuData file is not redistributed. Provenance
  including source SHA-256 is in `outputs/reports/module3_label_cache_provenance.json`.

## Other acknowledgments

Analyses in this course use public data from HuBMAP, GTEx, NCBI GEO (GSE150910), and
NASA OSDR / GeneLab (OSD-248). Cite those resources and their primary publications when
you reuse results. Course-held provenance is recorded in
`outputs/tables/module_data_sources.tsv` and `docs/references_verified.tsv`.
