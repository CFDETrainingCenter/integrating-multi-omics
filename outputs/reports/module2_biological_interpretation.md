# Module 2 -- Biological interpretation summary (DE + pathways)

_Auto-drafted from run at 2026-08-14T14:00:30.831325+00:00. Edit with learner conclusions -- this file does not assert biological consistency._

## Comparison design

- Grouping variable: `cell_class`
- DE method: `wilcoxon` (scanpy.tl.rank_genes_groups)
- Enrichment mode: `local_gmt` against `['MSigDB_Hallmark_2020', 'Reactome_2022', 'GO_Biological_Process_2023']`
- Ribosomal filter: `de.exclude_ribosomal_genes=True` removes RPS/RPL/MRPS/MRPL plus named exceptions UBA52 and FAU from reported marker tables, ORA inputs, and Wilcoxon prerank ranked lists (same effective enrichment flag)
- GTEx context: off in Module 2 (cross-resource comparison is Module 4)

## Top markers (per group)

- **endothelial**: HLA-E, CLDN5, EPAS1, PECAM1, EGFL7
- **epithelial**: SFTPC, SFTPA2, CYP2B7P, NAPSA, MUC1
- **immune**: EEF1A1, ZFP36L2, LAPTM5, CYBA, PTPRC
- **stromal**: GPX3, CALD1, CARMN, DST, LAMB1

## Over-representation enrichment (per group)

### Input list sizes (`n_input_genes`)

- **endothelial**: 100 genes submitted
- **epithelial**: 100 genes submitted
- **immune**: 100 genes submitted
- **stromal**: 100 genes submitted

Small lists can produce extreme or unexpected terms; treat ORA as hypothesis-generating and check gene membership before claiming biology.

- **endothelial**: Interferon Gamma Response; Endosomal/Vacuolar Pathway R-HSA-1236977; Interferon Alpha Response; Interferon Alpha/Beta Signaling R-HSA-909733; Interferon Signaling R-HSA-913531
- **epithelial**: Diseases Associated With Surfactant Metabolism R-HSA-5687613; Surfactant Metabolism R-HSA-5683826; Estrogen Response Early; Keratinization R-HSA-6805567; Diseases Of Metabolism R-HSA-5668914
- **immune**: Allograft Rejection; Immune System R-HSA-168256; Neutrophil Degranulation R-HSA-6798695; Innate Immune System R-HSA-168249; Immunoregulatory Interactions Between A Lymphoid And A non-Lymphoid Cell R-HSA-198933
- **stromal**: Epithelial Mesenchymal Transition; Extracellular Matrix Organization R-HSA-1474244; Elastic Fiber Assembly (GO:0048251); Extracellular Matrix Assembly (GO:0085029); Regulation Of Cell Migration (GO:0030334)

## Preranked GSEA (NES)

NES columns come from `gseapy.prerank` on Wilcoxon scores against local GMTs. Terms are listed for inspection; biological plausibility is the learner's task. Even after ribosomal-protein filtering, translation-adjacent Reactome/GO terms can still rank high for small secretory compartments -- treat marker genes as the identity check before claiming a pathway story.

- Reactome_2022.gmt__Immunoregulatory Interactions Between A Lymphoid And A non-Lymphoid Cell R-HSA-198933 [epithelial], NES=-2.325
- MSigDB_Hallmark_2020.gmt__Epithelial Mesenchymal Transition [immune], NES=-1.617
- GO_Biological_Process_2023.gmt__Positive Regulation Of Leukocyte Cell-Cell Adhesion (GO:1903039) [immune], NES=2.156
- GO_Biological_Process_2023.gmt__Positive Regulation Of Lymphocyte Proliferation (GO:0050671) [immune], NES=2.185
- GO_Biological_Process_2023.gmt__Positive Regulation Of T Cell Activation (GO:0050870) [immune], NES=2.192
- MSigDB_Hallmark_2020.gmt__Interferon Alpha Response [stromal], NES=-1.659
- MSigDB_Hallmark_2020.gmt__Allograft Rejection [stromal], NES=-1.660
- MSigDB_Hallmark_2020.gmt__Allograft Rejection [epithelial], NES=-1.831
- MSigDB_Hallmark_2020.gmt__Interferon Alpha Response [epithelial], NES=-1.823
- GO_Biological_Process_2023.gmt__Negative Regulation Of Protein Binding (GO:0032091) [stromal], NES=-1.932
- GO_Biological_Process_2023.gmt__B Cell Receptor Signaling Pathway (GO:0050853) [immune], NES=2.199
- Reactome_2022.gmt__Immunoregulatory Interactions Between A Lymphoid And A non-Lymphoid Cell R-HSA-198933 [immune], NES=2.234
- MSigDB_Hallmark_2020.gmt__Cholesterol Homeostasis [immune], NES=-1.617
- Reactome_2022.gmt__Immunoregulatory Interactions Between A Lymphoid And A non-Lymphoid Cell R-HSA-198933 [stromal], NES=-1.900
- GO_Biological_Process_2023.gmt__Positive Regulation Of B Cell Activation (GO:0050871) [immune], NES=2.281

## Interpretation prompts (learner task)

1. Which cell population or cluster was analyzed?
2. What are the top marker genes, and do they match known lung identities?
3. What pathways are enriched, and are they biologically coherent for that population?
4. What biological interpretation is supported by markers + pathways together?
5. What limitations remain (label transfer uncertainty, enrichment bias, small input lists)?

## Three-sentence results draft (edit me -- no auto biology claim)

Using the Module 2 integrated HuBMAP lung snRNA-seq object, we identified marker genes for `cell_class` groups with Wilcoxon rank-sum testing (scanpy) and scored pathways with gseapy (`local_gmt` ORA + prerank GSEA) against local Hallmark / Reactome / GO BP gene sets. Enriched terms are listed above with per-group `n_input_genes`; they have **not** been checked for biological plausibility in this auto-draft -- that check is the learner's task. Because labels are reference-mapped and enrichment is hypothesis-generating, follow-up validation should use independent cohorts or orthogonal assays.
