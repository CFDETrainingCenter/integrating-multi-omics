# Module 4 -- Cross-ecosystem comparison (methods + results draft)

_Auto-drafted at 2026-08-14T14:06:21.498263+00:00. Edit with learner conclusions._

## Framing

Module 4 is **cross-ecosystem**, not a merge. Four resources, four within-resource contrasts, plus a HuBMAP<->GTEx level bridge. Counts are **never combined** across species or studies. What travels is the direction and relative magnitude of a change (logFC / NES), never the native value.

GEO (GSE150910) is NIH but **not** CFDE. OSDR is NASA. GTEx and HuBMAP are CFDE core.

Methods assets can succeed while biology pairs stay near-null or uninterpretable: ERCC calibration and the GeneLab cross-check validate the pipeline; global NES is supporting only; shared-mechanism concordance is the headline.

## Resources and contrasts

| Ecosystem | Resource | Contrast | Native unit |
|---|---|---|---|
| NIH (not CFDE) | GSE150910 | H_DISEASE: IPF vs control | bulk counts |
| CFDE | GTEx v10 lung | H_AGING_GTEX: older vs younger (gene_reads) | bulk counts |
| NASA | OSD-248 | M_FLIGHT: ISS-T flight vs GC | bulk counts (mouse; RSEM rounded) |
| CFDE | HuBMAP | H_AGING_HUBMAP: older vs younger donors (2 vs 2 illustration) | summed counts-layer (rounded) |
| CFDE | HuBMAP / GTEx | level bridge only | linear CP10K vs TPM |

## Sample counts (from data)

- **GSE150910 / H_DISEASE / IPF**: n=103 (CHP excluded from teaching contrast)
- **GSE150910 / H_DISEASE / control**: n=103 (CHP excluded from teaching contrast)
- **OSD-248 / M_FLIGHT / ground**: n=10 (ISS-T ~60 Carcass; duration-matched GC (not hindlimb unloading))
- **OSD-248 / M_FLIGHT / flight**: n=10 (ISS-T ~60 Carcass; duration-matched GC (not hindlimb unloading))
- **OSD-248 / M_FLIGHT / flight_qa_below_threshold**: n=1 (Mmus_C57-6T_LNG_FLT_ISS-T_Rep7_F7)
- **GTEx / H_AGING_GTEX / younger_brackets**: n=100 (brackets=['20-29', '30-39']; samples=257; DE uses one sample/subject capped at 40/arm)
- **GTEx / H_AGING_GTEX / older_brackets**: n=269 (brackets=['60-69', '70-79']; samples=655; DE uses one sample/subject capped at 40/arm)
- **HuBMAP / H_AGING_HUBMAP / older**: n=2 (illustration only; donors=Donor_3,Donor_4; ages=56.8,52.96)
- **HuBMAP / H_AGING_HUBMAP / younger**: n=2 (illustration only; donors=Donor_1,Donor_2; ages=37.0,25.0)
- **HuBMAP / bridge / composition_weighted_genes**: n=42258 (linear-scale pseudobulk from Module 2 (levels only))
- **GTEx / bridge / lung_mean_tpm_genes**: n=59033 (mean TPM for levels/bridge only; age DE uses gene_reads)

## Ortholog loss (M_FLIGHT only)

- Source: Ensembl BioMart hsapiens_gene_ensembl (with_mmusculus_homolog) (2026-08-10)

- genes_in_human_contrast: 18401
- genes_in_mouse_contrast: 15429
- one2one_orthologs_available: 17041
- retained_after_intersection: 13577
- dropped_human_side: 4824
- dropped_mouse_side: 1852

## NES pairs (Hallmark + Reactome)

| Pair | Frame | n | Spearman rho | Sign agree |
|---|---|---:|---:|---:|
| H_DISEASE vs H_AGING_GTEX | biology_best_case | 1185 | -0.197 | 39.6% |
| H_DISEASE vs M_FLIGHT | cross_ecosystem | 1136 | -0.109 | 48.4% |
| H_AGING_GTEX vs M_FLIGHT | cross_ecosystem | 1158 | -0.022 | 46.3% |
| H_AGING_HUBMAP vs H_AGING_GTEX | method | 945 | 0.191 | 58.2% |

`H_DISEASE` vs `H_AGING_GTEX` is the **best-case** within-species pair (same tissue, both powered) -- not a guaranteed positive control. IPF is age-associated; it is not aging.
`H_AGING_HUBMAP` pairs only with `H_AGING_GTEX` (method question). State age-range and power differences before reading that result.

- ``spearman_p`` is not headline evidence: overlapping gene sets violate independence.
- Gene-set sensitivity table (Hallmark / Reactome / GO BP slices) ships as committed `module4_nes_sensitivity.tsv`. Default `module4.sensitivity.enabled: false` loads it; `true` re-preranks (high RAM).

## Pipeline calibration

- ERCC ExFold (flight=Mix1, GC=Mix2): n_obs=86/92; OLS slope=0.9888407236950889, intercept=-0.3398306233138027, R^2=0.8817060396863435 (Spearman secondary rho~0.9269457172513337).
- Mix is aligned with arm, so spike-ins calibrate and cannot normalize this contrast.
- `ERCC-` rows remain in the DE table for ExFold calibration; size factors are estimated on endogenous genes only (`exclude_ercc_from_size_factors`). ERCC read fraction is recorded by sample; excluding spike-ins from size factors does not move pathway passer `p_up` off zero.

## Confounds (one clause each)

- GTEx is postmortem: age is entangled with cause of death and ischemic time.
- OSD-248 is female C57BL/6NTac at 36 weeks: nothing here estimates an age effect.
- GSE150910 controls are not matched to GTEx donors.
- ISS-T flight carcasses were transported from orbit; ground carcasses were not (both held at -80; dissection months after euthanasia). Flight euthanasia dates 06-09 Feb vs GC 03-06 Feb; dissection dates identical; **no times of day recorded** (circadian genes often dominate non-ERCC ranks).
- H_DISEASE uses `~ condition` only; the published analysis adjusted for clinical covariates.
- OSD-248 ERCC spike-ins often dominate padj ranks because flight and GC used different ExFold mixes; volcanoes exclude ERCC by default.
- H_AGING_HUBMAP is 2 vs 2 illustration only; HuBMAP older ages sit inside GTEx middle brackets.
- HuBMAP integrated `counts` layer is not integer UMIs; donor sums are rounded before pydeseq2.
- No CFDE program supplies fibrotic human lung; that is why Module 4 goes to GEO (NIH, not CFDE).
- Reference levels are set explicitly (`Treatment(reference)` + DeseqStats contrast). Alphabetical factor order would have inverted three of four contrasts (`IPF` before `control`; `older` before `younger`).

## What this module teaches

1. Levels never travel (four resources, four unit systems).
2. Global NES correlations show these contrasts are not broadly equivalent (supporting line).
3. Shared-mechanism concordance names pathways that pass thresholds in the same direction, against chance; that is the headline.
4. Validate method before interpreting biology (ERCC, GeneLab; explicit DE polarity).
5. An underpowered contrast looks like nothing; recognizing that is a skill.

## Shared-mechanism concordance (headline)

Global NES Spearman is supporting only. The headline is pathways that pass FDR q<0.25 and |NES|>=1.5 in each contrast independently, with agreeing sign, compared to analytic and permutation nulls. All pathways ship as tier 3 (concordance only). Leading-edge Jaccard on this run spans roughly 0.21 to 0.48 with no natural tier 1/2 cut; do not invent one.

Suspect ribosome / translation / OxPhos terms are flagged, not excluded.

| Set | Obs | Unflagged | Perm mean | Obs/perm | perm p |
|---|---:|---:|---:|---:|---:|
| disease_aging | 5 | 5 | 3.01 | 1.66 | 0.169 |
| disease_flight | 0 | 0 | 0.42 | 0.00 | 1.000 |
| aging_flight | 6 | 3 | 2.48 | 2.42 | 0.025 |
| three_way | 0 | 0 | 0.03 | 0.00 | 1.000 |

- **disease_aging**: 5 mechanisms unflagged (5 including flagged); 5 concordant pathways (perm p=0.169).
- **disease_flight**: 0 mechanisms unflagged (0 including flagged); 0 concordant pathways (perm p=1.000).
- **aging_flight**: 3 mechanisms unflagged (5 including flagged); 6 concordant pathways (perm p=0.025).
- **three_way**: 0 mechanisms unflagged (0 including flagged); 0 concordant pathways (perm p=1.000).

Flagged representatives: Respiratory Electron Transport, ATP Synthesis By Chemiosmotic Coupling, Heat Production By Uncoupling Proteins R-HSA-163200; Oxidative Phosphorylation.

