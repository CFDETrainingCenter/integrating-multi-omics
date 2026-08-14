# Module 3 -- Donor_2 SNARE multiome interpretation

_Auto-drafted at 2026-08-14T14:00:55.288626+00:00. Edit with learner conclusions._

## Role in the course

This component reads two modalities measured in the same nucleus (paired RNA + ATAC). Modules 1-2 and 4 are multi-*resource* transcriptomics.

## Design

- Donor: `HBM994.NXZC.854` (Donor_2)
- Dataset: `HBM828.GPVG.252`
- Assay: `SNARE-seq2 [Salmon + ArchR + Muon]`
- Portal: https://portal.hubmapconsortium.org/browse/dataset/HBM828.GPVG.252
- Default path: interpret fitted MOFA (`multiome_mofa.hdf5`) plus a committed Azimuth/WNN label cache; full MuData is not a learner download (`load_h5mu`)

## Barcode pairing

- RNA barcodes: 21978
- ATAC barcodes: 21978
- Exact overlap: 21978
- Paired multiome: True
- Note: NOT independent pairing evidence from a barcode crosstab. Pairing is asserted by MOFA file construction: one shared sample barcode list is written once for rna + atac_cbg views. n_exact_overlap equals that list length by definition (MuData barcode crosstab skipped because load_h5mu=false).

## MOFA: modality-private factors (not a joint latent space)

- View `rna`: MOFA explains ~1.31 percent of total variance (stored on a percent scale in the HDF5, so 1.3128 is 1.31 percent and not a fraction above 1).
- View `atac_cbg`: MOFA explains ~0.91 percent of total variance (stored on a percent scale in the HDF5, so 0.9145 is 0.91 percent and not a fraction above 1).

### R^2 recomputation (from Z, W, Y in the HDF5)

- Recomputation matches stored/100 (percent scale); Gate 2 confirms convention before swapping published %.
- Details: `{"ok": true, "view": "rna", "rna_r2_recomputed": 0.013128170817511609, "rna_r2_stored": 1.3128374404452026, "ss_res": 74844399.34997773, "ss_tot": 75840040.35455936, "n_nan_in_Y": 0, "matches_stored": true, "matches_direct": false, "matches_percent": true, "scale_branch": "percent", "W_shape": [30, 8991], "Z_shape": [21978, 30], "Y_shape": [21978, 8991]}`

**Finding.** No factor explains more than 0.0004 percent of variance in its non-dominant view (HDF5 value 0.00042; units match stored totals = percent).
This MOFA space is **two stacked single-modality spaces**, not a joint shared space. Trailing factors beyond the configured active set are numerically dead.

### Interpretation question

What would a genuinely shared factor look like, and what would have to be true of the data for MOFA to find one? Candidates: view scaling, sparsity/dimensionality mismatch, and ATAC gene activity being a smoothed proxy rather than a measurement.

Factor dominance (active factors):

- `Factor1`: dominant=`atac_cbg` (rna=0.0002, atac_cbg=0.7950, nondominant=0.00024)
- `Factor2`: dominant=`rna` (rna=0.5163, atac_cbg=0.0002, nondominant=0.00018)
- `Factor3`: dominant=`rna` (rna=0.1715, atac_cbg=0.0003, nondominant=0.00026)
- `Factor4`: dominant=`rna` (rna=0.1456, atac_cbg=0.0002, nondominant=0.00022)
- `Factor5`: dominant=`rna` (rna=0.1281, atac_cbg=0.0002, nondominant=0.00019)
- `Factor6`: dominant=`rna` (rna=0.1118, atac_cbg=0.0002, nondominant=0.00019)
- `Factor7`: dominant=`rna` (rna=0.0797, atac_cbg=0.0002, nondominant=0.00019)
- `Factor8`: dominant=`atac_cbg` (rna=0.0004, atac_cbg=0.0585, nondominant=0.00035)
- `Factor9`: dominant=`atac_cbg` (rna=0.0004, atac_cbg=0.0553, nondominant=0.00042)
- `Factor10`: dominant=`rna` (rna=0.0431, atac_cbg=0.0001, nondominant=0.00013)
- `Factor11`: dominant=`rna` (rna=0.0368, atac_cbg=0.0001, nondominant=0.00012)
- `Factor12`: dominant=`rna` (rna=0.0254, atac_cbg=0.0001, nondominant=0.00006)
- `Factor13`: dominant=`rna` (rna=0.0253, atac_cbg=0.0001, nondominant=0.00007)
- `Factor14`: dominant=`atac_cbg` (rna=0.0001, atac_cbg=0.0237, nondominant=0.00013)
- `Factor15`: dominant=`rna` (rna=0.0216, atac_cbg=0.0001, nondominant=0.00007)
- `Factor16`: dominant=`rna` (rna=0.0171, atac_cbg=0.0000, nondominant=0.00005)
- `Factor17`: dominant=`rna` (rna=0.0139, atac_cbg=0.0000, nondominant=0.00005)
- `Factor18`: dominant=`rna` (rna=0.0031, atac_cbg=0.0000, nondominant=0.00004)


## Do the factors track cell identity?

- `atac_leiden`: strongest factor `Factor2` eta_squared=0.179
- `azimuth_label`: strongest factor `Factor6` eta_squared=0.562
- `leiden_wnn`: strongest factor `Factor15` eta_squared=0.903
- `rna_leiden`: strongest factor `Factor10` eta_squared=0.883


## Cluster concordance (ARI / NMI)

- atac_leiden vs azimuth_label: ARI=0.0117, NMI=0.0283 (k_a=30, k_b=34, n=21978)
- atac_clusters vs azimuth_label: ARI=0.0069, NMI=0.0106 (k_a=25, k_b=34, n=21978)
- leiden_wnn vs rna_leiden: ARI=0.2469, NMI=0.4430 (k_a=23, k_b=20, n=21978)
- leiden_wnn vs atac_leiden: ARI=0.0201, NMI=0.0398 (k_a=23, k_b=30, n=21978)
- rna_leiden vs azimuth_label: ARI=0.1801, NMI=0.2622 (k_a=20, k_b=34, n=21978)
- leiden_wnn vs azimuth_label: ARI=0.0655, NMI=0.2007 (k_a=23, k_b=34, n=21978)
- atac_leiden vs rna_leiden: ARI=0.0300, NMI=0.0409 (k_a=30, k_b=20, n=21978)


## ATAC QC vs clustering (eta squared)

- log10 nFrags vs `atac_leiden`: eta_squared=0.555 (transform=log10)
- log10 nFrags vs `rna_leiden`: eta_squared=0.102 (transform=log10)
- log10 ReadsInTSS vs `atac_leiden`: eta_squared=0.499 (transform=log10)
- log10 ReadsInTSS vs `rna_leiden`: eta_squared=0.174 (transform=log10)
- DoubletEnrichment vs `atac_leiden`: eta_squared=0.162 (transform=none)
- DoubletEnrichment vs `rna_leiden`: eta_squared=0.052 (transform=none)
- PromoterRatio vs `atac_leiden`: eta_squared=0.141 (transform=none)
- PromoterRatio vs `rna_leiden`: eta_squared=0.094 (transform=none)

## Gene-level bridge (MOFA feature space)

- Bridge metric: `feature_std_on_mofa_inputs`
- Shared genes: 1962
- Pearson r (log1p SD): -0.046
- Spearman rho (log1p SD): -0.04

## Focus markers

RNA and ATAC use **different units** -- interpret concordance, not equality.

- **AGER**: RNA mean=-0.000, ATAC activity mean=-0.000
- **SCGB1A1**: RNA mean=0.000, ATAC activity mean=0.000
- **LUM**: RNA mean=0.000, ATAC activity mean=-0.000
- **ACTA2**: RNA mean=0.000, ATAC activity mean=-0.000

### Markers not recovered in both modality spaces

- **SFTPC**: Not recovered in both modalities used by Module 3
- **SFTPA2**: Not recovered in both modalities used by Module 3
- **FOXJ1**: Not recovered in both modalities used by Module 3
- **PECAM1**: Not recovered in both modalities used by Module 3
- **CLDN5**: Not recovered in both modalities used by Module 3
- **PTPRC**: Not recovered in both modalities used by Module 3
- **NKG7**: Not recovered in both modalities used by Module 3
- **PDGFRB**: Not recovered in both modalities used by Module 3
- **EPCAM**: Documented omission: not present in both MOFA/MuData ATAC gene spaces (e.g. RNA-only epithelial marker)

## Focus markers (paired correlations)

Gene-level RNA<->ATAC correlation is weak everywhere in this dataset (largest absolute r typically << 0.2). Restricting to one cell type usually **shrinks** r by removing between-type variance -- that is the stricter test.

### Global (all nuclei)

- **AGER**: r~-0.006 (n=21978; rna_mean=-0.0000)
- **SCGB1A1**: r~0.002 (n=21978; rna_mean=0.0000)
- **LUM**: r~-0.008 (n=21978; rna_mean=0.0000)
- **ACTA2**: r~0.002 (n=21978; rna_mean=0.0000)

### Within expected Azimuth cell type(s)

- **AGER**: r~-0.031 (n=2210; rna_mean=0.2890; AT1)
- **SCGB1A1**: r~0.108 (n=37; rna_mean=4.8506; Transitional Club-AT2) -- **n_nuclei=37 below threshold 50; correlation not interpretable**
- **LUM**: r~-0.036 (n=1750; rna_mean=0.1990; Alveolar fibroblasts;Adventitial fibroblasts)
- **ACTA2**: r~0.048 (n=343; rna_mean=0.7856; Myofibroblasts;Smooth muscle;SM activated stress response)

## Interpretation prompts

1. What evidence shows this is true multiome rather than unpaired same-donor assays?
2. Why is this MOFA space not a joint latent space? Quote the max non-dominant R^2.
3. What would a shared factor look like, and what data properties might prevent one?
4. Why does within-cell-type correlation usually shrink relative to global?
5. When rna_mean==0 for a defining marker (e.g. NKG7 in NK/CD8), what hypotheses remain?

## Recorded layer probe (Q1 / Q11 / M12)

- RNA layers present: ['spliced', 'spliced_unspliced_sum', 'unspliced']
- Encoding / layout: CSR (`encoding-type=csr_matrix`); Q11 confirmed
- Preferred teaching layer when MuData is loaded: `spliced_unspliced_sum`
- NKG7 column sums (all nuclei): `{"spliced__ENSG00000105374.10": {"sum": 72.39434814453125, "nnz": 21, "max": 4.655576705932617, "mean": 0.00479973154142499}, "spliced_unspliced_sum__ENSG00000105374.10": {"sum": 29.0, "nnz": 26, "max": 2.0, "mean": 0.0019226943841204047}, "unspliced__ENSG00000105374.10": {"sum": 11.5, "nnz": 14, "max": 2.0, "mean": 0.0007624477730132639}}`
- **M12 finding:** Layer switch does not rescue NKG7 in Azimuth NK/CD8 nuclei. Global NKG7 is non-zero in all layers, so the within-type zero is a label-vs-expression discrepancy to report openly (not a silent nan).

## Three-sentence results draft (edit me)

In Donor_2 (`HBM994.NXZC.854`), SNARE-seq2 product `HBM828.GPVG.252` provides paired RNA and ATAC for the same nuclei. MOFA factors are modality-private: the maximum non-dominant R^2 is 0.0004 percent (units match percent-scale stored totals) -- two stacked single-modality spaces, not a joint embedding. Focus-gene RNA<->ATAC correlations are near zero globally and usually smaller within cell type; zeros and low-n groups are flagged explicitly rather than reported as bare nan.
