# Data description

This folder contains the public data files used by the analysis scripts.

## GSE150318

`GSE150318/GSE150318_counts.csv.gz`

Count matrix for the fish aging RNA-seq dataset. The script also downloads GEO metadata through `GEOquery` when needed.

## GSE79396

`GSE79396/GSE79396_series_matrix.txt.gz`

GEO series matrix used to build the subject-by-probe-by-time tensor for the vaccine response analysis.

`GSE79396/GPL13158.annot.gz`

GEO platform annotation file kept with the dataset for probe annotation checks.

## GSE79396 Cell 2017 supplement

The files under `GSE79396_Cell2017_supplement/` are a small subset of the supplemental package from Li et al. (Cell, 2017). They are used by the response-level tPLS analysis:

- `genetable_vax05_rma_ordered.txt`
- `formated2_IgG_NLSrenalyzed.txt`
- `formatted2_TFH.txt`
- `formatted_BTMs.txt`
- `VZV_age_gender2.txt`
- `BTM_for_GSEA_20131008.gmt`

Only the files required by the script were copied here.
