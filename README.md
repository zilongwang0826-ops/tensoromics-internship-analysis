# TensorOmics internship analyses

This repository contains the scripts, input data, saved results and report files for my internship project on tensor-based analysis of longitudinal omics data.

The work uses two public transcriptomic datasets:

- `GSE150318`: longitudinal RNA-seq data from fin biopsies of the short-lived killifish *Nothobranchius furzeri*.
- `GSE79396`: longitudinal human PBMC transcriptomic data after herpes zoster vaccination, together with selected response measurements from the associated Cell 2017 study.

## Repository structure

```text
data/       Public input files used by the scripts
scripts/    R and R Markdown analysis scripts
results/    Saved tables and figures produced by the analyses
report/     LaTeX report source, figures and compiled report PDF
```

## Scripts

The scripts are kept in the order in which they were used:

1. `scripts/01_GSE150318_fish_tPLSDA.Rmd`
   Fish aging analysis using tensor-based PLS-DA.
2. `scripts/02_GSE79396_vaccine_tPCA.R`
   Tensor PCA analysis of the vaccine transcriptome trajectories.
3. `scripts/03_GSE79396_response_tPLS.R`
   Response-level tPLS analysis for IgG and CXCR3+ Tfh outcomes.
4. `scripts/04_trajectory_analysis.R`
   Additional trajectory-level checks used for the final report figures and tables.

## Running the analyses

The scripts were run in R. Main packages include `tidyverse`, `data.table`, `GEOquery`, `ggplot2`, `tensorOmics`, `mixOmics`, `pls` and `ggrepel`.

Example commands:

```r
rmarkdown::render("scripts/01_GSE150318_fish_tPLSDA.Rmd")
source("scripts/02_GSE79396_vaccine_tPCA.R")
source("scripts/03_GSE79396_response_tPLS.R")
source("scripts/04_trajectory_analysis.R")
```

Some scripts may download GEO metadata if it is not already present locally.

## Data note

The data files are public research data from GEO and from the supplementary material of Li et al. (Cell, 2017). They are included here so that the analyses can be rerun with the same local inputs used for the report.

