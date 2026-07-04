# Scripts

This folder contains the analysis scripts used for the internship report.

## Main workflow

1. `01_GSE150318_fish_tPLSDA.Rmd`
   Builds the fish subject-by-gene-by-age tensor and applies tensor PLS-DA to compare short-lived and long-lived fish trajectories.

2. `02_GSE79396_vaccine_tPCA.R`
   Builds the vaccine subject-by-probe-by-time tensor and applies tensor PCA to compare young and elderly subjects across post-vaccination time points.

3. `03_GSE79396_response_tPLS.R`
   Uses transcriptome trajectory changes to study IgG and CXCR3+ Tfh response variables from the Li et al. Cell 2017 dataset.

4. `04_trajectory_analysis.R`
   Produces additional trajectory-level figures and summary tables used in the final report.

The scripts are meant to be run from the repository root.

