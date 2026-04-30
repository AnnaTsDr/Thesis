# CD4 T Cell Analysis Scripts

This folder contains scripts specifically focused on CD4+ T cell analysis across different conditions and age groups.

## Scripts

### Core CD4 Analysis
- **CD4_TMS_gsea.R**: Gene Set Enrichment Analysis for CD4 T cells from TMS data
  - Performs GSEA on CD4 T cell subsets
  - Identifies enriched pathways in different CD4 subsets
  - Compares pathway activity across age groups

- **Merge_CD4_Human.R**: Merges human and mouse CD4 T cell data
  - Cross-species comparison of CD4 T cells
  - Identifies conserved CD4 T cell signatures
  - Enables translational insights

### CD4 Subset Analysis
- **cd4_analizis_before_merge.R**: CD4 analysis before dataset merging
  - Pre-merge quality control
  - Initial clustering and characterization
  - Subset identification

- **cd4_big_data_combined.R**: Analysis of combined large CD4 dataset
  - Handles large-scale CD4 T cell data
  - Integrates multiple CD4 datasets
  - Comprehensive subset analysis

- **cd4t0n6.R**: CD4 T cell analysis at time points 0 and 6
  - Time-series analysis of CD4 T cells
  - Dynamic changes in CD4 subsets
  - Temporal gene expression patterns

- **cd4t24.R**: CD4 T cell analysis at time point 24
  - Extended time-point analysis
  - Long-term CD4 T cell dynamics

- **cd4t6.R**: CD4 T cell analysis at time point 6
  - Intermediate time-point analysis
  - Mid-term CD4 T cell changes

## CD4 T Cell Subsets Analyzed

The scripts analyze various CD4 T cell subsets including:
- Naive CD4 T cells
- Effector memory (TEM) CD4 T cells
- Regulatory T cells (Tregs: aTregs and rTregs)
- Cytotoxic CD4 T cells
- Exhausted CD4 T cells
- Naive_ISG15 CD4 T cells

## Key Analyses

- Identification of age-related changes in CD4 subsets
- Differential gene expression between CD4 subsets
- Pathway analysis of CD4 T cell function
- Cross-species conservation analysis
- Temporal dynamics of CD4 T cells

## Usage

1. Ensure CD4 T cell data is available in `../../data/processed/`
2. Run scripts in order of analysis complexity
3. Start with basic subset analysis before advanced integrative analysis

## Dependencies

- Seurat
- clusterProfiler
- dplyr
- ggplot2
- clustree
- Other specialized packages

## Output

Results are saved to:
- Processed CD4 data in `../../data/processed/`
- CD4-specific results in `../../data/results/`
- CD4 figures in `../../figures/`
