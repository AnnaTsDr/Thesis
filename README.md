# Immunaging Thesis Repository

This repository contains the complete analysis pipeline and results for a thesis investigating immune system aging and senescence using single-cell RNA sequencing data from the Tabula Muris Senis (TMS) project.

## Overview

This research focuses on understanding how aging affects the immune system, with particular emphasis on:
- T cell subset dynamics across different age groups
- Senescence-associated gene expression patterns
- Tissue-specific immune aging signatures
- CD4+ T cell heterogeneity in aging
- Supercentenarians as a model of exceptional aging

## Repository Structure

```
Immunaging_Thesis/
├── README.md
├── .gitignore
├── Immunaging.Rproj
├── data/
│   ├── raw/              # Raw data files (age-specific RDS files, etc.)
│   ├── processed/        # Processed analysis objects (clustered, normalized, etc.)
│   └── results/          # CSV results, markers, differential expression tables
├── scripts/
│   ├── main_analysis/    # Core analysis scripts (TMS_merge, integrated analysis)
│   ├── tumor_analysis/   # Immunotherapy/tumor analysis (Anti_PD-L1, chemo, bladder cancer)
│   ├── cd4_analysis/     # CD4 T cell specific analysis (GSEA, human data)
│   ├── senescence_analysis/ # Senescence-related analysis (DEGs, markers)
│   ├── centenarians_analysis/ # Centenarians/supercentenarians analysis
│   ├── clustering/       # Cluster-specific analysis scripts
│   └── visualization/    # Plotting and visualization scripts
├── figures/              # Generated plots, heatmaps, clustree images
├── presentations/        # PowerPoint presentations
├── documentation/       # Thesis documents, methods, materials
├── external_data/        # External datasets and references
└── logs/                # Analysis logs
```

## Key Analysis Scripts

### Main Analysis
- **TMS_merge.R**: Merges tissue-specific datasets from Tabula Muris Senis project
- **TMS_R_Analysis.R**: Main analysis pipeline for merged TMS data
- **Senescence_integrated.R**: Integrated senescence analysis across tissues

### CD4 T Cell Analysis
- **CD4_TMS_gsea.R**: Gene Set Enrichment Analysis for CD4 T cells
- **Merge_CD4_Human.R**: Merges human and mouse CD4 T cell data

### Tumor/Immunotherapy Analysis
- **Anti_PD_L1_*.R**: Analysis of anti-PD-L1 treatment effects
- **Chemo_tumor.R**: Chemotherapy treatment analysis
- **Bladdet_cancer.R**: Bladder cancer specific analysis

### Senescence Analysis
- **Senescence_deg.R**: Differential expression analysis of senescence markers
- **DEG_analysis_senescent_vs_nonsenescent.R**: DEG analysis between senescent and non-senescent cells
- **AverageExpression_Senescence_vs_nonsenescence.R**: Average expression analysis

### Centenarians Analysis
- **Supercentenarians_T_cells_clusters.R**: T cell clustering in supercentenarians
- **Supercentenerians.R**: General supercentenarians analysis

## Dependencies

This analysis requires the following R packages:
- Seurat (single-cell analysis)
- dplyr (data manipulation)
- ggplot2 (visualization)
- monocle3 (trajectory analysis)
- clustree (clustering visualization)
- EnhancedVolcano (volcano plots)
- clusterProfiler (functional enrichment)
- And others (see individual scripts for specific requirements)

## Usage

### Running the Main Analysis
1. Ensure all required packages are installed
2. Set working directory to the repository root
3. Run scripts in order:
   - Start with `scripts/main_analysis/TMS_merge.R`
   - Follow with tissue-specific or cell-type-specific analyses as needed

### Data Requirements
- Raw data files should be placed in `data/raw/`
- Processed data will be automatically saved to `data/processed/`
- Results will be saved to `data/results/`

## Key Findings

### T Cell Aging
- Identification of distinct CD4+ T cell subsets that change with age
- Characterization of naive, effector memory, and regulatory T cell dynamics
- Discovery of age-associated cytotoxic CD4+ T cells

### Senescence Patterns
- Tissue-specific senescence signatures
- Conservation of senescence markers across age groups
- Identification of novel senescence-associated genes

### Supercentenarians
- Unique immune profile in exceptional aging
- Comparison with regular aging patterns
- Potential protective mechanisms

## Data Sources

- **Tabula Muris Senis (TMS)**: Main single-cell RNA-seq dataset
- **GEO Datasets**: External validation datasets (GSE141540, GSE149652, GSE157277)
- **Human Data**: CD4 T cell data for cross-species comparison

## Contact

For questions about this analysis, please contact the thesis author.

## License

This repository contains research data and analysis scripts. Please contact the author for permission to use or modify any materials.

## Citation

If you use this analysis or data in your research, please cite the corresponding thesis and the Tabula Muris Senis consortium paper.
