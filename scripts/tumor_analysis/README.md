# Tumor/Immunotherapy Analysis Scripts

This folder contains analysis scripts related to tumor immunology and immunotherapy responses.

## Scripts

### Anti-PD-L1 Treatment Analysis
- **Anti_PD_L1_A_tumor.R**: Analysis of anti-PD-L1 treatment in tumor model A
- **Anti_PD_L1_B_tumor.R**: Analysis of anti-PD-L1 treatment in tumor model B
- **Anti_PD_L1_C_tumor.R**: Analysis of anti-PD-L1 treatment in tumor model C
- **Anti_PD_L1_D_tumor.R**: Analysis of anti-PD-L1 treatment in tumor model D

These scripts analyze:
- Tumor microenvironment changes
- Immune cell infiltration
- Treatment response markers
- Differential gene expression

### Chemotherapy Analysis
- **Chemo_tumor.R**: Analysis of chemotherapy treatment effects on tumors
  - Evaluates chemotherapy impact on immune cells
  - Identifies chemotherapy response signatures
  - Compares with immunotherapy responses

### Bladder Cancer Analysis
- **Bladdet_cancer.R**: Bladder cancer specific analysis
  - Tumor-specific immune profiling
  - Identification of bladder cancer biomarkers
  - Comparison with aging-related immune changes

### Untreated Controls
- **Untreated_A_tumor.R**: Control analysis for tumor model A
- **Untreated_B_tumor.R**: Control analysis for tumor model B
- **tumor.R**: General tumor analysis script

## Usage

Each script can be run independently for specific tumor models:
1. Ensure required data files are available
2. Set appropriate working directory
3. Run specific tumor model script

## Key Analyses

- Clustering of tumor-infiltrating immune cells
- Differential expression between treated and untreated
- Identification of treatment-responsive cell populations
- Pathway analysis of treatment effects

## Dependencies

- Seurat
- dplyr
- ggplot2
- EnhancedVolcano
- Other standard R packages

## Output

Results are saved to:
- Processed data objects in `../../data/processed/`
- Analysis results in `../../data/results/`
- Figures in `../../figures/`
