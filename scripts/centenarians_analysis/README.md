# Centenarians Analysis Scripts

This folder contains analysis scripts focused on supercentenarians and exceptional aging.

## Scripts

### Supercentenarians Analysis
- **Supercentenarians_T_cells_clusters.R**: T cell clustering in supercentenarians
  - Identifies T cell subsets in supercentenarians
  - Compares with regular aging patterns
  - Characterizes unique immune features

- **Supercentenerians.R**: General supercentenarians analysis
  - Comprehensive analysis of supercentenarian immune profiles
  - Comparison across different age groups
  - Identification of protective immune mechanisms

## Research Focus

This analysis examines:
- **Exceptional Aging**: Individuals who live to 100+ years
- **Immune Profile**: Characterization of immune system in extreme old age
- **Protective Mechanisms**: Identification of factors that promote longevity
- **Comparison**: Supercentenarians vs young vs regular old adults

## Key Comparisons

- Supercentenarians vs Young adults
- Supercentenarians vs Old adults (60-80 years)
- Supercentenarians vs Centenarians (100-109 years)
- Identification of unique signatures in exceptional aging

## Analyses Performed

- T cell subset distribution in supercentenarians
- Differential gene expression vs other age groups
- Pathway analysis of longevity-associated genes
- Identification of protective immune signatures
- Senescence burden in exceptional aging

## Usage

1. Ensure supercentenarian data is available in `../../data/processed/`
2. Run clustering analysis to identify cell subsets
3. Perform comparative analysis with other age groups
4. Conduct pathway and functional analysis

## Dependencies

- Seurat
- dplyr
- ggplot2
- clusterProfiler
- Other standard analysis packages

## Output

Results are saved to:
- Supercentenarian-specific data in `../../data/processed/`
- Comparative analysis results in `../../data/results/`
- Supercentenarian figures in `../../figures/`

## Significance

Understanding supercentenarian immune profiles provides insights into:
- Successful aging mechanisms
- Potential therapeutic targets for age-related diseases
- Natural models of extreme longevity
- Protective immune signatures
