# Senescence Analysis Scripts

This folder contains scripts focused on cellular senescence analysis across different tissues and age groups.

## Scripts

### Differential Expression Analysis
- **Senescence_deg.R**: Differential expression analysis of senescence markers
  - Identifies genes differentially expressed in senescent cells
  - Compares senescence across age groups
  - Tissue-specific senescence signatures

- **DEG_analysis_senescent_vs_nonsenescent.R**: DEG analysis between senescent and non-senescent cells
  - Direct comparison of senescent vs non-senescent populations
  - Identifies key senescence-associated genes
  - Statistical analysis of senescence markers

### Marker Analysis
- **Senescence_pos_markers_between_age_groups.R**: Positive senescence markers between age groups
  - Identifies markers that increase with age
  - Compares young vs old senescence patterns
  - Age-dependent senescence signatures

- **Senescence_pos_markers_senescent_vs_nonsenescent.R**: Positive markers for senescent vs non-senescent cells
  - Cell-autonomous senescence markers
  - Identifies definitive senescence signatures
  - Validates senescence classification

### Expression Analysis
- **AverageExpression_Senescence_vs_nonsenescence.R**: Average expression analysis
  - Calculates average expression of senescence genes
  - Compares expression levels across conditions
  - Provides quantitative senescence metrics

### Pathway Analysis
- **senescence_gsea.R**: Gene Set Enrichment Analysis for senescence
  - Pathway analysis of senescence genes
  - Identifies biological processes involved in senescence
  - Reactome and KEGG pathway analysis

## Senescence Markers Analyzed

Key senescence markers examined include:
- p16 (CDKN2A)
- p21 (CDKN1A)
- p53 (TP53)
- SA-beta-gal associated genes
- Inflammatory markers (IL6, IL1B, TNF)
- Other senescence-associated secretory phenotype (SASP) genes

## Key Analyses

- Identification of tissue-specific senescence patterns
- Age-dependent changes in senescence markers
- Conservation of senescence pathways across tissues
- Quantification of senescence burden
- Pathway analysis of senescence processes

## Usage

1. Ensure senescence-annotated data is available
2. Run DEG analysis scripts first
3. Follow with marker and pathway analysis
4. Use expression analysis for validation

## Dependencies

- Seurat
- clusterProfiler
- ReactomePA
- dplyr
- ggplot2
- EnhancedVolcano
- Other enrichment analysis packages

## Output

Results are saved to:
- Senescence DEG tables in `../../data/results/`
- Pathway analysis results in `../../data/results/`
- Senescence figures in `../../figures/`
