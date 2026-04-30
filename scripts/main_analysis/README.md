# Main Analysis Scripts

This folder contains the core analysis scripts for the immunaging thesis project.

## Scripts

### Core Analysis Pipeline
- **TMS_merge.R**: Merges tissue-specific datasets from the Tabula Muris Senis (TMS) project
  - Reads clustered and annotated tissue datasets
  - Performs integrated analysis across 16 tissues
  - Saves merged dataset for downstream analysis

- **TMS_R_Analysis.R**: Main analysis pipeline for merged TMS data
  - Performs dimensionality reduction (PCA, t-SNE, UMAP)
  - Identifies highly variable genes
  - Conducts clustering analysis

- **TMS_analysis_without_spliting.R**: Alternative analysis approach without tissue splitting
  - Provides comparative analysis methodology
  - Useful for different normalization strategies

- **Senescence_integrated.R**: Integrated senescence analysis across tissues
  - Identifies senescence-associated genes
  - Compares senescence patterns across age groups
  - Performs pathway analysis

### Utility Scripts
- **AssessNodesAdaptedfromSeurat.R**: Node assessment for clustering
- **Create_table_for_correlation.R**: Creates correlation tables
- **DEG_overlap_percent.R**: Calculates differential expression overlap percentages
- **cell_ontology_class.R**: Cell ontology classification
- **add_zero_count_CD4_subsets.R**: Handles zero counts in CD4 subset analysis
- **alternative_analysis.R**: Alternative analysis methodologies
- **Spleen_reclustering.R**: Re-clustering of spleen data
- **Tissues_TMS.R**: Tissue-specific TMS analysis
- **Tissues_TMS_Senescence_Reactome.R**: Reactome pathway analysis for senescence
- **United_Seurat_monocle.R**: Combined Seurat and Monocle analysis
- **create_reference.R**: Creates reference datasets
- **Hezis_data_new_ref.R**: Reference data processing
- **mouse_to_human_genes.R**: Mouse to human gene conversion
- **pearsons_corr.R**: Pearson correlation analysis
- **pearsons_corr_2.R**: Additional correlation analysis
- **ctl_merge.R**: Control dataset merging
- **Utilities_1.R**: General utility functions
- **Venn_Diagram.R**: Venn diagram generation
- **Venn_Diagram_2.R**: Alternative Venn diagram approach
- **Venn_Diagram_Senescence_age_groups_DGE.R**: Venn diagrams for senescence DGE across age groups

## Usage

Run scripts in the following order:
1. Start with `TMS_merge.R` to create the merged dataset
2. Run `TMS_R_Analysis.R` for main analysis
3. Use specialized scripts for specific analyses

## Dependencies

All scripts require:
- Seurat package
- dplyr
- ggplot2
- Other specialized packages (see individual scripts)

## Output

Results are saved to:
- `../../data/processed/` for processed data objects
- `../../data/results/` for analysis results
- `../../figures/` for generated plots
