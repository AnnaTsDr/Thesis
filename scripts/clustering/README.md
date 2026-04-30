# Clustering Analysis Scripts

This folder contains individual cluster analysis scripts for detailed examination of specific cell clusters identified in the main analysis.

## Scripts

### Individual Cluster Analysis
- **Cluster_1.R through Cluster_68.R**: Individual analysis scripts for each cluster
  - Each script analyzes a specific cell cluster
  - Identifies cluster-specific marker genes
  - Performs functional annotation of clusters
  - Examines cluster distribution across conditions

## Clustering Context

These clusters were identified from the main TMS merged dataset analysis, specifically from T cell clustering that identified 47 distinct T cell subsets, which were then further analyzed.

## Cluster Types Analyzed

The scripts cover various cluster types including:
- T cell subsets (CD4, CD8, regulatory, etc.)
- Different activation states
- Tissue-specific clusters
- Age-associated clusters
- Functional subsets (naive, memory, effector, etc.)

## Analysis Performed for Each Cluster

For each cluster, the scripts typically perform:
1. **Marker Identification**: Find cluster-specific marker genes
2. **Differential Expression**: Compare cluster vs other clusters
3. **Functional Annotation**: Identify biological functions
4. **Pathway Analysis**: Enriched pathways in the cluster
5. **Age Distribution**: How cluster changes with age
6. **Tissue Distribution**: Which tissues contain the cluster

## Usage

1. Run the specific cluster script of interest
2. Scripts can be run independently
3. Results are saved with cluster-specific naming
4. Use for detailed characterization of specific cell populations

## Dependencies

- Seurat
- dplyr
- ggplot2
- clusterProfiler
- Other standard analysis packages

## Output

Results are saved to:
- Cluster-specific marker files in `../../data/results/`
- Cluster analysis objects in `../../data/processed/`
- Cluster figures in `../../figures/`

## File Naming Convention

Results follow the pattern:
- `Cluster_[number]` for cluster-specific outputs
- Includes markers, expression data, and visualizations

## Integration with Main Analysis

These cluster scripts complement the main analysis by:
- Providing detailed view of individual clusters
- Enabling focused investigation of specific populations
- Supporting validation of main clustering results
- Facilitating downstream functional analysis
