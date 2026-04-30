# Figures Directory

This directory contains all figures and plots generated during the immunaging analysis.

## Figure Types

This folder includes various types of visualizations:

### Dimensionality Reduction Plots
- UMAP plots showing cell distributions
- t-SNE plots for cluster visualization
- PCA plots for variance explanation

### Expression Visualizations
- Heatmaps showing gene expression patterns
- Violin plots for gene expression distributions
- Feature plots showing expression on UMAP
- Dot plots for comparative expression

### Clustering Visualizations
- Clustree plots showing clustering at different resolutions
- Dendrograms showing cluster relationships
- Cluster tree plots

### Differential Expression
- Volcano plots showing DE results
- MA plots for expression comparisons
- Venn diagrams showing gene overlap

### Population Analysis
- Pie charts showing cell type proportions
- Bar plots for subset comparisons
- Stacked bar plots for composition

### Pathway Analysis
- Pathway enrichment plots
- Reactome pathway visualizations
- GSEA plots

### Specialized Figures
- Pseudotime trajectory plots
- Correlation heatmaps
- Senescence-specific visualizations
- Age-related change plots

## File Formats

Figures are saved in various formats:
- `.png` - Standard raster format (300 DPI)
- `.pdf` - Vector format for publications
- `.tiff` - High-quality raster format
- `.jpeg` - Compressed raster format

## Naming Conventions

Figures follow descriptive naming:
- `[analysis_type]_[description]_[optional_info].format`
- Examples: `UMAP_CD4_T_cells.png`, `Volcano_senescence_DEGs.pdf`

## Organization

Figures are organized by:
- Analysis type (clustering, DEG, pathway, etc.)
- Cell type (CD4, CD8, general immune cells)
- Comparison (age groups, tissues, conditions)
- Figure type (UMAP, heatmap, volcano, etc.)

## Usage

### Including in Publications
1. Choose appropriate format (PDF for vector, PNG/TIFF for raster)
2. Check resolution requirements
3. Ensure colorblind accessibility
4. Include figure legends in manuscript

### Figure Modification
- Original R scripts are in `scripts/visualization/`
- Modify scripts to regenerate figures with different parameters
- Use vector formats when possible for easy modification

## Best Practices

1. **Resolution**: Use at least 300 DPI for publications
2. **Color**: Use colorblind-friendly palettes
3. **Labels**: Ensure all axes and legends are clearly labeled
4. **Consistency**: Use consistent styling across related figures
5. **Backup**: Keep original scripts to regenerate figures
6. **Organization**: Use descriptive names for easy identification

## Figure Quality

All figures are generated at publication quality:
- High resolution (300 DPI or higher)
- Clear typography
- Appropriate figure sizes
- Professional color schemes
- Publication-ready formatting

## Related Files

- Generation scripts: `../scripts/visualization/`
- Analysis results: `../data/results/`
- Manuscript figures: `../documentation/`
