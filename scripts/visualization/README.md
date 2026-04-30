# Visualization Scripts

This folder contains scripts for generating plots, figures, and visualizations for the immunaging analysis.

## Scripts

### General Plotting
- **Plot_generator.R**: Main plot generation script
  - Generates standard analysis plots
  - Creates publication-quality figures
  - Handles multiple plot types

- **Plot_generator_2.R**: Alternative plot generator
  - Additional visualization options
  - Different styling approaches
  - Supplementary figure generation

- **Plot_generator_3.R**: Extended plot generator
  - Advanced visualization techniques
  - Complex multi-panel figures
  - Custom plot layouts

- **Plot_generator_4.R**: Additional plot generator
  - Specialized visualizations
  - Alternative figure formats
  - Custom plot configurations

### Specific Visualizations
- **Pie_Ploting.R**: Pie chart generation
  - Creates pie charts for cell type proportions
  - Shows subset distributions
  - Compares proportions across conditions

- **DoMultiBarHeatmap.R**: Multi-bar heatmap generation
  - Creates complex heatmaps with multiple bars
  - Shows expression patterns across groups
  - Displays comparative expression data

## Visualization Types

The scripts generate various types of plots including:
- **UMAP/t-SNE plots**: Dimensionality reduction visualizations
- **Heatmaps**: Gene expression patterns
- **Volcano plots**: Differential expression results
- **Violin plots**: Gene expression distributions
- **Feature plots**: Gene expression on UMAP
- **Dot plots**: Expression across groups
- **Bar plots**: Cell type proportions
- **Pie charts**: Subset distributions
- **Clustree plots**: Clustering resolution trees
- **Venn diagrams**: Gene overlap analysis

## Usage

1. Ensure analysis results are available
2. Run appropriate plot generator script
3. Specify output directory (defaults to `../../figures/`)
4. Customize plot parameters as needed

## Customization

Most scripts allow customization of:
- Color schemes
- Figure dimensions
- Font sizes
- Legend placement
- Title and label text
- File formats (PNG, PDF, TIFF, etc.)

## Output

All figures are saved to:
- `../../figures/` by default
- Various formats depending on script settings
- Publication-ready quality (300 DPI or higher)

## Dependencies

- ggplot2
- cowplot
- patchwork
- ComplexHeatmap
- ggrepel
- EnhancedVolcano
- Other visualization packages

## Best Practices

1. Check figure resolution before publication
2. Ensure colorblind-friendly color schemes
3. Include appropriate legends and labels
4. Maintain consistent styling across related figures
5. Save in multiple formats for flexibility
