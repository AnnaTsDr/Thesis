# External Data Directory

This directory contains external datasets, reference files, and third-party data used in the immunaging analysis.

## Contents

### Public Datasets
- Tabula Muris Senis (TMS) data
- GEO dataset downloads
- Reference genomes
- Annotation files

### Reference Materials
- Gene ontology databases
- Pathway databases
- Cell type annotations
- Species conversion files

### Analysis Support
- External analysis tools
- Comparison datasets
- Validation data
- Literature datasets

## Subdirectories

### Centenarians
- Supercentenarians-specific external data
- Comparison datasets for exceptional aging
- Related publication data

### GSE Datasets
- **GSE141540_RAW**: Raw data from GEO dataset
- **GSE149652**: Additional GEO dataset
- **GSE157277_RAW**: Raw data from another GEO study
- Other downloaded datasets for validation

### Tabula Muris Senis
- Official TMS data files
- Raw count matrices
- Processed annotations
- Metadata files

## Data Sources

### Primary Sources
- **Tabula Muris Senis Consortium**: Main single-cell dataset
- **Gene Expression Omnibus (GEO)**: Public gene expression data
- **Reactome**: Pathway database
- **Gene Ontology**: Functional annotation
- **Ensembl/BioMart**: Gene annotation and conversion

### Reference Files
- **HOM_MouseHumanSequence.txt**: Mouse-human gene conversion
- **Mus_musculus.GRCm38.84.gtf**: Mouse genome annotation
- **cell_ontology_class_functional_annotation**: Cell type annotations
- **Senescensce_genes.gmt**: Senescence gene sets

## File Formats

### Data Files
- `.h5ad`: AnnData format (Python single-cell)
- `.h5seurat`: Seurat HDF5 format
- `.rds`: R data files
- `.csv/.tsv`: Tabular data
- `.gtf`: Genome annotation
- `.gmt`: Gene set format

### Compressed Files
- `.tar.gz`: Compressed archives
- `.zip`: ZIP archives

## Usage Guidelines

### Data Access
```r
# Read external data
data <- readRDS("external_data/filename.rds")
data <- read.csv("external_data/filename.csv")
```

### Data Integration
- Use appropriate conversion tools
- Check coordinate systems and versions
- Validate data quality
- Document data sources

## Data Citation

When using external data, always cite:
- Original publication
- Database or repository
- Accession numbers
- Download date
- Data version

## Data Quality

### Validation
- Check data integrity
- Verify file formats
- Validate annotations
- Cross-check with metadata

### Documentation
- Record data source
- Note any processing steps
- Document transformations
- Track data versions

## Storage Considerations

- External data can be large
- Consider compression for storage
- Use appropriate file formats
- Document data organization
- Regular backups of irreplaceable data

## Licensing and Usage

- Respect data usage licenses
- Follow citation requirements
- Acknowledge data sources
- Comply with usage restrictions
- Share derived data appropriately

## Related Directories

- Raw project data: `../data/raw/`
- Analysis scripts: `../scripts/`
- Documentation: `../documentation/`

## Data Management

### Download Records
- Source URL
- Download date
- File checksums
- Accession numbers
- Version information

### Processing Records
- Conversion steps
- Quality control measures
- Integration methods
- Any modifications made

## Best Practices

1. **Citation**: Always cite external data sources
2. **Documentation**: Keep detailed records of data sources
3. **Validation**: Verify data quality before use
4. **Backup**: Maintain backups of downloaded data
5. **Version Control**: Track data versions
6. **Compliance**: Follow usage licenses and restrictions
