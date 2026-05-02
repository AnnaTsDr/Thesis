# Thesis

This repository contains data and analysis files for my thesis project.

## External Data

### Mouse Genome Reference
- **Mus_musculus.GRCm38.84.gtf** (857MB) - Mouse genome annotation file (GRCm38 release 84)
  - Tracked with Git LFS
  - Source: Ensembl

### Tabula Muris Senis Data
The following Tabula Muris Senis datasets are used in this project:

- **tabula-muris-senis-droplet-processed-official-annotations.h5ad** (7.7GB)
- **tabula-muris-senis-droplet-official-raw-obj.h5seurat** (7.5GB)
- **tabula-muris-senis-facs-official-raw-obj.h5ad** (2.3GB)
- **tabula-muris-senis-droplet-official-raw-obj.h5ad** (3.8GB)

**Note:** These files exceed Git LFS's 2GB size limit and are stored locally in the `external_data/` directory. They are not tracked in this repository.

### Download Tabula Muris Senis Data
The Tabula Muris Senis project provides a comprehensive resource for cell biology with detailed molecular and cell-type specific portraits of aging.

**Official Repository:** https://github.com/czbiohub-sf/tabula-muris-senis

**Paper:** [Tabula Muris Senis: Ageing cell atlas of Mus musculus](https://www.biorxiv.org/content/10.1101/661728v2)

To download the data files, visit the official Tabula Muris Senis repository or the project's data portal.

### Supercentenarians Data (Centenarians)
Single-cell RNA sequencing data from supercentenarians studying cytotoxic CD4 T cells.

**Research Paper:** [Single-cell transcriptomics reveals expansion of cytotoxic CD4 T cells in supercentenarians](https://doi.org/10.1073/pnas.1907883116)
- Hashimoto K, et al. Proc Natl Acad Sci U S A. 2019 Nov 26;116(48):24242-24251

**Data Source:** GEO Accession GSE141540
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141540

**Note:** Large files (>2GB) are stored locally in `external_data/Centenarians/` and are not tracked in this repository due to Git LFS size limits.

### Bladder Cancer CD4+ T Cells Data
Single-cell RNA and TCR sequencing data from human bladder cancer patients studying intratumoral CD4+ T cell cytotoxicity.

**Research Paper:** [Intratumoral CD4+ T cells mediate anti-tumor cytotoxicity in human bladder cancer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149652)
- Published in PNAS (2020)

**Data Source:** GEO Accession GSE149652
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149652

**Note:** All files from this dataset are tracked with Git LFS as they are under the 2GB size limit.

## Directory Structure

```
.
├── external_data/          # Large data files
├── documentation/          # Project documentation
├── figures/               # Generated figures and plots
├── logs/                  # Analysis logs
├── presentations/         # Presentation materials
└── scripts/               # Analysis scripts
```

## License

This project uses data from Tabula Muris Senis, which is available under a BSD-3-Clause license.
