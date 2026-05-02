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
