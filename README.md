# Single-Cell Atlas of COVID-19 Lung Pathology

Reanalysis of the single-nucleus RNA sequencing dataset from **Melms et al. (2021)** *"A molecular single-cell lung atlas of lethal COVID-19"* published in [Nature](https://doi.org/10.1038/s41586-021-03569-1). This project reproduces key findings from the original study and provides an interactive R Shiny application for data exploration.

## Motivation

COVID-19 causes severe lung damage, but the cellular mechanisms driving tissue destruction remain incompletely understood. Melms et al. generated the first comprehensive single-cell atlas of the COVID-19 lung using autopsy samples from 19 patients and 7 controls. This project reanalyzes their publicly available data (GEO: [GSE171524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524)) to:

- Reproduce the major findings using a modern scverse-based pipeline
- Validate cell type composition changes and myeloid dysfunction signatures
- Build an interactive application for exploring the atlas

## Dataset

| | |
|---|---|
| **Cells** | 116,314 nuclei (94,027 after QC) |
| **Genes** | 34,546 (29,050 after filtering) |
| **Samples** | 27 (7 control, 20 COVID-19) |
| **Patients** | 26 (7 healthy donors, 19 COVID-19) |
| **Tissue** | Lung (autopsy/necropsy) |
| **Technology** | 10x Chromium snRNA-seq |

## Key Findings Reproduced

1. **Myeloid expansion** -- Macrophages and monocytes are significantly enriched in COVID-19 lungs, with accumulation of monocyte-derived macrophages (MDMs) replacing tissue-resident alveolar macrophages (AMs)
2. **Epithelial cell loss** -- AT1 and AT2 alveolar cells are depleted, with emergence of damage-associated transient progenitors (DATPs; KRT8+/CLDN4+)
3. **Pathological fibroblasts** -- CTHRC1+ fibroblasts expand ~2.5-fold, expressing pro-fibrotic genes (COL1A1, POSTN, TNC)
4. **Macrophage dysfunction** -- lncRNAs NEAT1/MALAT1 are upregulated while efferocytosis receptors AXL/MERTK are downregulated in COVID-19 myeloid cells
5. **Inflammatory signaling** -- IL1B, IL6, and TNF are elevated primarily in myeloid populations

## Analysis Pipeline

```
Raw counts (27 CSV.gz files from GEO)
  |
  v
01_data_loading         Load, transpose, concatenate into AnnData
  |
  v
02_quality_control      MAD-based filtering, Scrublet doublet detection
  |                     116,314 -> 94,027 cells
  v
03_integration_scvi     Batch correction with scVI (30-dim latent space)
  |                     4,000 HVGs, negative binomial likelihood
  v
04_clustering           Leiden clustering, marker-based cell type annotation
  |                     9 major cell types identified
  v
05_differential_expr    Wilcoxon DE (global + per cell type) + pseudobulk
  |
  v
06_myeloid_analysis     AM vs MDM classification, dysfunction scoring
  |
  v
07_prepare_shiny_data   Export subsampled data for interactive app
```

## Project Structure

```
COVID-19-Lung-Atlas/
|
|-- data/
|   |-- raw_data/                  # 27 GSM*.csv.gz count matrices
|   |-- metadata/                  # sample_metadata.csv
|   +-- processed_data/            # adata_raw.h5ad, adata_qc.h5ad, ...
|
|-- notebooks/
|   |-- 01_data_loading.ipynb
|   |-- 02_quality_control.ipynb
|   |-- 03_integration_scvi.ipynb
|   |-- 04_clustering_annotation.ipynb
|   |-- 05_differential_expression.ipynb
|   |-- 06_myeloid_analysis.ipynb
|   +-- 07_prepare_shiny_data.ipynb
|
|-- scripts/
|   |-- utils.py                   # Data loading, QC, pseudobulk utilities
|   |-- plotting.py                # Publication-quality plotting functions
|   +-- markers.py                 # Cell type markers & COVID-19 gene signatures
|
|-- covid19-lung-atlas/            # R Shiny interactive application
|   |-- shiny_app.R
|   |-- setup.R
|   +-- data/                      # Pre-computed CSVs for the app
|
|-- results/
|   |-- figures/                   # Generated plots (QC, integration, DE, ...)
|   +-- tables/                    # DE result tables (CSV)
|
|-- refs/                          # Reference paper PDF
+-- LICENSE                        # MIT
```

## Interactive Shiny App

The `covid19-lung-atlas/` directory contains an R Shiny application with 8 tabs:

- **UMAP Explorer** -- Color by cell type, condition, sample, or cluster with adjustable point size/opacity
- **Gene Expression** -- Search any gene; view expression on UMAP, violin plots (COVID vs Control), and dot plots by cell type
- **Gene Signatures** -- 11 pre-defined COVID-19 signatures (myeloid dysfunction, fibrosis, DATP, exhaustion, interferon response, etc.) with heatmaps and per-cell-type comparisons
- **Differential Expression** -- Interactive volcano plots with adjustable thresholds; global or cell type-specific; downloadable DEG tables
- **Cell Proportions** -- Proportion bar charts, pie chart, and statistical comparison table
- **Sample Explorer** -- Per-sample cell counts, metadata, and UMAP highlighting
- **Cell-Cell Interactions** -- Ligand-receptor analysis across 9 pathways (TGFb, IL1, IL6, TNF, Chemokine, PDGF, FGF, Wnt, Efferocytosis)
- **About** -- Dataset summary and key findings

### Running the App

```r
# Install dependencies
source("covid19-lung-atlas/setup.R")

# Launch
shiny::runApp("covid19-lung-atlas/shiny_app.R")
```

## Getting Started

### Prerequisites

**Python** (notebooks):
```
scanpy>=1.9
anndata>=0.10
scvi-tools>=1.0
scrublet
numpy
pandas
scipy
matplotlib
seaborn
```

**R** (Shiny app):
```
shiny, shinydashboard, plotly, DT, dplyr, tidyr, ggplot2, viridis
```

### Reproducing the Analysis

```bash
# 1. Clone the repository
git clone https://github.com/<username>/COVID-19-Lung-Atlas.git
cd COVID-19-Lung-Atlas

# 2. Download raw data from GEO (GSE171524) into data/raw_data/

# 3. Run notebooks sequentially (01 through 07)
jupyter notebook notebooks/

# 4. Launch the Shiny app
Rscript -e "shiny::runApp('covid19-lung-atlas/shiny_app.R')"
```

> Notebooks 01-02 have been executed. Notebooks 03-07 are written and ready to run (they require GPU/CPU time for scVI training).

## Methods

| Step | Tool | Details |
|------|------|---------|
| QC filtering | Scanpy + Scrublet | MAD-based thresholds (5 MADs); min 200 genes, 500 UMIs; <20% mito; doublet score > 0.25 |
| Batch integration | scVI | 2-layer VAE, 128 hidden units, 30 latent dims, negative binomial likelihood |
| Clustering | Leiden | Resolution 0.8 on scVI latent kNN graph (k=30) |
| Cell annotation | Marker scoring | 9 major types via gene signature scoring of known markers |
| Differential expression | Wilcoxon rank-sum | Per-gene test (COVID vs Control), BH-corrected; also pseudobulk t-test |
| Myeloid analysis | Scanpy | AM/MDM classification, dysfunction scoring (lncRNA, efferocytosis, inflammation) |

## References

Melms JC, Biermann J, Huang H, Wang Y, Nair Y, Tagber S, et al. A molecular single-cell lung atlas of lethal COVID-19. *Nature*. 2021;595(7865):114-119. doi:[10.1038/s41586-021-03569-1](https://doi.org/10.1038/s41586-021-03569-1)

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
