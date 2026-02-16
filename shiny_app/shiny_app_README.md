# COVID-19 Lung Atlas Interactive Explorer

An R Shiny application for exploring single-cell RNA sequencing data from COVID-19 lung tissue.

## Data Source

**Melms et al. (2021)** "A molecular single-cell lung atlas of lethal COVID-19" *Nature*

- GEO Accession: [GSE171524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524)
- DOI: [10.1038/s41586-021-03569-1](https://doi.org/10.1038/s41586-021-03569-1)

## Dataset Summary

- **94,027 cells** from 27 samples
- **19 COVID-19 patients**, 7 healthy controls
- **9 major cell types** identified
- Single-nucleus RNA sequencing (snRNA-seq)

## Features

### 1. UMAP Explorer
- Interactive UMAP visualization
- Color by cell type, condition, sample, or cluster
- Split view by condition

### 2. Gene Expression
- Visualize gene expression on UMAP
- Violin plots comparing COVID vs Control
- Dot plots showing expression across cell types
- Key genes: NEAT1, MALAT1, AXL, KRT8, CTHRC1, IL1B, IL6, etc.

### 3. Differential Expression
- Interactive volcano plots
- Filter by cell type (global or cell type-specific)
- Adjustable significance thresholds
- Searchable DEG table

### 4. Cell Proportions
- Bar charts comparing proportions by condition
- Pie chart of overall distribution
- Statistical comparison table

### 5. Cell-Cell Interactions
- Ligand-receptor interaction analysis
- Filter by pathway (TGFb, IL1, Chemokine, etc.)
- Visualization of interaction changes

## Installation

### Requirements
```r
install.packages(c("shiny", "shinydashboard", "plotly", "DT", "dplyr", "ggplot2", "viridis", "tidyr"))
```

### Running the App
```r
# Set working directory to shiny_app folder
setwd("path/to/shiny_app")

# Run the app
shiny::runApp()
```

Or from RStudio, open `app.R` and click "Run App".

## Data Files

- `data/umap_sampled.csv` - UMAP coordinates and cell metadata (20,000 cells)
- `data/gene_expression.csv` - Expression matrix for key genes
- `data/de_results.csv` - Global differential expression results
- `data/de_celltype.csv` - Cell type-specific DE results
- `data/interactions.csv` - Cell-cell interaction analysis

## Key Findings

1. **Myeloid cells** significantly increased in COVID-19 (p < 0.0001)
2. **Epithelial cells** (AT1/AT2) decreased (p = 0.0007)
3. **Pathological fibroblasts** (CTHRC1+) expanded 2.5x
4. **NEAT1/MALAT1** upregulated in macrophages
5. **Impaired efferocytosis** (AXL/MERTK signaling decreased)
6. **Fibrosis signatures** elevated (COL1A1, POSTN, TNC)

## Citation

If you use this data or application, please cite:

```
Melms JC, Biermann J, Huang H, et al. A molecular single-cell lung atlas of lethal COVID-19. 
Nature. 2021;595(7865):114-119. doi:10.1038/s41586-021-03569-1
```
