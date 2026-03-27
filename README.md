# scAutoTune

**scAutoTune** is an R package for automated optimization of feature selection and clustering parameters in **Seurat**-based single-cell RNA-seq analysis.

Single-cell RNA-seq workflows often rely on manual choices such as the number of highly variable genes (HVGs), the number of principal components (PCs), and clustering resolution. These parameters can strongly affect downstream clustering and interpretation. **scAutoTune** provides a reproducible framework to evaluate these choices quantitatively and identify high-performing parameter settings.

## Features

- Estimates an appropriate number of principal components across HVG settings
- Sweeps combinations of HVG counts and clustering resolutions
- Evaluates clustering performance using silhouette scores
- Computes modularity during parameter sweeps
- Uses GAM-based smoothing to identify optimal parameter regions
- Supports optional batch-aware workflows with Harmony integration
- Returns publication-ready visual outputs, including heatmaps and UMAPs

## Installation

Install from GitHub:

```r
# install.packages("devtools")
devtools::install_github("shihaisun-scott/scAutoTune")

# load the package (assuming Seurat is already installed)
library(scAutoTune)
library(Seurat)
```

# Main Functions
autotune_find_pcs()
- calculates variance explained per PC across multiple HVG settings
- identifies optimum number of PCs via elbow plot

autotune_find_features_resolution()
- 

# Contact
Shi Hai Sun
shihaisun.scott@gmail.com
