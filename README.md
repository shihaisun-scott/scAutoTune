
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scAutoTune

scAutoTune helps find optimal parameters (PCs, HVGs, clustering
resolution) for single-cell transcriptomics data.

## Installation

# install.packages(“devtools”)

devtools::install_github(“yourusername/scAutoTune”)

## Example

library(scAutoTune) autotune_sc(“example.h5ad”)
