install.packages(c(
  "Seurat",
  "R.utils",
  "ggplot2",
  "dplyr",
  "cluster",
  "mgcv",
  "igraph",
  "cowplot",
  "gridExtra",
  "Matrix",
  "devtools",
	"harmony",
  "BiocManager,"
  "Remotes"))

BiocManager::install("glmGamPoi")
devtools::install_github("satijalab/seurat-wrappers")
devtools::install_github("shihaisun-scott/scAutoTune")

# check if installed
library(Seurat)
library(scAutoTune)
