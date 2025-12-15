cran_pkgs <- c(
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
  "BiocManager"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  BiocManager::install("glmGamPoi")
}

if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
  devtools::install_github("satijalab/seurat-wrappers")
}

if (!requireNamespace("scAutoTune", quietly = TRUE)) {
  devtools::install_github("shihaisun-scott/scAutoTune")
}
