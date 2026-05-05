# Default Seurat workflow

library(Seurat)
library(ggplot2)

default_pcs <- 30

# pbmc
obj_pbmc <- readRDS("data/pbmc_cluster0.rds")

obj_default <- obj_pbmc

obj_default <- NormalizeData(obj_default)
obj_default <- FindVariableFeatures(obj_default)   # default: 2000 features
obj_default <- ScaleData(obj_default)
obj_default <- RunPCA(obj_default)
obj_default <- FindNeighbors(obj_default, dims = 1:default_pcs)
obj_default <- FindClusters(obj_default)           # default resolution = 0.8
obj_default <- RunUMAP(obj_default, dims = 1:default_pcs)

p_default <- DimPlot(
  obj_default,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) +
  ggtitle("Default Seurat parameters")

p_default
