
# Load Libraries and Seurat object
library(Seurat)
library(ggplot2) 
library(dplyr) 
library(cluster) 
library(mgcv) 
library(igraph) 
library(cowplot) 
library(gridExtra) 
library(Matrix) 
library(harmony) 
library(SeuratWrappers) 
library(scAutoTune)

set.seed(67)

obj_pbmc <- readRDS("data/pbmc_cluster0.rds")
obj_lgnK <- readRDS("data/K_cells_seurat.rds")

# Running scAutoTune (1) – Optimal number of PCs
steps <- seq(500, 6000, by = 500)
max_pc <- 30
pc_pbmc <- autotune_find_pcs(obj_pbmc, feature_steps = steps, max_pcs = max_pc)
pc_pbmc$suggested_pcs
pc_pbmc$plot

# Running scAutoTune (2) – Parameter sweep
pcs <- 5
features <- seq(500, 3000, by = 50)
resolution <- seq(0.2, 1.0, by = 0.02)
k_val <- 10
filename_sweep <- "results/pbmc_sweep.csv"
filename_umap <- "results/pbmc_heatmap_umap.pdf"

results <- autotune_find_features_resolution(
  obj_pbmc, 
  n_pcs = 5, 
  nfeatures_range = features, 
  resolutions = resolution,
  k_val = k_val, 
  output_csv = filename_sweep,
  output_pdf = filename_umap)


