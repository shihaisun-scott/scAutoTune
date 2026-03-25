
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

obj_lgnK <- readRDS("data/K_cells_seurat.rds")

# Running scAutoTune (1) – Optimal number of PCs
pc_lgnK <- autotune_find_pcs(obj_lgnK,
                              batch_var = "donor",
                              feature_steps = seq(500, 6000, by = 500),
                              max_pcs = 20)

# Running scAutoTune (2) – Parameter sweep
sweep <- autotune_find_features_resolution(
  obj_lgnK, 
  batch_var = "donor",
  n_pcs = 5, 
  nfeatures_range = seq(500, 3000, by = 50), 
  resolutions = seq(0.2, 1.0, by = 0.02),
  k_val = 10, 
  output_csv = "k_sweep.csv",
  output_pdf = "k_heatmap_umap.pdf")

# check features
sweep$chosen_params
pcs <-sweep$chosen_params$n_pcs
nfeatures <- sweep$chosen_params$nfeatures
res <- sweep$chosen_params$resolution




