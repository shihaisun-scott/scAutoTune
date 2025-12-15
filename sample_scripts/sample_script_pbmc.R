# load packages
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

# set random seed
set.seed(67)

# load seurat object
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
filename_sweep <- "pbmc_sweep.csv"
filename_umap <- "pbmc_heatmap_umap.pdf"

results <- autotune_find_features_resolution(
  obj_pbmc, 
  n_pcs = 5, 
  nfeatures_range = features, 
  resolutions = resolution,
  k_val = k_val, 
  output_csv = filename_sweep,
  output_pdf = filename_umap)

# Extracting and visualizing the results  
results$results_df 
print(results$heatmap)
print(results$umap_cluster_plot)
results$chosen_params


# Downstream analysis
pcs <-results$chosen_params$n_pcs
nfeatures <- results$chosen_params$nfeatures
res <- results$chosen_params$resolution

pcs <- 5
nfeatures <- 2750
res <- 0.3

obj_pbmc <- SCTransform(object = obj_pbmc, assay = "RNA", variable.features.n = nfeatures, return.only.var.genes = FALSE)
obj_pbmc <- RunPCA(object = obj_pbmc, npcs = 30)
obj_pbmc <- FindNeighbors(obj_pbmc, dims = 1:pcs)
obj_pbmc <- FindClusters(obj_pbmc, resolution = res)
obj_pbmc <- PrepSCTFindMarkers(obj_pbmc,  assay = "SCT")
obj_pbmc.markers <- FindAllMarkers(
  obj_pbmc, assay = "SCT",
  only.pos = TRUE)
significant_markers <- obj_pbmc.markers %>%
  filter(p_val_adj < 0.05 )%>%
  group_by(cluster)%>%
  slice_max(order_by = avg_log2FC, n = 10)
hm <- DoHeatmap(obj_pbmc, features = significant_markers$gene)
print(hm)
dp <- DotPlot(obj_pbmc,
               features = significant_markers$gene,
               group.by = "seurat_clusters",
               assay = "SCT")
print(dp)


