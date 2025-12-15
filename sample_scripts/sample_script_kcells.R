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


# load k cells seurat object
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
sweep$chosen_params


# Downstream analysis
pcs <-sweep$chosen_params$n_pcs
nfeatures <- sweep$chosen_params$nfeatures
res <- sweep$chosen_params$resolution

pcs <- 5
nfeatures <- 600
res <- 0.2

obj_lgnK[["RNA"]] <- split(obj_lgnK[["RNA"]], f = obj_lgnK$donor)
obj_lgnK <- SCTransform(object = obj_lgnK, assay = "RNA", variable.features.n = nfeatures, return.only.var.genes = FALSE)
obj_lgnK <- RunPCA(object = obj_lgnK, npcs = 30)
obj_lgnK <- IntegrateLayers(object = obj_lgnK, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony")
obj_lgnK <- RunUMAP(obj_lgnK, reduction = "harmony", dims = 1:pcs)
obj_lgnK <- FindNeighbors(obj_lgnK, reduction = "harmony", dims = 1:pcs)
obj_lgnK <- FindClusters(obj_lgnK, resolution = res)
obj_lgnK <- PrepSCTFindMarkers(obj_lgnK,  assay = "SCT")
obj_lgnK.markers <- FindAllMarkers(
  obj_lgnK, assay = "SCT",
  only.pos = TRUE)
significant_markers <- obj_lgnK.markers %>%
  filter(p_val_adj < 0.05 )%>%
  group_by(cluster)%>%
  slice_max(order_by = avg_log2FC, n = 10)
hm <- DoHeatmap(obj_lgnK, features = significant_markers$gene)
print(hm)
dp <- DotPlot(obj_lgnK,
              features = significant_markers$gene,
              group.by = "seurat_clusters",
              assay = "SCT")
print(dp)

# check donor integration
umap_donor <- DimPlot(obj_lgnK, group.by = "donor", alpha = 0.5)
print(umap_donor)

# save plots
save_plots <- plot_grid(hm, dp, umap_donor, ncol = 3)
savename <- paste0("plots_kcells.pdf")
pdf(savename, height = 5, width = 20)
print(save_plots)
dev.off()
