# Extracting and visualizing the results  
results$results_df 
print(results$heatmap)
print(results$umap_cluster_plot)
results$chosen_params

library(Seurat)
pcs <-results$chosen_params$n_pcs
nfeatures <- results$chosen_params$nfeatures
res <- results$chosen_params$resolution

obj_pbmc <- SCTransform(object = obj_pbmc, 
                        assay = "RNA",
                        variable.features.n = nfeatures,
                        return.only.var.genes = FALSE)
obj_pbmc <- RunPCA(object = obj_pbmc, npcs = 30)
obj_pbmc <- FindNeighbors(obj_pbmc, dims = 1:pcs)
obj_pbmc <- FindClusters(obj_pbmc, resolution = res)
obj_pbmc <- PrepSCTFindMarkers(obj_pbmc,  assay = "SCT")
obj_pbmc.markers <- FindAllMarkers(obj_pbmc,
                                   assay = "SCT", only.pos = TRUE)
significant_markers <- obj_pbmc.markers %>%
  filter(p_val_adj < 0.05 )%>%
  group_by(cluster)%>%
  slice_max(order_by = avg_log2FC, n = 10)
hm <- DoHeatmap(obj_pbmc, features = significant_markers$gene)
dp <- DotPlot(obj_pbmc,
              features = significant_markers$gene,
              group.by = "seurat_clusters",
              assay = "SCT")

