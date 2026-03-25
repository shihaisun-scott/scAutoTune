# Extracting and visualizing the results  
obj_lgnK[["RNA"]] <- split(obj_lgnK[["RNA"]],
                            f = obj_lgnK$donor) 
obj_lgnK <- SCTransform(object = obj_lgnK,
                         assay = "RNA",
                         variable.features.n = nfeatures,
                         return.only.var.genes = FALSE) 
obj_lgnK <- RunPCA(object = obj_lgnK, npcs = 30) 
obj_lgnK <- IntegrateLayers(object = obj_lgnK,
                             method = HarmonyIntegration,
                             orig.reduction = "pca",
                             new.reduction = "harmony")
obj_lgnK <- RunUMAP(obj_lgnK,
                     reduction = "harmony",
                     dims = 1:pcs) 
obj_lgnK <- FindNeighbors(obj_lgnK,
                           reduction = "harmony",
                           dims = 1:pcs) 
obj_lgnK <- FindClusters(obj_lgnK, resolution = res) 
obj_lgnK <- PrepSCTFindMarkers(obj_lgnK,  assay = "SCT")
obj_lgnK.markers <- FindAllMarkers(obj_lgnK, 
                                    assay = "SCT",
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
