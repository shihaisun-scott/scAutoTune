#' HVG × resolution sweep with optional Harmony batch correction + GAM optimization.
#'
#' @param seurat_obj A Seurat object already loaded in R.
#' @param n_pcs Number of PCs to use in reduction.
#' @param nfeatures_range Range of HVGs to sweep.
#' @param resolutions Clustering resolutions to sweep.
#' @param k_val GAM smoothing df.
#' @param batch_var OPTIONAL metadata column (e.g., "donor") for Harmony integration.
#' @param output_csv Optional CSV file path.
#' @param output_pdf Optional PDF path for heatmap + UMAP.
#'
#' @return A list containing:
#'   \itemize{
#'     \item results_df
#'     \item chosen_params
#'     \item heatmap
#'     \item umap_cluster_plot
#'     \item umap_label_plot
#'   }
#'
#' @export
autotune_find_features_resolution <- function(
    seurat_obj,
    n_pcs = 10,
    nfeatures_range = seq(500, 6000, by = 500),
    resolutions = seq(0.02, 1.2, by = 0.02),
    k_val = 10,
    batch_var = NULL,
    output_csv = NULL,
    output_pdf = NULL
) {

  if (!inherits(seurat_obj, "Seurat"))
    stop("Input must be a Seurat object")

  expr <- Seurat::GetAssayData(seurat_obj, slot = "counts")
  all_results <- list()

  # -------------------------------------------------------------
  # Sweep over nfeatures × resolution
  # -------------------------------------------------------------
  for (nf in sort(unique(as.integer(nfeatures_range)))) {

    message("Sweeping HVGs = ", nf)

    obj <- Seurat::CreateSeuratObject(expr)

    # ---------------------- HARMONY BRANCH ----------------------
    if (!is.null(batch_var)) {
      obj[[batch_var]] <- seurat_obj[[batch_var]][colnames(obj)]
      obj[["RNA"]] <- split(obj[["RNA"]], f = obj[[batch_var]])
    }

    obj <- Seurat::SCTransform(
      obj, variable.features.n = nf,
      return.only.var.genes = FALSE, verbose = FALSE
    )
    obj <- Seurat::RunPCA(obj, npcs = max(n_pcs, 30), verbose = FALSE)

    if (!is.null(batch_var)) {
      obj <- Seurat::IntegrateLayers(
        obj,
        method = HarmonyIntegration,
        orig.reduction = "pca",
        new.reduction = "harmony",
        assay = "SCT",
        verbose = FALSE
      )
      obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
      reduction_used <- "harmony"
    } else {
      reduction_used <- "pca"
    }

    obj <- Seurat::RunUMAP(obj, reduction = reduction_used, dims = 1:n_pcs, verbose = FALSE)
    obj <- Seurat::FindNeighbors(obj, reduction = reduction_used, dims = 1:n_pcs, verbose = FALSE)

    emb <- Seurat::Embeddings(obj, reduction_used)[, 1:n_pcs]
    dmat <- as.matrix(dist(emb))

    # ---------------- RESOLUTION LOOP ----------------
    for (res in sort(unique(as.numeric(resolutions)))) {

      obj <- Seurat::FindClusters(obj, resolution = res, verbose = FALSE)
      clusters <- as.integer(Idents(obj))

      if (length(unique(clusters)) > 1) {
        sil <- cluster::silhouette(clusters, dmat)
        sil_score <- mean(sil[, 3])
        snn <- obj@graphs$SCT_snn
        ig <- igraph::graph_from_adjacency_matrix(
          as.matrix(snn), weighted = TRUE, mode = "undirected"
        )
        modularity_score <- igraph::modularity(ig, clusters)
      } else {
        sil_score <- NA
        modularity_score <- NA
      }

      all_results[[length(all_results) + 1]] <- data.frame(
        nfeatures = nf,
        n_pcs = n_pcs,
        resolution = res,
        silhouette = sil_score,
        modularity = modularity_score
      )
    }
  }

  results_df <- do.call(rbind, all_results)

  if (!is.null(output_csv))
    write.csv(results_df, output_csv, row.names = FALSE)

  # -------------------------------------------------------------
  # Fit GAM over silhouette surface
  # -------------------------------------------------------------
  df <- results_df[!is.na(results_df$silhouette), ]

  gam_model <- mgcv::gam(
    silhouette ~ s(nfeatures, k = k_val) + s(resolution, k = k_val),
    data = df
  )

  df$pred <- predict(gam_model, df)
  best <- df[which.max(df$pred), ]

  chosen_params <- list(
    nfeatures  = best$nfeatures,
    n_pcs      = best$n_pcs,
    resolution = best$resolution
  )

  # -------------------------------------------------------------
  # Heatmap
  # -------------------------------------------------------------
  grid <- expand.grid(
    nfeatures = nfeatures_range,
    resolution = resolutions
  )
  grid$pred <- predict(gam_model, grid)

  heatmap <- ggplot2::ggplot(
    grid,
    aes(factor(nfeatures), factor(resolution), fill = pred)
  ) +
    geom_tile() +
    geom_point(
      data = best,
      aes(factor(nfeatures), factor(resolution)),
      color = "red", size = 3
    ) +
    scale_fill_viridis_c() +
    theme_minimal() +
    labs(
      title = "GAM-smoothed silhouette heatmap",
      x = "HVGs",
      y = "Resolution",
      fill = "Silhouette"
    )

  # -------------------------------------------------------------
  # Rebuild FINAL object at optimal parameters
  # -------------------------------------------------------------
  obj_final <- Seurat::CreateSeuratObject(expr)

  if (!is.null(batch_var)) {
    obj_final[[batch_var]] <- seurat_obj[[batch_var]][colnames(obj_final)]
    obj_final[["RNA"]] <- split(obj_final[["RNA"]], f = obj_final[[batch_var]])
  }

  obj_final <- Seurat::SCTransform(
    obj_final, variable.features.n = chosen_params$nfeatures,
    return.only.var.genes = FALSE, verbose = FALSE
  )
  obj_final <- Seurat::RunPCA(obj_final, npcs = max(chosen_params$n_pcs, 30), verbose = FALSE)

  if (!is.null(batch_var)) {
    obj_final <- Seurat::IntegrateLayers(
      obj_final,
      method = HarmonyIntegration,
      orig.reduction = "pca",
      new.reduction = "harmony",
      assay = "SCT",
      verbose = FALSE
    )
    obj_final[["RNA"]] <- JoinLayers(obj_final[["RNA"]])
    red_final <- "harmony"
  } else {
    red_final <- "pca"
  }

  obj_final <- Seurat::RunUMAP(obj_final, reduction = red_final, dims = 1:chosen_params$n_pcs)
  obj_final <- Seurat::FindNeighbors(obj_final, reduction = red_final, dims = 1:chosen_params$n_pcs)
  obj_final <- Seurat::FindClusters(obj_final, resolution = chosen_params$resolution)
  obj_final$cluster_id <- Idents(obj_final)

  umap_cluster_plot <- DimPlot(obj_final, group.by = "cluster_id") +
    ggplot2::ggtitle(
      sprintf("UMAP (PC=%d, HVG=%d, Res=%.2f)",
              chosen_params$n_pcs,
              chosen_params$nfeatures,
              chosen_params$resolution)
    ) +
    theme_minimal()

  umap_label_plot <- NULL
  if ("cell_label" %in% colnames(seurat_obj@meta.data)) {
    obj_final$cell_label <- seurat_obj$cell_label[colnames(obj_final)]
    umap_label_plot <- DimPlot(obj_final, group.by = "cell_label") +
      theme_minimal() + labs(title = "UMAP — true labels")
  }

  # PDF
  if (!is.null(output_pdf)) {
    cowplot::ggsave2(
      output_pdf,
      plot = cowplot::plot_grid(
        heatmap, umap_cluster_plot,
        if (!is.null(umap_label_plot)) umap_label_plot else ggplot(),
        ncol = 3, rel_widths = c(1.3, 1, 1)
      ),
      width = 16, height = 7
    )
  }

  return(list(
    results_df        = results_df,
    chosen_params     = chosen_params,
    heatmap           = heatmap,
    umap_cluster_plot = umap_cluster_plot,
    umap_label_plot   = umap_label_plot
  ))
}
