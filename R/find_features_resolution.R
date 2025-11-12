#' Sweep HVGs (nfeatures) and clustering resolution; pick optimum via GAM + plot
#'
#' Loads a Seurat .rds, sweeps over `nfeatures` and `resolution` using the
#' given `n_pcs`, computes silhouette & modularity, smooths the surface with a
#' GAM, selects optimal parameters, and returns plots + results.
#'
#' @param seurat_path Path to a Seurat .rds file.
#' @param n_pcs Integer, number of PCs to use in neighbors/UMAP. Default 10.
#' @param nfeatures_range Integer vector of HVG counts to evaluate.
#' @param resolutions Numeric vector of resolutions to evaluate.
#' @param k_val Integer, GAM degrees of freedom (per smoother). Default 10.
#' @param output_csv Optional path to write the sweep table (.csv).
#' @param output_pdf Optional path to write a combined figure (.pdf).
#' @return A list with:
#'   \itemize{
#'     \item results_df: sweep table (nfeatures, n_pcs, resolution, silhouette, modularity)
#'     \item selections_df: summary table with all selection methods
#'     \item chosen_params: list (nfeatures, n_pcs, resolution)
#'     \item heatmap: ggplot heatmap (GAM-smoothed)
#'     \item umap_cluster_plot: ggplot UMAP colored by clusters
#'     \item umap_label_plot: ggplot UMAP by `cell_label` if present (else NULL)
#'   }
#' @export
autotune_find_features_resolution <- function(
    seurat_path,
    n_pcs = 10,
    nfeatures_range = seq(500, 6000, by = 500),
    resolutions     = seq(0.02, 1.2,  by = 0.02),
    k_val = 10,
    output_csv = NULL,
    output_pdf = NULL
) {
  # --- 1) Load Seurat object ---
  if (!file.exists(seurat_path)) stop("File not found: ", seurat_path)
  message("Loading Seurat object: ", seurat_path)
  obj0 <- readRDS(seurat_path)
  if (!inherits(obj0, "Seurat")) stop("Input file must be a Seurat object")
  expr <- Seurat::GetAssayData(obj0, slot = "counts")

  # --- 2) Sweep over nfeatures × resolution ---
  all_results <- list()
  for (nf in sort(unique(as.integer(nfeatures_range)))) {
    message(sprintf("Processing nfeatures: %d", nf))

    obj <- Seurat::CreateSeuratObject(expr)
    obj <- Seurat::SCTransform(
      obj, verbose = FALSE,
      return.only.var.genes = FALSE,
      variable.features.n = nf
    )
    obj <- Seurat::RunPCA(obj, npcs = max(n_pcs, 30), verbose = FALSE)
    obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = 1:n_pcs, verbose = FALSE)
    obj <- Seurat::FindNeighbors(obj, reduction = "pca", dims = 1:n_pcs, verbose = FALSE)
    pca_mat <- Seurat::Embeddings(obj, "pca")
    dmat <- as.matrix(stats::dist(pca_mat[, 1:n_pcs, drop = FALSE]))

    for (res in sort(unique(as.numeric(resolutions)))) {
      obj <- Seurat::FindClusters(obj, resolution = res, verbose = FALSE)
      clusters <- as.integer(Seurat::Idents(obj))

      if (length(unique(clusters)) > 1) {
        sil <- cluster::silhouette(clusters, dmat)
        sil_score <- mean(sil[, 3], na.rm = TRUE)
        snn_graph <- obj@graphs$SCT_snn
        ig <- igraph::graph_from_adjacency_matrix(
          as.matrix(snn_graph), mode = "undirected", weighted = TRUE
        )
        modularity_score <- igraph::modularity(ig, clusters)
      } else {
        sil_score <- NA_real_
        modularity_score <- NA_real_
      }

      all_results[[length(all_results) + 1]] <- data.frame(
        nfeatures  = nf,
        n_pcs      = n_pcs,
        resolution = res,
        silhouette = sil_score,
        modularity = modularity_score
      )
    }
  }
  results_df <- do.call(rbind, all_results)

  if (!is.null(output_csv)) {
    utils::write.csv(results_df, output_csv, row.names = FALSE)
    message("Results saved to: ", output_csv)
  }

  # --- 3) Fit GAM and find optimum (smoothed) ---
  df <- results_df[!is.na(results_df$silhouette), , drop = FALSE]
  df_fit <- df[, c("nfeatures","resolution","silhouette")]
  gam_model <- mgcv::gam(
    silhouette ~ s(nfeatures, k = k_val) + s(resolution, k = k_val),
    data = df_fit
  )
  df$predicted <- as.numeric(stats::predict(gam_model, newdata = df))
  best_gam <- dplyr::slice_max(df, predicted, n = 1, with_ties = FALSE) |>
    dplyr::mutate(method = "smoothed_gam")

  chosen_params <- list(
    nfeatures  = as.integer(best_gam$nfeatures[[1]]),
    n_pcs      = as.integer(best_gam$n_pcs[[1]]),
    resolution = as.numeric(best_gam$resolution[[1]]),
    method     = "smoothed_gam"
  )

  # --- 4) Build smoothed heatmap ---
  grid <- expand.grid(
    nfeatures  = sort(unique(results_df$nfeatures)),
    resolution = sort(unique(results_df$resolution))
  )
  grid$silhouette_pred <- as.numeric(stats::predict(gam_model, newdata = grid))
  df_hm <- dplyr::left_join(grid, results_df, by = c("nfeatures","resolution"))
  df_hm$silhouette_combined <- ifelse(
    is.na(df_hm$silhouette), df_hm$silhouette_pred, df_hm$silhouette
  )

  heatmap <- ggplot2::ggplot(
    df_hm,
    ggplot2::aes(x = factor(nfeatures), y = factor(round(resolution, 2)))
  ) +
    ggplot2::geom_tile(ggplot2::aes(fill = silhouette_combined)) +
    ggplot2::geom_point(
      data = best_gam,
      ggplot2::aes(x = factor(nfeatures), y = factor(round(resolution, 2))),
      color = "red", size = 3, shape = 21, fill = "white", stroke = 1
    ) +
    ggplot2::scale_fill_viridis_c(option = "D", na.value = "grey90") +
    ggplot2::labs(
      title = "GAM-smoothed silhouette heatmap (red = optimum)",
      x = "Number of features",
      y = "Resolution (rounded)",
      fill = "Silhouette"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )

  # --- 5) Rebuild object with chosen parameters & plot UMAP ---
  obj_final <- Seurat::CreateSeuratObject(expr)
  obj_final <- Seurat::SCTransform(
    obj_final, verbose = FALSE,
    return.only.var.genes = FALSE,
    variable.features.n = chosen_params$nfeatures
  )
  obj_final <- Seurat::RunPCA(obj_final, npcs = max(chosen_params$n_pcs, 30), verbose = FALSE)
  obj_final <- Seurat::RunUMAP(obj_final, dims = 1:chosen_params$n_pcs, verbose = FALSE)
  obj_final <- Seurat::FindNeighbors(obj_final, dims = 1:chosen_params$n_pcs, verbose = FALSE)
  obj_final <- Seurat::FindClusters(obj_final, resolution = chosen_params$resolution, verbose = FALSE)
  obj_final$cluster_id <- Seurat::Idents(obj_final)

  umap_cluster_plot <- Seurat::DimPlot(obj_final, group.by = "cluster_id") +
    ggplot2::ggtitle(sprintf(
      "Clusters — GAM optimum (PCs=%d, HVGs=%d, res=%.2f)",
      chosen_params$n_pcs, chosen_params$nfeatures, chosen_params$resolution
    )) +
    ggplot2::theme_minimal()

  umap_label_plot <- NULL
  if ("cell_label" %in% colnames(obj0@meta.data)) {
    if (all(colnames(obj_final) %in% colnames(obj0))) {
      obj_final$cell_label <- obj0@meta.data[colnames(obj_final), "cell_label", drop = TRUE]
      umap_label_plot <- Seurat::DimPlot(obj_final, group.by = "cell_label") +
        ggplot2::ggtitle("Labels (from input)") +
        ggplot2::theme_minimal()
    }
  }

  # --- 6) Optionally save a combined figure ---
  if (!is.null(output_pdf)) {
    cowplot::ggsave2(
      output_pdf,
      plot = cowplot::plot_grid(
        heatmap,
        umap_cluster_plot,
        if (is.null(umap_label_plot)) ggplot2::ggplot() + ggplot2::theme_void() else umap_label_plot,
        ncol = 3, rel_widths = c(1.5, 2, 2)
      ),
      width = 16, height = 8, units = "in"
    )
    message("Figure saved to: ", output_pdf)
  }

  list(
    results_df        = results_df,
    chosen_params     = chosen_params,
    heatmap           = heatmap,
    umap_cluster_plot = umap_cluster_plot,
    umap_label_plot   = umap_label_plot
  )
}
