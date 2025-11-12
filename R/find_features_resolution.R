#' Sweep for optimal number of features and clustering resolution
#'
#' Performs a grid search over variable feature counts and clustering resolutions
#' to identify settings that maximize silhouette and modularity scores.
#'
#' @param seurat_path Path to a Seurat .rds file.
#' @param n_pcs Number of principal components to use. Default = 10.
#' @param nfeatures_range Numeric vector of feature counts to test (e.g. seq(500, 6000, by = 500)).
#' @param resolutions Numeric vector of clustering resolutions to test (e.g. seq(0.02, 1.2, by = 0.02)).
#' @param output_csv Optional path to save the results as a CSV file.
#' @return A data frame with columns: nfeatures, resolution, silhouette, modularity.
#' @export
autotune_find_features_resolution <- function(
    seurat_path,
    n_pcs = 10,
    nfeatures_range = seq(500, 6000, by = 500),
    resolutions = seq(0.02, 1.2, by = 0.02),
    output_csv = NULL
) {
  # --- 1. Load Seurat object ---
  if (!file.exists(seurat_path))
    stop("File not found: ", seurat_path)
  message("Loading Seurat object: ", seurat_path)
  seurat_obj <- readRDS(seurat_path)
  if (!inherits(seurat_obj, "Seurat"))
    stop("Input file must be a Seurat object")

  # --- 2. Use full object expression matrix ---
  expr <- Seurat::GetAssayData(seurat_obj, slot = "counts")

  # --- 3. Initialize results list ---
  all_results <- list()

  # --- 4. Parameter sweep ---
  for (nfeatures in nfeatures_range) {
    message(sprintf("Processing nfeatures: %d", nfeatures))

    # Preprocess for this feature setting
    obj <- Seurat::CreateSeuratObject(expr)
    obj <- Seurat::SCTransform(
      obj, verbose = FALSE,
      return.only.var.genes = FALSE,
      variable.features.n = nfeatures
    )
    obj <- Seurat::RunPCA(obj, npcs = max(n_pcs, 30), verbose = FALSE)
    obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = 1:n_pcs, verbose = FALSE)
    obj <- Seurat::FindNeighbors(obj, reduction = "pca", dims = 1:n_pcs, verbose = FALSE)

    for (res in resolutions) {
      message(sprintf("  Resolution = %.2f", res))
      obj <- Seurat::FindClusters(obj, resolution = res, verbose = FALSE)
      clusters <- as.integer(Seurat::Idents(obj))

      if (length(unique(clusters)) > 1) {
        # Silhouette score
        pca_mat <- Seurat::Embeddings(obj, "pca")
        dmat <- as.matrix(dist(pca_mat[, 1:n_pcs]))
        sil <- cluster::silhouette(clusters, dmat)
        sil_score <- mean(sil[, 3], na.rm = TRUE)

        # Modularity from SNN graph
        snn_graph <- obj@graphs$SCT_snn
        ig <- igraph::graph_from_adjacency_matrix(
          as.matrix(snn_graph),
          mode = "undirected",
          weighted = TRUE
        )
        modularity_score <- igraph::modularity(ig, clusters)
      } else {
        sil_score <- NA
        modularity_score <- NA
      }

      all_results[[length(all_results) + 1]] <- data.frame(
        nfeatures = nfeatures,
        n_pcs = n_pcs,
        resolution = res,
        silhouette = sil_score,
        modularity = modularity_score
      )
    }
  }

  # --- 5. Combine results ---
  results_df <- do.call(rbind, all_results)

  # --- 6. Save results if requested ---
  if (!is.null(output_csv)) {
    utils::write.csv(results_df, output_csv, row.names = FALSE)
    message("Results saved to: ", output_csv)
  }

  # --- 7. Return results ---
  return(results_df)
}
