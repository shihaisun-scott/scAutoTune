#' Explore PCs on the sample PBMC data
#'
#' Uses the bundled pbmc_clustered.rds sample dataset to compute
#' variance-explained curves for PCs at different numbers of features,
#' fits a smooth curve, and suggests an "elbow" (optimal PCs).
#'
#' @param feature_steps Integer vector of feature counts to try.
#'   Defaults to seq(500, 4000, by = 500).
#' @param max_pcs Maximum number of PCs to compute.
#' @return A list with:
#'   \itemize{
#'     \item pc_df: data frame of variance per PC and feature count
#'     \item elbow_df: data frame for the chosen feature set with fitted curve
#'     \item suggested_pcs: integer suggested PC cut
#'     \item plot: ggplot object of the overlap / elbow plot
#'   }
#' @export
run_pc_elbow_sample <- function(
    feature_steps = seq(500, 6000, by = 500),
    max_pcs = 40
) {
  # 1) Locate and load the sample data
  sample_path <- system.file("sample_data", "pbmc_clustered.rds", package = "scAutoTune")
  if (sample_path == "") {
    stop("Sample data not found. Make sure pbmc_clustered.rds is in inst/sample_data.")
  }
  obj <- readRDS(sample_path)

  # 2) Subset cluster 0 and remove specific labels
  cluster0_cells <- Seurat::WhichCells(obj, idents = "0")
  cluster0_cells <- cluster0_cells[obj@meta.data[cluster0_cells, "cell_label"] != "cytotoxic_t"]
  cluster0_cells <- cluster0_cells[obj@meta.data[cluster0_cells, "cell_label"] != "cd14_monocytes"]

  expr <- Seurat::GetAssayData(obj[, cluster0_cells], slot = "counts")
  cell_labels_cluster0 <- obj$cell_label[cluster0_cells]

  # 3) Compute variance explained for each feature step
  feature_steps <- sort(unique(feature_steps))
  pc_variances <- list()

  for (f in feature_steps) {
    message("Processing ", f, " features")

    obj_tmp <- Seurat::CreateSeuratObject(expr)
    obj_tmp$cell_label <- cell_labels_cluster0
    obj_tmp <- Seurat::SCTransform(
      obj_tmp,
      verbose = FALSE,
      return.only.var.genes = FALSE,
      variable.features.n = f
    )
    obj_tmp <- Seurat::RunPCA(obj_tmp, npcs = max_pcs, verbose = FALSE)

    sdev <- obj_tmp[["pca"]]@stdev
    var_explained <- sdev^2 / sum(sdev^2)

    pc_variances[[as.character(f)]] <- var_explained
  }

  # 4) Build a long data.frame
  pc_df <- do.call(rbind, lapply(names(pc_variances), function(f) {
    data.frame(
      PC = seq_along(pc_variances[[f]]),
      Variance = pc_variances[[f]],
      Features = as.numeric(f)
    )
  }))

  # 5) Choose one feature count curve to analyze for elbow (e.g. max features)
  chosen_features <- max(pc_df$Features)
  elbow_df <- pc_df[pc_df$Features == chosen_features, ]

  # Fit a smooth curve (GAM) to Variance ~ PC
  # (you already load mgcv in your script)
  # Fit smoother curve
  gam_fit <- mgcv::gam(Variance ~ s(PC, k = 4, bs = "cs"), data = elbow_df, method = "REML")
  elbow_df$Fitted <- as.numeric(predict(gam_fit, newdata = elbow_df))

  # Derivative and flattening detection
  d1 <- diff(elbow_df$Fitted)
  start_idx <- 3
  search_d1 <- d1[start_idx:length(d1)]
  ref <- median(abs(search_d1[1:3]))
  gain_threshold <- 0.05 * ref
  idx_flat <- which(abs(search_d1) < gain_threshold)[1] + start_idx
  suggested_pcs <- if (!is.na(idx_flat)) idx_flat else max_pcs

  # 6) Overlapping elbow plot for all feature sets
  pc_df$Features <- factor(pc_df$Features)

  overlap_plot <- ggplot2::ggplot(pc_df, ggplot2::aes(
    x = PC,
    y = Variance,
    group = Features,
    colour = Features
  )) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_vline(
      xintercept = suggested_pcs,
      linetype = "dashed",
      colour = "red"
    ) +
    ggplot2::labs(
      title = paste0(
        "Overlapping Elbow Plots (chosen features: ", chosen_features,
        "; suggested PCs: ", suggested_pcs, ")"
      ),
      x = "Principal Component",
      y = "Proportion of Variance Explained"
    ) +
    ggplot2::theme_minimal()

  # Optional: overlay the fitted curve for the chosen feature count
  overlap_plot <- overlap_plot +
    ggplot2::geom_line(
      data = elbow_df,
      ggplot2::aes(x = PC, y = Fitted),
      inherit.aes = FALSE,
      linetype = "dotted",
      linewidth = 0.8,
      colour = "black"
    )

  # 7) Return results; user can manually choose PCs based on suggested_pcs
  list(
    pc_df = pc_df,
    elbow_df = elbow_df,
    suggested_pcs = suggested_pcs,
    plot = overlap_plot
  )
}
