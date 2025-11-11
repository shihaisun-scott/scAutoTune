#' Estimate the optimal number of PCs for a Seurat object
#'
#' Loads a Seurat object from an .rds file, computes variance explained curves
#' for varying numbers of features, and uses a LOESS smoother plus
#' maximum-distance ("elbow") method to suggest an optimal number of PCs.
#'
#' @param seurat_path Path to a Seurat .rds file.
#' @param feature_steps Integer vector of feature counts to try.
#'   Defaults to seq(500, 6000, by = 500).
#' @param max_pcs Maximum number of PCs to compute. Default is 40.
#' @param span Smoothing parameter for LOESS fit. Smaller = tighter fit.
#' @return A list containing:
#'   \itemize{
#'     \item pc_df: variance explained data frame for all feature counts
#'     \item elbow_df: variance + fitted curve for chosen feature count
#'     \item suggested_pcs: suggested number of PCs (+1 adjustment)
#'     \item plot: ggplot object with overlapping elbow curves
#'   }
#' @export
autotune_find_pcs <- function(
    seurat_path,
    feature_steps = seq(500, 6000, by = 500),
    max_pcs = 40,
    span = 0.2
) {
  # --- 1. Load Seurat object ---
  if (!file.exists(seurat_path)) stop("File not found: ", seurat_path)
  message("Loading Seurat object: ", seurat_path)
  obj <- readRDS(seurat_path)
  if (!inherits(obj, "Seurat")) stop("Input file must contain a Seurat object.")

  # --- 2. Compute variance explained for each feature step ---
  expr <- Seurat::GetAssayData(obj, slot = "counts")
  pc_variances <- list()

  for (f in sort(unique(feature_steps))) {
    message("Processing ", f, " features")

    obj_tmp <- Seurat::CreateSeuratObject(expr)
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

  pc_df <- do.call(rbind, lapply(names(pc_variances), function(f) {
    data.frame(
      PC = seq_along(pc_variances[[f]]),
      Variance = pc_variances[[f]],
      Features = as.numeric(f)
    )
  }))

  # --- 3. Choose largest feature count for elbow detection ---
  chosen_features <- max(pc_df$Features)
  elbow_df <- pc_df[pc_df$Features == chosen_features, ]

  # --- 4. Fit loess smoother and detect elbow (+1) ---
  lo_fit <- stats::loess(Variance ~ PC, data = elbow_df, span = span)
  elbow_df$Fitted <- as.numeric(predict(lo_fit, newdata = elbow_df))

  x <- elbow_df$PC
  y <- elbow_df$Fitted
  line_vec <- c(x[length(x)] - x[1], y[length(y)] - y[1])
  line_len <- sqrt(sum(line_vec^2))
  distances <- abs((y[length(y)] - y[1]) * x -
                     (x[length(x)] - x[1]) * y +
                     x[length(x)] * y[1] - y[length(y)] * x[1]) / line_len
  suggested_pcs <- x[which.max(distances)] + 1
  suggested_pcs <- min(suggested_pcs, max_pcs)

  # --- 5. Build overlap plot ---
  pc_df$Features <- factor(pc_df$Features)
  overlap_plot <- ggplot2::ggplot(pc_df, ggplot2::aes(
    x = PC, y = Variance, colour = Features, group = Features
  )) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_vline(xintercept = suggested_pcs,
                        linetype = "dashed", colour = "red") +
    ggplot2::labs(
      title = paste0("Elbow plots (suggested PCs: ", suggested_pcs, ")"),
      x = "Principal Component", y = "Proportion of Variance Explained"
    ) +
    ggplot2::theme_minimal()

  # --- 6. Return results ---
  list(
    pc_df = pc_df,
    elbow_df = elbow_df,
    suggested_pcs = suggested_pcs,
    plot = overlap_plot
  )
}
