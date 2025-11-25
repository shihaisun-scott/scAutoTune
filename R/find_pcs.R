#' Determine variance explained per PC across HVG settings
#' and detect elbows per HVG curve (no smoothing).
#'
#' Workflow (for each `feature_steps` value):
#'   1. Optionally split RNA assay by batch (`batch_var`)
#'   2. SCTransform with given number of variable features
#'   3. Run PCA
#'   4. Compute variance explained per PC
#'   5. Detect elbow on each variance curve (geometric method)
#'
#' @param seurat_obj A Seurat object already loaded into R.
#' @param feature_steps Vector of HVG numbers to test.
#' @param max_pcs Maximum number of PCs to compute.
#' @param span Ignored (kept for backwards compatibility; no smoothing used).
#' @param batch_var OPTIONAL metadata column name for batch splitting
#'        (e.g., "donor"). If provided, the RNA assay is split by this
#'        variable prior to SCTransform.
#'
#' @return A list containing:
#'   \itemize{
#'     \item pc_df         – all variance explained curves
#'                           (one row per PC per HVG)
#'     \item elbow_df      – one row per HVG, with the elbow PC and variance
#'     \item suggested_pcs – median elbow across HVG curves (not plotted)
#'     \item plot          – ggplot object with overlapping curves and
#'                           stacked elbow lines
#'   }
#'
#' @export
autotune_find_pcs <- function(
    seurat_obj,
    feature_steps = seq(500, 6000, by = 500),
    max_pcs = 40,
    batch_var = NULL
) {
  # ---------------------- INPUT CHECKS ----------------------
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (!is.null(batch_var) &&
      !batch_var %in% colnames(seurat_obj@meta.data)) {
    stop("`batch_var` ('", batch_var, "') not found in seurat_obj@meta.data")
  }

  pc_variances <- list()
  feature_steps <- sort(unique(feature_steps))

  # ---------------------- MAIN LOOP OVER HVGs ----------------------
  for (f in feature_steps) {
    message("Processing HVGs = ", f)

    obj_tmp <- seurat_obj

    # ---------------------- SPLIT BY BATCH (if requested) ----------------------
    if (!is.null(batch_var)) {
      # align batch column to current object
      batch_vec <- seurat_obj@meta.data[colnames(obj_tmp), batch_var, drop = TRUE]
      obj_tmp[[batch_var]] <- batch_vec

      # split RNA assay into per-batch layers
      # uses Seurat's split.Assay method via base::split()
      obj_tmp[["RNA"]] <- base::split(
        obj_tmp[["RNA"]],
        f = obj_tmp[[batch_var, drop = TRUE]]
      )
    }

    # ---------------------- SCTransform on (possibly split) RNA ----------------------
    obj_tmp <- Seurat::SCTransform(
      object                = obj_tmp,
      verbose               = FALSE,
      return.only.var.genes = FALSE,
      variable.features.n   = f
    )

    # ---------------------- PCA ----------------------
    obj_tmp <- Seurat::RunPCA(
      obj_tmp,
      npcs    = max_pcs,
      verbose = FALSE
    )

    # ---------------------- VARIANCE PER PC ----------------------
    sdev <- obj_tmp[["pca"]]@stdev
    sdev <- sdev[is.finite(sdev)]

    if (length(sdev) < 2) {
      warning("Skipping HVGs = ", f, " because PCA returned < 2 PCs.")
      next
    }

    var_explained <- sdev^2 / sum(sdev^2)
    pc_variances[[as.character(f)]] <- var_explained
  }

  # ------------------------ COMBINE INTO DATA FRAME ------------------------
  if (length(pc_variances) == 0) {
    stop("No valid PC variance vectors were computed.")
  }

  pc_df <- do.call(
    rbind,
    lapply(names(pc_variances), function(f) {
      v <- pc_variances[[f]]
      data.frame(
        PC       = seq_along(v),
        Variance = as.numeric(v),
        Features = as.numeric(f)
      )
    })
  )

  # Make Features a factor for nicer plotting
  pc_df$Features <- factor(pc_df$Features)

  # ------------------------ ELBOW PER HVG (NO FITTING) ------------------------
  split_by_feat <- split(pc_df, pc_df$Features)
  elbow_rows <- list()

  for (nm in names(split_by_feat)) {
    df <- split_by_feat[[nm]]

    if (nrow(df) < 3) {
      warning("Skipping Features = ", nm, ": too few PCs for elbow detection.")
      next
    }

    x <- df$PC
    y <- df$Variance

    # if all variance identical, no elbow
    if (length(unique(y)) == 1) {
      warning("Flat variance curve at Features = ", nm, "; skipping elbow.")
      next
    }

    # normalize to [0, 1]
    x_min <- min(x)
    x_max <- max(x)
    y_min <- min(y)
    y_max <- max(y)

    if (x_max == x_min || y_max == y_min) {
      warning("Degenerate scaling at Features = ", nm, "; skipping elbow.")
      next
    }

    x_n <- (x - x_min) / (x_max - x_min)
    y_n <- (y - y_min) / (y_max - y_min)

    # line from first to last point
    start <- c(x_n[1], y_n[1])
    end   <- c(x_n[length(x_n)], y_n[length(y_n)])
    v     <- end - start
    v_len <- sqrt(sum(v^2))

    if (v_len == 0) {
      warning("Flat start–end line at Features = ", nm, "; skipping elbow.")
      next
    }

    # perpendicular distance to the line
    dist <- abs(
      v[2] * (x_n - start[1]) -
        v[1] * (y_n - start[2])
    ) / v_len

    elbow_idx <- which.max(dist)
    elbow_pc  <- x[elbow_idx]

    # grab that row from the original df
    elbow_row <- df[df$PC == elbow_pc, , drop = FALSE]
    elbow_rows[[nm]] <- elbow_row
  }

  if (length(elbow_rows) == 0) {
    stop("No elbows could be determined from any HVG curve.")
  }

  elbow_df <- do.call(rbind, elbow_rows)

  # consensus PC across all curves (not plotted, but returned)
  suggested_pcs <- stats::median(elbow_df$PC)
  suggested_pcs <- as.integer(round(suggested_pcs))

  # ------------------------ PLOT: OVERLAPPING CURVES + STACKED ELBOWS ------------------------
  p <- ggplot2::ggplot(
    pc_df,
    ggplot2::aes(
      x     = PC,
      y     = Variance,
      group = Features,
      color = Features
    )
  ) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::geom_vline(
      data     = elbow_df,
      mapping  = ggplot2::aes(xintercept = PC, color = Features),
      linetype = "dashed",
      alpha    = 0.7
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Per-HVG elbow positions",
      x     = "Principal Component",
      y     = "Proportion of Variance Explained"
    )

  # ------------------------ RETURN ------------------------
  return(list(
    pc_df         = pc_df,
    elbow_df      = elbow_df,
    suggested_pcs = suggested_pcs,
    plot          = p
  ))
}
