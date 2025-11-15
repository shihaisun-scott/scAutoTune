#' Determine optimal number of PCs using LOESS elbow detection
#' with optional Harmony batch correction.
#'
#' @param seurat_obj A Seurat object already loaded into R.
#' @param feature_steps Vector of HVG numbers to test.
#' @param max_pcs Maximum number of PCs to compute.
#' @param span LOESS smoothing span.
#' @param batch_var OPTIONAL metadata column name for batch correction
#'        (e.g., "donor"). If provided, Harmony integration is used.
#'
#' @return A list containing:
#'   \itemize{
#'     \item pc_df        – all variance explained curves
#'     \item elbow_df     – curve for chosen HVG count
#'     \item suggested_pcs – LOESS-based elbow estimate
#'     \item plot         – ggplot object
#'   }
#'
#' @export
autotune_find_pcs <- function(
    seurat_obj,
    feature_steps = seq(500, 6000, by = 500),
    max_pcs = 40,
    span = 0.25,
    batch_var = NULL
) {
  if (!inherits(seurat_obj, "Seurat"))
    stop("Input must be a Seurat object")

  expr <- Seurat::GetAssayData(seurat_obj, slot = "counts")
  pc_variances <- list()
  feature_steps <- sort(unique(feature_steps))

  for (f in feature_steps) {
    message("Processing HVGs = ", f)

    obj <- Seurat::CreateSeuratObject(expr)

    # ---------------------- BATCH CORRECTION ----------------------
    if (!is.null(batch_var)) {
      obj[[batch_var]] <- seurat_obj[[batch_var]][colnames(obj)]
      obj[["RNA"]] <- split(obj[["RNA"]], f = obj[[batch_var]])
    }

    obj <- Seurat::SCTransform(
      obj, variable.features.n = f,
      return.only.var.genes = FALSE, verbose = FALSE
    )
    obj <- Seurat::RunPCA(obj, npcs = max_pcs, verbose = FALSE)

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
      pcs <- obj[["harmony"]]@stdev
    } else {
      pcs <- obj[["pca"]]@stdev
    }

    pc_variances[[as.character(f)]] <- pcs^2 / sum(pcs^2)
  }

  # ------------------------ DATAFRAME ------------------------
  pc_df <- do.call(rbind,
                   lapply(names(pc_variances), function(f){
                     data.frame(
                       PC = seq_along(pc_variances[[f]]),
                       Variance = pc_variances[[f]],
                       Features = as.numeric(f)
                     )
                   }))

  chosen_features <- max(pc_df$Features)
  elbow_df <- pc_df[pc_df$Features == chosen_features, ]

  lo <- stats::loess(Variance ~ PC, data = elbow_df, span = span)
  elbow_df$Fitted <- predict(lo)

  d <- diff(elbow_df$Fitted)
  suggested_pcs <- which.min(d) + 1

  p <- ggplot2::ggplot(pc_df, aes(PC, Variance, color = Features)) +
    geom_line(alpha = 0.6) +
    geom_line(data = elbow_df, aes(PC, Fitted),
              color = "black", linetype = "dotted") +
    geom_vline(xintercept = suggested_pcs, color = "red") +
    theme_minimal() +
    labs(title = sprintf("Suggested PCs = %d", suggested_pcs))

  return(list(
    pc_df = pc_df,
    elbow_df = elbow_df,
    suggested_pcs = suggested_pcs,
    plot = p
  ))
}
