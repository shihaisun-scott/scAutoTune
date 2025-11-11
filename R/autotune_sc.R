#' Auto-tune single-cell parameters (stub)
#'
#' Minimal placeholder so the package builds. Later, this will read
#' an h5ad/h5seurat and search for optimal n_pcs, n_hvg, and resolution.
#'
#' @param h5_path Character path to the input file (not used yet).
#' @param batch_key Optional metadata column for batch correction (not used yet).
#' @return A named list with three numbers: n_pcs_opt, n_hvg_opt, resolution_opt.
#' @examples
#' res <- autotune_sc(tempfile())
#' str(res)
#' @export
autotune_sc <- function(h5_path, batch_key = NULL) {
  # placeholder return values so users can install and test the package
  list(
    n_pcs_opt = 30L,
    n_hvg_opt = 2000L,
    resolution_opt = 0.8
  )
}
