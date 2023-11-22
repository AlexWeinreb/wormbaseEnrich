



#' Install the tissue enrichment analysis Python package
#'
#' @param ... Additional arguments to be passed to `reticulate::py_install()`
#' @param envname Name for a Python virtual environment in which to perform install
#'
#' @export
#'
install_tea <- function(..., envname = "r-tea"){
  reticulate::py_install("tissue_enrichment_analysis", envname = envname, ...)
}
