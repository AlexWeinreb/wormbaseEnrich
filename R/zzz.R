# global reference to tea (will be initialized in .onLoad)
tea <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to tea
  # if(! reticulate::virtualenv_exists("r-tea")){
  #   if(! reticulate::py_module_available("tissue_enrichment_analysis")){
  #     warning("No installation of the tea package detected, run the `install_tea()` function.")
  #     return()
  #   }
  #   tea <<- reticulate::import("tissue_enrichment_analysis", delay_load = TRUE)
  # } else{
  #   reticulate::use_virtualenv("r-tea", required = FALSE)
  #   tea <<- reticulate::import("tissue_enrichment_analysis", delay_load = TRUE)
  # }
  tea <<- reticulate::import("tissue_enrichment_analysis", delay_load = TRUE)
}
