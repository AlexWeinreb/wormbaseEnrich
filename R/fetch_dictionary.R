#' Fetch a dictionary
#'
#' Download a dictionary file from Wormbase
#'
#'
#' @param analysis the type of analysis (tissue, phenotype, go)
#'
#' @returns A data.frame containing the dictionary of interest, one row per gene, one column per term
#'
#' @export
fetch_dictionary <- function(analysis = c("tissue", "phenotype", "go")) {
  analysis <- tolower(analysis)
  if(length(analysis) != 1L){
    stop("Only provide one type of analysis")
  }
  analysis <- match.arg(analysis)



  url_base = 'http://caltech.wormbase.org/TissueEnrichmentAnalysis/'

  filename <- switch (analysis,
    "tissue" = "anatomy_dict.csv",
    "phenotype" = "phenotype_dict.csv",
    "go" = "go_dict.csv"
  )

  utils::read.csv(paste0(url_base, filename), check.names = FALSE)

}


