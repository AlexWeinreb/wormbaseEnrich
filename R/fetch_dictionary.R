#' Fetch a dictionary
#'
#' If `analysis` isn't specified, fetches the tissue dictionary.
#'
#'
#' @param analysis the type of analysis (tissue, phenotype, Gene Ontology)
#'
#' @section Output:
#' A dataframe containing the dictionary of interest
#'
#' @export
fetch_dictionary <- function(analysis = c("tissue", "phenotype", "go")) {
  analysis <- tolower(analysis)
  analysis <- match.arg(analysis)
  tea$fetch_dictionary(
    analysis = analysis
  )
}


