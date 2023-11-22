#' Execute complete enrichment analysis (hypergeometric test, BH correction).
#'
#'
#' @param gene_list a list of non-redundant Wormbase IDs.
#' @param dictionary dictionary created with `fetch_dictionary()`.
#' @param alpha significance threshold, defaults to 0.05.
#' @param filename if not NULL, file name to save results as a csv.
#'
#'
#' @returns a data.frame containing significantly enriched tissues
#'
#' @export
enrichment_analysis <- function(gene_list, tissue_df, alpha = 0.05, filename = NULL) {
  if(! is.null(filename)){
    if(!(is.character(filename) && nchar(filename) > 0)){
      stop("If filename is not NULL, it should be a valid file name.")
    }
    save <- TRUE
  } else{
    save <- FALSE
    filename <- ""
  }
  tea$enrichment_analysis(
    gene_list = gene_list,
    tissue_df = tissue_df,
    alpha = alpha,
    aname = filename,
    save = save,
    show = FALSE
  )
}

