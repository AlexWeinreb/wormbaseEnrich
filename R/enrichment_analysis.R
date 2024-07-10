#' Execute complete enrichment analysis (hypergeometric test, BH correction).
#'
#'
#' @param gene_list a list of Wormbase IDs.
#' @param dictionary dictionary created with `fetch_dictionary()`.
#' @param alpha significance threshold.
#' @param background_genes list of background genes
#'
#'
#' @returns a data.frame containing significantly enriched tissues
#'
#' @export
enrichment_analysis <- function(gene_list, dictionary, alpha = 0.05, background_genes = NULL) {



  if(! is.null(background_genes)){

    background_ids_not_found <- background_genes[! background_genes %in% dictionary$wbid]
    if(length(background_ids_not_found) > 0){
      warning(length(background_ids_not_found),
              " background gene(s) not in dictionary, will be ignored: ",
              paste0(utils::head(background_ids_not_found), collapse = ", "))

      background_genes <- setdiff(background_genes, background_ids_not_found)
    }

    dictionary <- dictionary[dictionary$wbid %in% background_genes,]
  }



  ids_not_found <- gene_list[! gene_list %in% dictionary$wbid]
  if(length(ids_not_found) > 0){

    warning(length(ids_not_found),
            " gene ID(s) not found in dictionary, will be ignored: ",
            paste0(utils::head(ids_not_found), collapse = ", "))

    gene_list <- setdiff(background_genes, ids_not_found)
  }


  if(anyDuplicated(gene_list) > 0){
    warning("Duplicate entries in gene_list were removed")
    gene_list <- unique(gene_list)
  }


  nb_terms <- ncol(dictionary) - 1L

  term_full <- names(dictionary[-1])
  term_matches <- regexec("([[:print:]]+) ((?:WBbt|WBPhenotype|GO)\\:[0-9]{7}$)", term_full)
  term_split <- regmatches(term_full, term_matches)


  res <- tibble::tibble(term_name = vapply(term_split, \(.x) .x[[2]], character(1L)),
                        term_id   = vapply(term_split, \(.x) .x[[3]], character(1L)))





  # q: nb of white (in tissue) balls (terms) picked (in list)
  # m: nb of white balls: nb of terms in tissue
  # n: nb of black balls: nb of terms not in tissue
  # k: nb of balls drawn: nb of terms in list
  #
  # Note, following the Python version, we count a gene each time it appears
  # in a tissue, so the same gene is counted several times.
  # From paper: "Although a user will input gene IDs, we test the number of occurrences
  # of a term within the gene list, so a single gene can contribute to multiple terms"


  nb_total <- sum(dictionary[-1])
  nb_listed_in_all_tissues <- sum(dictionary[dictionary$wbid %in% gene_list, -1])


  res$expected <- vapply(dictionary[-1],
                         \(ref){
                           in_tissue <- sum(ref)
                           # mean = k*m/(m+n)
                           nb_listed_in_all_tissues*in_tissue/(nb_total)

                         },
                         FUN.VALUE = double(1L))




  res$observed <- vapply(dictionary[-1],
                         \(ref){
                           intersect(gene_list, dictionary$wbid[ref == 1L]) |> length()

                         },
                         FUN.VALUE = integer(1L))

  res$enrichment_fc <- res$observed / res$expected


  res$p_value <- vapply(dictionary[-1],
                        \(ref){

                          in_ref_and_list <- intersect(gene_list, dictionary$wbid[ref == 1L]) |> length()
                          in_tissue <- sum(ref)

                          stats::phyper(q = in_ref_and_list,
                                        m = in_tissue,
                                        n = nb_total - in_tissue,
                                        k = nb_listed_in_all_tissues,
                                        lower.tail = FALSE)
                        },
                        FUN.VALUE = double(1L))



  res$FDR <- stats::p.adjust(res$p_value, method = "BH")

  res <- res[res$FDR <= alpha & res$observed > 0 & res$enrichment_fc > 1, ]

  res[order(res$enrichment_fc, decreasing = TRUE),]
}











