#' Execute complete enrichment analysis (hypergeometric test, BH correction).
#'
#'
#' @param gene_list a list of Wormbase IDs.
#' @param dictionary dictionary created with `fetch_dictionary()`.
#' @param alpha significance threshold.
#'
#'
#' @returns a data.frame containing significantly enriched tissues
#'
#' @export
enrichment_analysis <- function(gene_list, dictionary, alpha = 0.05) {


  if(!all(gene_list %in% dictionary$wbid)){
    ids_not_found <- gene_list[! gene_list %in% dictionary$wbid]
    stop(length(ids_not_found), " gene ID not found: ", paste0(head(ids_not_found), collapse = ", "))
  }

  if(anyDuplicated(gene_list) > 0){
    warning("Duplicate entries in gene_list were removed")
    gene_list <- unique(gene_list)
  }


  nb_terms <- ncol(dictionary) - 1L

  res <- tibble::tibble( term = names(dictionary[-1]) ) |>
    tidyr::separate_wider_regex(term,
                                patterns = c(term_name = "[[:print:]]+",
                                             " ",
                                             term_id = "(?:WBbt|WBPhenotype|GO)\\:[0-9]{7}$"))



  # q: nb of white (in tissue) balls picked (in list)
  # m: nb of white balls: nb of genes in tissue
  # n: nb of black balls: nb of genes not in tissue
  # k: nb of balls drawn: nb of genes in list
  #
  # Note, following the Python version, we count a gene each time it appears
  # in a tissue, so the same gene is counted several times.
  # From paper: "Although a user will input gene IDs, we test the number of occurrences
  # of a term within the gene list, so a single gene can contribute to multiple terms"

  # tot_not_in_list <- setdiff(dictionary$wbid, gene_list) |> length()
  # tot_in_ref <- sum(ref)

  nb_total <- sum(dict[-1])
  nb_listed_in_all_tissues <- sum(dict[dict$wbid %in% gene_list, -1])


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

                          phyper(q = in_ref_and_list,
                                 m = in_tissue,
                                 n = nb_total - in_tissue,
                                 k = nb_listed_in_all_tissues,
                                 lower.tail = FALSE)
                        },
                        FUN.VALUE = double(1L))



  res$q_value <- stats::p.adjust(res$p_value, method = "BH")
  res
}

res |> filter(paste(term_name, term_id) %in% res_py2$term) |> as.data.frame()


term_full <- names(dictionary[-1])
reg



stringr::str_detect(res$term[1:3],
                    pattern = "[a-zA-Z0-9]+ (?:WBbt|WBPhenotype|GO)\\:[0-9]{7}")

dict_go |>
  tidyr::pivot_longer(-wbid,
                      names_to = "term") |>
  dplyr::select(-c(wbid, value)) |>
  unique() |>
  tidyr::separate_wider_regex(term,
                              patterns = c(term2 = "[[:print:]]+", " ",
                                           term_id = "(?:WBbt|WBPhenotype|GO)\\:[0-9]{7}$"))
res |>
  tidyr::separate_wider_regex(term,
                              patterns = c(term2 = "[[:print:]]+", " ",
                                           term_id = "(?:WBbt|WBPhenotype|GO)\\:[0-9]{7}$"))













