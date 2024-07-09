# library(testthat)

# compare using dictionaries downloaded on 2024-07-09
# The results were manually compared to the online tool (match, except FDR),
# future tests can simply compare to these hardcoded values


# prepare data
dict_tissue <- readr::read_csv( testthat::test_path("dict_tissue.csv"),
                                show_col_types = FALSE ) |>
  as.data.frame()

dict_pheno <- readr::read_csv( testthat::test_path("dict_pheno.csv"),
                               show_col_types = FALSE ) |>
  as.data.frame()

dict_go <- readr::read_csv( testthat::test_path("dict_go.csv"),
                            show_col_types = FALSE ) |>
  as.data.frame()

gene_list <- c("WBGene00007913", "WBGene00006890", "WBGene00007802", "WBGene00000449",
               "WBGene00009386", "WBGene00021950", "WBGene00012670", "WBGene00004301")

background_list <- c("WBGene00007913", "WBGene00006890", "WBGene00007802", "WBGene00000449",
                     "WBGene00009386", "WBGene00021950", "WBGene00012670", "WBGene00004301",
                     "WBGene00010957", "WBGene00010958", "WBGene00010959", "WBGene00010960",
                     "WBGene00000829", "WBGene00010962", "WBGene00010963", "WBGene00010964",
                     "WBGene00010965", "WBGene00010967", "WBGene00017926", "WBGene00004493",
                     "WBGene00020297", "WBGene00022046", "WBGene00004427", "WBGene00000677",
                     "WBGene00004497", "WBGene00019537", "WBGene00021420", "WBGene00022089",
                     "WBGene00021088", "WBGene00001999", "WBGene00015248", "WBGene00000611",
                     "WBGene00004419", "WBGene00017925", "WBGene00017121", "WBGene00004432",
                     "WBGene00004494", "WBGene00022235", "WBGene00004473", "WBGene00004477",
                     "WBGene00019900", "WBGene00016790", "WBGene00004471", "WBGene00019510",
                     "WBGene00017075", "WBGene00008149", "WBGene00010556", "WBGene00002344",
                     "WBGene00004474", "WBGene00011480", "WBGene00008505", "WBGene00004480",
                     "WBGene00003920", "WBGene00004448", "WBGene00004492", "WBGene00004424",
                     "WBGene00004410", "WBGene00013385", "WBGene00004430", "WBGene00004443",
                     "WBGene00001558", "WBGene00013311", "WBGene00004487", "WBGene00045146",
                     "WBGene00012768", "WBGene00013418", "WBGene00002005", "WBGene00008707",
                     "WBGene00004498", "WBGene00000263", "WBGene00022336", "WBGene00021248",
                     "WBGene00002025", "WBGene00004491", "WBGene00022170", "WBGene00006725",
                     "WBGene00004428", "WBGene00002267", "WBGene00004414", "WBGene00006537",
                     "WBGene00006959", "WBGene00004470", "WBGene00000380", "WBGene00003915",
                     "WBGene00004469", "WBGene00006519", "WBGene00000229", "WBGene00022599",
                     "WBGene00004481", "WBGene00004435", "WBGene00017166", "WBGene00004472",
                     "WBGene00018562", "WBGene00004482", "WBGene00002083", "WBGene00001168",
                     "WBGene00004483", "WBGene00004450", "WBGene00004490", "WBGene00004417",
                     "WBGene00003024", "WBGene00004420", "WBGene00004449", "WBGene00018846",
                     "WBGene00004433", "WBGene00015755", "WBGene00016250", "WBGene00018418")

test_that("Warnings as expected",{


  expect_no_condition( enrichment_analysis(gene_list, dict_tissue) )

  gene_list_duplicates <- c(gene_list, "WBGene00007802")
  expect_warning( enrichment_analysis(gene_list_duplicates, dict_tissue))

  gene_list_wrong_names <- c(gene_list, "not_a_real_name")
  expect_warning( enrichment_analysis(gene_list_wrong_names, dict_tissue))

})


test_that("Tissue enrichment as expected",{

  res <- enrichment_analysis(gene_list, dict_tissue)

  expect_equal(
    res[1,] |>
      lapply(setNames, nm = NULL),

    list(term_name = "ABplpapp",
         term_id = "WBbt:0006420",
         expected = 0.01906973,
         observed = 1,
         enrichment_fc = 52.439118,
         p_value = 0.0001733859,
         FDR = 0.01284194),
    tolerance = .00001
  )


  res_backgr <- enrichment_analysis(gene_list, dict_tissue, background_genes = background_list)

  expect_equal(
    res_backgr[1,] |>
      lapply(setNames, nm = NULL),

    list(term_name = "AVL",
         term_id = "WBbt:0003843",
         expected = 0.04971638,
         observed = 1,
         enrichment_fc = 20.11409,
         p_value = 0,
         FDR = 0),
    tolerance = .00001
  )

})




test_that("Phenotype enrichment as expected",{

  gene_list_pheno <- intersect(gene_list, dict_pheno$wbid)
  backg_list_pheno <- intersect(background_list, dict_pheno$wbid)


  # with this short list, no pheno enrichement unless background list used
  res <- enrichment_analysis(gene_list_pheno, dict_pheno)

  expect_identical(
    res |>
      lapply(setNames, nm = NULL),

    list(term_name = character(),
         term_id = character(),
         expected = double(),
         observed = integer(),
         enrichment_fc = double(),
         p_value = double(),
         FDR = double())
  )



  res_backgr <- enrichment_analysis(gene_list_pheno, dict_pheno, background_genes = backg_list_pheno)

  expect_equal(
    res_backgr[1,] |>
      lapply(setNames, nm = NULL),

    list(term_name = "male nervous system morphology variant",
         term_id = "WBPhenotype:0001358",
         expected = 0.02747792,
         observed = 1,
         enrichment_fc = 36.39286,
         p_value = 0.000182107,
         FDR = 0.0006480866),
    tolerance = .00001
  )

})



test_that("GO enrichment as expected",{


  gene_list_go <- intersect(gene_list, dict_go$wbid)
  backg_list_go <- intersect(background_list, dict_go$wbid)


  # with this short list, online need to use FDR cutoff at 0.1 instead of 0.05
  res <- enrichment_analysis(gene_list_go, dict_go)

  expect_equal(
    res[1,] |>
      lapply(setNames, nm = NULL),

    list(term_name = "steroid hydroxylase activity",
         term_id = "GO:0008395",
         expected = 0.01927075,
         observed = 1,
         enrichment_fc = 51.89211,
         p_value = 0.0001742513,
         FDR = 0.02627938),
    tolerance = .00001
  )



  res_backgr <- enrichment_analysis(gene_list_go, dict_go, background_genes = backg_list_go)

  expect_equal(
    res_backgr[1,] |>
      lapply(setNames, nm = NULL),

    list(term_name = "peptidase inhibitor activity",
         term_id = "GO:0030414",
         expected = 0.04434072,
         observed = 1,
         enrichment_fc = 22.55263,
         p_value = 0,
         FDR = 0),
    tolerance = .00001
  )

})





