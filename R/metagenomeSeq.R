#' Run MetagenomeSeq Analysis
#'
#' MetagenomeSeq is designed to determine features (such as Operational
#' Taxanomic Unit (OTU), species, genes, etc.) that are differentially abundant
#' between two or more groups of multiple samples
#' (https://www.cbcb.umd.edu/software/metagenomeSeq). It address the effects of
#' both normalization and under-sampling of microbial communities on disease
#' association detection and the testing of feature correlations. For this
#' function, features consist on the abundance of gene families or metabolic
#' pathways of multiple samples, obtained using HUMAnN3
#' (https://huttenhower.sph.harvard.edu/humann).
#'
#' The function takes a SummarizedExperiment object containing gene/ pathway
#' abundance values as the assay and metadata as the column data, as well as
#' several parameter options to run the analysis, using a zero-inflated
#' Log-Normal mixture model for each feature (`fitFeatureModel()`).It converts
#' the SummarizedExperiment object to a metagenomeSeq object, checks if the
#' count data is normalized, and if not, calculates the normalization statistics
#' to apply cumulative sum scaling normalization. It constructs a model matrix,
#' fits the model to the data and returns the results of the analysis as a data
#' frame (tibble).
#'
#' @param se A SummarizedExperiment object containing normalized gene/pathway
#'   abundance values.
#' @param variable A character string specifying the name of the variable of
#'   interest to include in the model, matching the desired column name of the
#'   sample metadata.
#' @param mod An optional model matrix for the count distribution. If NULL,
#'   model matrix will be created using the provided variable.
#' @param coef An integer specifying the coefficient of interest to grab log
#'   fold-changes. Default is 2.
#' @param B An integer specifying the number of bootstraps to perform for the
#'   fit. Default is 1.
#' @param szero A logical value indicating whether to shrink zero component
#'   parameters. Default is FALSE.
#' @param spos A logical value indicating whether to shrink positive component
#'   parameters. Default is FALSE.
#'
#' @return A data frame (tibble) containing the results of the MetagenomeSeq
#'   analysis.
#' @export
#' @autoglobal
#' @tests
#' # Read HUMAnN3 & MetagenomeSeq Analysis
#' da_result <- read_humann(
#'  file_path = system.file("extdata", "All_genefam_rpk_kegg.tsv", package = "microfunk"),
#'  metadata = system.file("extdata", "ex_meta.csv", package = "microfunk") ) %>%
#'  run_metagenomeseq(variable = "ARM")
#'
#'  # Test name of features returned
#'  f <- c("K03300", "K00863", "K19130", "K16951", "K07488", "K00135", "K05522",
#'       "K06015", "K07776", "K06196")
#'  res <- as.vector(da_result$function_id)
#'  testthat::expect_equal(f, res)
#'
#'  # Test p-values
#'  da_result %>%
#'    dplyr::pull(adjPvalues) %>%
#'    mean() %>%
#'    round(3) %>%
#'    testthat::expect_equal(0.147)
#'
#' @examples
#' # Read HUMAnN3 & MetagenomeSeq Analysis
#' da_result <- read_humann(
#'   file_path = system.file("extdata", "All_genefam_cpm_kegg.tsv", package = "microfunk"),
#'   metadata = system.file("extdata", "ex_meta.csv", package = "microfunk")
#'  ) %>% run_metagenomeseq(variable = "ARM")
#'
#' da_result
run_metagenomeseq <- function(se,
                              variable,
                              mod = NULL,
                              coef = 2,
                              B = 1,
                              szero = FALSE,
                              spos = FALSE) {

  mr_exp <- .create_MRexp(se = se)

  p <- 0
  if (!.is_normalized(se)) { p <- metagenomeSeq::cumNormStatFast(mr_exp) }
  mr_exp <- metagenomeSeq::cumNorm(mr_exp, p = p)

  if (is.null(mod)) {
    mod <-
      stats::as.formula(glue::glue("~ 1 + { variable }")) %>%
      stats::model.matrix(data = Biobase::pData(mr_exp))
  }

  metagenomeSeq::fitFeatureModel(
    obj = mr_exp,
    mod = mod,
    coef = coef,
    B = B,
    szero = szero,
    spos = spos
  ) %>%
    metagenomeSeq::MRfulltable() %>%
    tibble::as_tibble(rownames = "function_id")
}

#' Create an MRexperiment Object
#'
#' This function creates an MRexperiment object from a normalized
#' SummarizedExperiment object containing gene family or pathway abundance
#' values and associated metadata.
#'
#' @param se A SummarizedExperiment object containing gene family or pathway
#'   abundance data.
#'
#' @return An MRexpriment object containing the gene family or pathway abundance
#'   values as `counts` and their associated metadata as `pheno`.
#' @autoglobal
#' @keywords internal
#' @noRd
.create_MRexp <- function(se) {

  counts <-
    SummarizedExperiment::assays(se)$humann %>%
    tibble::as_tibble(rownames = "function_id") %>%
    dplyr::filter(!stringr::str_detect(function_id, "[|]")) %>%
    tibble::column_to_rownames(var = "function_id") %>%
    data.frame(check.names = FALSE)

  pheno <- se %>%
    SummarizedExperiment::colData() %>%
    data.frame() %>%
    Biobase::AnnotatedDataFrame()

  mr_exp <- metagenomeSeq::newMRexperiment(counts, phenoData = pheno)
}

#' Check if SummarizedExperiment object is normalized
#'
#' This function checks whether the assay contained in a SummarizedExperiment
#' object is normalized to CPM (counts per million)
#'
#' @param se A SummarizedExperiment object containing the assay data to be
#'   checked for normalization.
#'
#' @return A logical value (TRUE/FALSE) indicating whether the assay is
#'   normalized using the CPM method.
#' @autoglobal
#' @keywords internal
#' @noRd
.is_normalized <- function(se) {
  method <-
    SummarizedExperiment::assays(se)$humann %>%
    tibble::as_tibble(rownames = "function_id") %>%
    .extract_norm_method()

  ifelse(method == "cpm", TRUE, FALSE)
}
