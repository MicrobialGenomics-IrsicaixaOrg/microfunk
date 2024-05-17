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
#' The function takes a SummarizedExperiment object containing normalized gene/
#' pathway abundance values as the assay and metadata as the column data, as
#' well as several parameter options to run the analysis, using a zero-inflated
#' Gaussian mixture model as default (`fitZig()`). To use a zero-inflated
#' Log-Normal mixture model for each feature (`fitFeatureModel()`), the argument
#' `model`should be modified (`model = "fitFeatureModel"`).It returns the
#' results of the analysis as a data frame (tibble).
#'
#' @param se A SummarizedExperiment object containing normalized gene/pathway
#'   abundance values.
#' @param variable A character string specifying the name of the variable of
#'   interest to include in the model, matching the desired column name of the
#'   sample metadata.
#' @param mod document...
#' @param coef document...
#' @param B document ..
#' @param szero document
#' @param spos document
#'
#' @return A data frame (tibble) containing the results of the MetagenomeSeq
#'   analysis.
#' @export
#' @autoglobal
#' @tests
#' # Read HUMAnN3
#' data <- read_humann(
#'   file_path = system.file("extdata", "All_genefam_cpm_kegg.tsv", package = "microfunk"),
#'   metadata = system.file("extdata", "ex_meta.csv", package = "microfunk")
#'  )
#'
#' # Invalid model
#' expect_error(run_metagenomeseq(se = data, model = "fit", variable = "ARM"))
#'
#' # MetagenomeSeq Analysis (FitZig)
#' da_fitzig <- run_metagenomeseq(se = data, variable = "ARM")
#'
#' # Test name of features returned
#' n_zig <- c("K07488", "K09805", "K00863", "K02977", "K11905", "K00263", "K00737",
#'            "K11904", "K09961", "K22902")
#' f_zig <- as.vector(da_fitzig$function_id)
#' testthat::expect_equal(n_zig, f_zig)
#'
#' # Test p-values
#' da_fitzig %>%
#'  dplyr::pull(adjPvalues) %>%
#'  mean() %>%
#'  round(3) %>%
#' testthat::expect_equal(0.027)
#'
#' # MetagenomeSeq Analysis (FitFeatureModel)
#' da_fitfeature <- run_metagenomeseq(se = data,
#'                    model = "fitFeatureModel", variable = "ARM")
#'
#' # Test name of features returned
#' n_feature <- c("K03300", "K16951", "K00135", "K07488", "K06015", "K06196",
#'                "K01908", "K19130", "K01060", "K09858")
#' f_feature <- as.vector(da_fitfeature$function_id)
#' testthat::expect_equal(n_feature, f_feature)
#'
#' # Test p-values
#' da_fitfeature %>%
#'  dplyr::pull(adjPvalues) %>%
#'  mean() %>%
#'  round(3) %>%
#' testthat::expect_equal(0.258)
#'
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

.is_normalized <- function(se) {
  method <-
    SummarizedExperiment::assays(se)$humann %>%
    tibble::as_tibble(rownames = "function_id") %>%
    .extract_norm_method()

  ifelse(method == "cpm", TRUE, FALSE)
}
