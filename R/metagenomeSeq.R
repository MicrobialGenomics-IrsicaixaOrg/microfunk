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
#' @param model A character string specifying the model to use. Options are
#'   "fitZig" (default) or "fitFeatureModel".
#' @param zeroMod The model to account for the change in the number of features
#'   observed as a linear effect of the depth of coverage. Default is NULL.
#' @param useMixedModel Estimate the correlation between duplicate features or
#'   replicates using duplicateCorrelation.
#' @param useCSSoffset Whether to include the default scaling parameters in the
#'   model. Default is FALSE.
#' @param max_significance The p-value threshold for significance. Default is
#'   0.05.
#'
#' @return A data frame (tibble) containing the results of the MetagenomeSeq
#'   analysis.
#' @export
#' @autoglobal
#' @tests
#' # Read HUMAnN3 & MetagenomeSeq Analysis (FitZig)
#' da_fitzig <- read_humann(
#'   file_path = system.file("extdata", "All_genefam_cpm_kegg.tsv", package = "microfunk"),
#'   metadata = system.file("extdata", "ex_meta.csv", package = "microfunk")
#'  ) %>% run_metagenomeseq(variable = "ARM")
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
#' # Read HUMAnN3 & MetagenomeSeq Analysis (FitFeatureModel)
#' da_fitfeature <- read_humann(
#'   file_path = system.file("extdata", "All_genefam_cpm_kegg.tsv", package = "microfunk"),
#'   metadata = system.file("extdata", "ex_meta.csv", package = "microfunk")
#'  ) %>% run_metagenomeseq(model = "fitFeatureModel", variable = "ARM")
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
                              variable = NULL,
                              model = "fitZig",
                              zeroMod = NULL,
                              useMixedModel = FALSE,
                              useCSSoffset = TRUE,
                              max_significance = 0.05){
  if(model == "fitZig"){
    .run_fitzig(se = se,
                variable = variable,
                zeroMod = zeroMod,
                useMixedModel = useMixedModel,
                useCSSoffset = useCSSoffset,
                max_significance = max_significance)
  } else if(model == "fitFeatureModel"){
    .run_fitfeaturemodel(se = se,
                         variable = variable,
                         max_significance = max_significance)
  } else{cli::cli_abort(c("x" = "Invalid model requested.",
                     "i" = "Available models are 'fitZig' and 'fitFeatureModel'."))
  }
}


#' Run MetagenomeSeq Fit Zig Model
#'
#' This function fits a zero-inflated Gaussian model using MetagenomeSeq's
#' `fitZig()` function. It takes as input a SummarizedExperiment object
#' containing gene/pathway abundance values and associated metadata. The
#' function performs data normalization using cumulative sum scaling (CSS) and
#' fits a zero-inflated Gaussian model to identify features that are
#' differentially abundant according to the specified variable of interest and a
#' normalization factor obtained from the data.
#'
#' @param se A SummarizedExperiment object containing normalized gene/pathway
#'   abundance values.
#' @param variable A character string specifying the name of the variable of
#'   interest to include in the model, matching the desired column name of the
#'   sample metadata.
#' @param zeroMod The model to account for the change in the number of features
#'   observed as a linear effect of the depth of coverage.
#' @param useCSSoffset Whether to include the default scaling parameters in the
#'   model.
#' @param max_significance The p-value threshold for significance.
#'
#' @return A data frame (tibble) containing the results of the feature-level
#'   model fitting, including an extra column indicating significance.
#' @autoglobal
#' @keywords internal
#' @noRd
.run_fitzig <- function(se,
                        variable,
                        zeroMod,
                        useMixedModel,
                        useCSSoffset,
                        max_significance){

  mr_exp <- .create_MRexp(se = se)

  p <- metagenomeSeq::cumNormStatFast(mr_exp)
  mr_obj <- metagenomeSeq::cumNorm(mr_exp, p = p)
  v <- Biobase::pData(mr_obj)[variable][[1]]

  norm_factor <- metagenomeSeq::normFactors(mr_obj)
  norm_factor <- log2(norm_factor / stats::median(norm_factor) + 1)

  mod_matrix <- stats::model.matrix(~1 + v + norm_factor)

  fit <- metagenomeSeq::fitZig(obj = mr_obj,
                               mod = mod_matrix,
                               useCSSoffset = useCSSoffset,
                               zeroMod = zeroMod,
                               useMixedModel = useMixedModel,
                               control = metagenomeSeq::zigControl(
                                 verbose = FALSE,
                                 maxit = 100,
                                 dfMethod = "modified"
                               )) %>%
    metagenomeSeq::MRfulltable() %>%
    tibble::as_tibble(rownames = "function_id") %>%
    dplyr::mutate(signif = ifelse(adjPvalues < max_significance, TRUE, FALSE))

}


#' Run MetagenomeSeq Fit Feature Model
#'
#' This function fits a feature-level model using MetagenomeSeq's
#' `fitFeatureModel()` function. It takes as input a SummarizedExperiment object
#' containing gene/pathway abundance values and associated metadata. The
#' function performs data normalization using cumulative sum scaling (CSS) and
#' fits a zero-inflated Log-Normal mixture model to identify features that are
#' differentially abundant according to the specified variable of interest.
#'
#' @param se A SummarizedExperiment object containing normalized gene/pathway
#'   abundance values.
#' @param variable A character string specifying the name of the variable of
#'   interest to include in the model, matching the desired column name of the
#'   sample metadata.
#' @param max_significance The p-value threshold for significance.
#'
#' @return A data frame (tibble) containing the results of the feature-level
#'   model fitting, including an extra column indicating significance.
#' @autoglobal
#' @keywords internal
#' @noRd
.run_fitfeaturemodel <- function(se,
                                 variable,
                                 max_significance){

  mr_exp <- .create_MRexp(se = se)

  p <- metagenomeSeq::cumNormStatFast(mr_exp)
  mr_obj <- metagenomeSeq::cumNorm(mr_exp, p = p)

  v <- Biobase::pData(mr_obj)[variable][[1]]
  pd <- Biobase::pData(mr_obj)

  mod_matrix <- stats::model.matrix(~1 + v, data = pd)

  fit <- metagenomeSeq::fitFeatureModel(obj = mr_obj,
                                        mod = mod_matrix) %>%
    metagenomeSeq::MRfulltable() %>%
    tibble::as_tibble(rownames = "function_id") %>%
    dplyr::mutate(signif = ifelse(adjPvalues < max_significance, TRUE, FALSE))

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
.create_MRexp <- function(se){

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

