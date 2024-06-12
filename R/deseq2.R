#' Run DESeq Analysis
#'
#' This function performs a DESeq analysis
#' (https://bioconductor.org/packages/release/bioc/html/DESeq2.html) on count
#' data, which tests for differential expression based on a model using the
#' Negative Binomial (Gamma-Poisson) Distribution. The function takes gene
#' family/pathway abundance data for multiple samples and associated metadata
#' contained in a SummarizedExperiment object, runs the analysis according to
#' several parameter options and applies shrinkage to log2-fold changes in order
#' to improve estimate precision. Finally, it returns a tibble of results
#' indicating significant changes between the specified experimental conditions.
#'
#' @param se A SummarizedExperiment object containing normalized abundance
#'   values as the assay and metadata as the column data.
#' @param factor A character string specifying the factor/categorical variable
#'   to use for differential expression analysis. This should be a column name
#'   in the colData of the SummarizedExperiment object.
#' @param design A formula. the formula expresses how the counts for each gene
#'   depend on the variables in colData. Many R formula are valid, including
#'   designs with multiple variables, e.g., ~ group + condition, and designs
#'   with interactions, e.g., ~ genotype + treatment + genotype:treatment. See
#'   results for a variety of designs and how to extract results tables.
#' @param test Either "Wald" or "LRT", which will then use either Wald
#'   significance tests (defined by nbinomWaldTest), or the likelihood ratio
#'   test on the difference in deviance between a full and reduced model formula
#'   (defined by nbinomLRT).
#' @param fitType Either "parametric", "local", "mean", or "glmGamPoi" for the
#'   type of fitting of dispersions to the mean intensity.
#' @param betaPrior Whether or not to put a zero-mean normal prior on the
#'   non-intercept coefficients. The default is set to FALSE, and shrunken LFCs
#'   are obtained afterwards using lfcShrink.
#' @param type A character string specifying the type of shrinkage to apply to
#'   log2-fold changes. Options include "normal", "apeglm", "ashr", etc. Default
#'   is "normal".
#' @param max_significance A numeric value specifying the maximum adjusted
#'   p-value (padj) threshold for considering a feature as significant. Default
#'   is 0.05.
#' @param log2FC A numeric value specifying the minimum absolute log2-fold
#'   change threshold for considering a feature as significant. Default is 0.
#' @param quiet Whether to print messages description
#'
#' @return A tibble containing log2-fold changes, p-values, p-adjusted values
#'   and significance indicators for each feature.
#' @export
#' @autoglobal
#' @examples
#' # Def data paths
#' metadata <-
#'   system.file("extdata", "ex_meta.csv", package = "microfunk")
#' file_path <-
#'   system.file("extdata", "All_genefam_rpk_kegg.tsv", package = "microfunk")
#'
#' # Read HUMAnN3 & DESeq2 Analysis
#' da_result <-
#'   read_humann(file_path, metadata) %>%
#'   run_deseq2(factor = "ARM")
#'
#' # Test number of significant associations
#' da_result %>%
#'   dplyr::filter(signif == TRUE) %>%
#'   nrow() %>%
#'   testthat::expect_equal(41)
#'
#' # Test P-values
#' da_result %>%
#'   dplyr::pull(pvalue) %>%
#'   mean() %>%
#'   round(3) %>%
#'   testthat::expect_equal(0.555)
#'
#' # Test P-adjusted values
#' da_result %>%
#'   dplyr::pull(padj) %>%
#'   mean(na.rm = TRUE) %>%
#'   round(3) %>%
#'   testthat::expect_equal(0.857)
#'
#' @examples
#' # Read HUMAnN3 & DESeq2 Analysis
#' da_result <- read_humann(
#'   file_path = system.file("extdata", "All_genefam_rpk_kegg.tsv", package = "microfunk"),
#'   metadata = system.file("extdata", "ex_meta.csv", package = "microfunk")
#' ) %>% run_deseq2(factor = "ARM")
#'
#' da_result
run_deseq2 <- function(se,
                       factor,
                       design = NULL,
                       test = "Wald",
                       fitType = "local",
                       betaPrior = FALSE,
                       type = "ashr",
                       max_significance = 0.05,
                       log2FC = 0,
                       quiet = TRUE) {


  if (type %in% c("ashr", "apeglm")) {
    tested <- try(find.package(type), silent = TRUE)
    if (methods::is(tested, "try-error")) {
      ins <- paste0('BiocManger::install("', type,'")') %>% cli::col_blue()
      cli::cli_abort(c(
        "i" = glue::glue('{cli::col_blue(type)} packege is not installed'),
        "*" = glue::glue('Start a clean R session and then run: {ins}')
      ))
    }
  }

  cts <-
    SummarizedExperiment::assays(se)$humann %>%
    tibble::as_tibble(rownames = "function_id") %>%
    dplyr::filter(!stringr::str_detect(function_id, "[|]")) %>%
    tibble::column_to_rownames(var = "function_id") %>%
    as.matrix() %>%
    round()

  coldata <-
    SummarizedExperiment::colData(se) %>%
    as.data.frame() %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(factor), as.factor))

  if (is.null(design)) { design <- stringr::str_c("~", factor) }
  if (!inherits(design, "formula")) { design <- stats::as.formula(design)}

  da_result <- DESeq2::DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design
  ) %>%
    DESeq2::DESeq(
      test = test,
      fitType = fitType,
      betaPrior = betaPrior,
      quiet = quiet
    )

  cont <- coldata %>% dplyr::pull(factor) %>% levels()
  lfcs <- DESeq2::lfcShrink(
    dds = da_result,
    contrast = c(factor, cont),
    type = type,
    quiet = quiet
  ) %>%
    tibble::as_tibble(rownames = "function_id") %>%
    dplyr::mutate(
      effect = log2FoldChange,
      signif = ifelse(
        padj < max_significance & abs(log2FoldChange) >= log2FC,
        TRUE,
        FALSE
      )
    )

  lfcs
}


