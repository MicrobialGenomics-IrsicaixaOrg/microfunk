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
#'
#' @return A tibble containing log2-fold changes, p-values, p-adjusted values
#'   and significance indicators for each feature.
#' @export
#' @autoglobal
#' @tests
#' # Def data paths
#' metadata <- system.file("extdata", "ex_meta.csv", package = "microfunk")
#' file_path <-
#'   system.file("extdata", "All_genefam_rpk_kegg.tsv", package = "microfunk")
#'
#' # Read HUMAnN3 & DESeq2 Analysis
#' da_result <-
#'   read_humann(file_path, metadata) %>%
#'   run_deseq2(factor = "ARM",
#'              test = "Wald",
#'              fitType = "local")
#'
#' # Test number of significant associations
#' da_result %>%
#'   dplyr::filter(signif == TRUE) %>%
#'   nrow() %>%
#'   testthat::expect_equal(40)
#'
#' # Test P-values
#' da_result %>%
#'   dplyr::pull(pvalue) %>%
#'   mean() %>%
#'   round(3) %>%
#'   testthat::expect_equal(0.557)
#'
#' # Test P-adjusted values
#' da_result %>%
#'   dplyr::pull(padj) %>%
#'   mean(na.rm = TRUE) %>%
#'   round(3) %>%
#'   testthat::expect_equal(0.858)
#'
#' @examples
#' # Read HUMAnN3 & DESeq2 Analysis
#' da_result <- read_humann(
#'   file_path = system.file("extdata", "All_genefam_rpk_kegg.tsv", package = "microfunk"),
#'   metadata = system.file("extdata", "ex_meta.csv", package = "microfunk")
#' ) %>% run_deseq2(factor = "ARM",
#'                  test = "Wald",
#'                  fitType = "local")
#' da_result
run_deseq2 <- function(se,
                       factor,
                       test,
                       fitType,
                       betaPrior = FALSE,
                       type = "normal",
                       max_significance = 0.05,
                       log2FC = 0,
                       quiet = TRUE){

  cts <- SummarizedExperiment::assays(se)$humann %>%
    tibble::as_tibble(rownames = "function_id") %>%
    dplyr::filter(!stringr::str_detect(function_id, "[|]")) %>%
    tibble::column_to_rownames(var = "function_id") %>%
    as.matrix() %>%
    round()

  coldata <- SummarizedExperiment::colData(se) %>% as.data.frame()

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                        colData = coldata,
                                        design = stats::as.formula(
                                          stringr::str_c("~", factor)))
  da_result <- DESeq2::DESeq(dds)

  cont <- coldata %>% dplyr::select(all_of(factor)) %>% unique() %>% as.vector()

  lfcs <- DESeq2::lfcShrink(
    dds = da_result,
    contrast = c(factor, cont[[1]]),
    type = type,
    quiet = TRUE
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


