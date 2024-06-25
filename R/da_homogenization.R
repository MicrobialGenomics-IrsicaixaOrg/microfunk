#' Homogenize Differential Abundance Results
#'
#' This internal function takes differential abundance analysis results in a
#' tibble format from either DESeq2 or MaAsLin2 and homogenizes them into a
#' consistent format. The function first determines the type of input
#' (DESeq2/MaAsLin2) based on the column names and then calls for another
#' internal function to create a new tibble with homogenized column names. If
#' the input type is unrecognized, an error is raised.
#'
#' @param da_result A tibble containing the differential abundance results.
#'   The data frame should have columns corresponding to either MaAsLin2 or
#'   DESeq2 output formats.
#'
#' @return A homogenized tibble with a consistent format regardless of the
#'   original DA result format.
#' @autoglobal
#' @tests
#' # Def file paths
#' file_path = system.file("extdata", "reduced_genefam_cpm_kegg.tsv", package = "microfunk")
#' metadata = system.file("extdata", "reduced_meta.csv", package = "microfunk")
#'
#' # Read HUMAnN3 & DA Analysis
#' da_deseq <- read_humann(file_path, metadata) %>%
#'  run_deseq2(factor = "ARM")
#' da_maaslin <- read_humann(file_path, metadata) %>%
#'  run_maaslin2(fixed_effects = "ARM")
#'
#'  # Test column names
#'  n <- c("feature", "metadata", "da_method", "coefficient", "log2FC", "stderr",
#'         "lfc_stderr", "p_value", "q_value", "padj_value", "signif")
#'
#'  testthat::expect_equal(n, colnames(da_deseq))
#'  testthat::expect_equal(n, colnames(da_maaslin))
#'
#'  # Check unknown input type
#'  da_unknown <- tibble::tibble(x = 1:10, y = 1:10)
#'  testthat::expect_error(.da_homogenization(da_unknown))
.da_homogenization <- function(da_result){

  # Classify input type
  maaslin <- c("feature", "metadata", "value", "coef", "stderr", "pval", "name",
               "qval", "N", "N.not.zero", "signif")
  deseq <- c("function_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue",
             "padj", "signif", "metadata")

  input_type <-  dplyr::case_when(
    all(names(da_result) %in% maaslin) ~ "MaAsLin",
    all(names(da_result) %in% deseq) ~ "DESeq",
    TRUE ~ "Unknown"
  )

  # Create homogenized table
  result <- if (input_type == "MaAsLin") {
    .maaslin_homogenize(tbl = da_result)
  } else if (input_type == "DESeq") {
    .deseq_homogenize(tbl = da_result)
  } else {
    cli::cli_abort(c(
      "x" = "Unknown input type.",
      "i" = "Please provide DA results of DESeq2 or MaAsLin2 analysis."))
  }

  result

}


#' Homogenize MaAsLin2 Differential Abundance Results
#'
#' This internal function homogenizes MaAsLin2 differential abundance (DA)
#' results into a consistent format with specific column names. The function
#' processes the input tibble and ensures that the output tibble has a
#' consistent structure with columns relevant for downstream analysis. It fills
#' certain columns with placeholder values as they are not applicable to
#' MaAsLin2 results.
#'
#' @param tbl A tibble containing the MaAsLin2 DA analysis results.
#'
#' @return A tibble with a homogenized format.
#' @autoglobal
.maaslin_homogenize <- function(tbl){

  n <- nrow(tbl)
  h_tbl <- tibble::tibble(
    feature = tbl %>% dplyr::pull(feature),
    metadata = tbl %>% dplyr::pull(metadata),
    da_method = rep("maaslin2", n),
    coefficient = tbl %>% dplyr::pull(coef),
    log2FC = rep("-", n),
    stderr = tbl %>% dplyr::pull(stderr),
    lfc_stderr = rep("-", n),
    p_value = tbl %>% dplyr::pull(pval),
    q_value = tbl %>% dplyr::pull(qval),
    padj_value = rep("-", n),
    signif = tbl %>% dplyr::pull(signif)
  )
}



#' Homogenize DESeq2 Differential Abundance Results
#'
#' This internal function homogenizes DESeq2 differential abundance (DA)
#' results into a consistent format with specific column names. The function
#' processes the input tibble and ensures that the output tibble has a
#' consistent structure with columns relevant for downstream analysis. It fills
#' certain columns with placeholder values as they are not applicable to
#' DESeq2 results.
#'
#' @param tbl A tibble containing the DESeq2 DA analysis results.
#'
#' @return A tibble with a homogenized format.
#' @autoglobal
.deseq_homogenize <- function(tbl){

  n <- nrow(tbl)
  h_tbl <- tibble::tibble(
    feature = tbl %>% dplyr::pull(function_id),
    metadata = tbl %>% dplyr::pull(metadata),
    da_method = rep("deseq2", n),
    coefficient = rep("-", n),
    log2FC = tbl %>%  dplyr::pull(log2FoldChange),
    stderr = rep("-", n),
    lfc_stderr = tbl %>% dplyr::pull(lfcSE),
    p_value = tbl %>% dplyr::pull(pvalue),
    q_value = rep("-", n),
    padj_value = tbl %>% dplyr::pull(padj),
    signif = tbl %>% dplyr::pull(signif)
  )
}
