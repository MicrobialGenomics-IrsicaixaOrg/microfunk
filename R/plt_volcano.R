#' Generate a Volcano Plot for MaAsLin2 Analysis Results
#'
#' This function generates a volcano plot based on the results of a differential
#' abundance analysis in order to visualize the differential abundance of
#' microbial gene families or pathways (features) associated with specific fixed
#' effects.
#'
#' The x-axis represents the coefficient/effect size (for MaAsLin2 DA) or
#' log2FoldChange (for DESeq2 DA) of each feature, indicating the magnitude and
#' direction of change associated with the chosen categorical variable. The
#' y-axis represents the statistical significance of the differences (negative
#' logarithm of the p-values).
#'
#' Each point on the plot corresponds to a gene family or pathway, with its
#' position determined by both its effect size and statistical significance.
#' Only features with a q-value < 0.05 (for MaAsLin2) or p-adjusted value < 0.05
#' (for DESeq2)  are deemed significant and therefore colored, in red if its
#' abundance increases or blue if its abundance decreases.
#'
#' @param da_result A tibble containing the result of the differential analysis.
#'
#' @return A volcano plot visualizing differential abundance of microbial
#'   features.
#' @export
#' @autoglobal
#' @tests
#'  # Def data paths
#'   metadata <- system.file("extdata", "ex_meta.csv", package = "microfunk")
#'   file_path <- system.file("extdata", "All_genefam_cpm_kegg.tsv", package = "microfunk")
#'
#'  # Read HUMAnN3 & MaAsLin2 Analysis
#'  da_maaslin <- read_humann(file_path, metadata) %>%
#'    run_maaslin2(fixed_effects = "ARM")
#'
#'  # Read HUMAnN3 & DESeq2 Analysis
#'  da_deseq <- read_humann(file_path, metadata) %>%
#'    run_deseq2(factor = "ARM")
#'
#'  # Volcano Plot
#'  plt1 <- plt_volcano(da_maaslin)
#'  plt2 <- plt_volcano(da_deseq)
#'
#'  # Check ggplot object
#'  testthat::expect_true("ggplot" %in% class(plt1))
#'  testthat::expect_true("ggplot" %in% class(plt2))
#'
#'  # Check number of layers
#'  testthat::expect_equal(length(plt1$layers), 4)
#'  testthat::expect_equal(length(plt2$layers), 4)
#'
#'  # Check labels
#'  labels1 <- c("coefficient", "-log10 p-value", "coefficient sign",
#'                "q-value", "yintercept", "xintercept", "feature")
#'  labels2 <- c("log2FoldChange", "-log10 p-value", "log2FC sign",
#'                "p-adjusted value", "yintercept", "xintercept", "feature")
#'
#'  testthat::expect_equal(labels1, unname(unlist(plt1$labels)))
#'  testthat::expect_equal(labels2, unname(unlist(plt2$labels)))
#'
#' # Check unknown input type
#' da_unknown <- da_maaslin %>%
#'  dplyr::mutate(da_method = "unknown")
#' testthat::expect_error(plt_volcano(da_unknown))
#'
#' @examples
#' # Def data paths
#' metadata <- system.file("extdata", "ex_meta.csv", package = "microfunk")
#' file_path <- system.file("extdata", "All_genefam_cpm_kegg.tsv", package = "microfunk")
#'
#' # Read HUMAnN3 & MaAsLin2 Analysis
#' da_result <-
#'   read_humann(file_path, metadata) %>%
#'   run_maaslin2(fixed_effects = "ARM")
#'
#' # Volcano Plot
#' plt_volcano(da_result)
plt_volcano <- function(da_result){

  input_type <- da_result %>%
    dplyr::pull(da_method) %>%
    unique()

  if(input_type == "maaslin2"){
    .maaslin_volcano(da_result = da_result)
  } else if(input_type == "deseq2"){
    .deseq_volcano(da_result = da_result)
  } else {
    cli::cli_abort(c(
      "x" = "Unknown input type.",
      "i" = "Please provide DA results of DESeq2 or MaAsLin2 analysis."))
  }

}

#' Plot MaAsLin2 Volcano Plot
#'
#' This function generates a volcano plot for MaAsLin2 differential abundance
#' analysis results.
#'
#' @param da_result A tibble containing the results of the MaAsLin2 DA analysis.
#'
#' @return A ggplot2 object representing the volcano plot.
.maaslin_volcano <- function(da_result){

  lower <- stats::quantile(da_result$coefficient, 0.01)
  upper <- stats::quantile(da_result$coefficient, 0.99)

  max_value <- max(abs(lower), abs(upper))

  da_result %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = coefficient, y = -log10(p_value))) +

    ggplot2::geom_point(
      ggplot2::aes(shape = factor(ifelse(q_value < 0.05, "s", "n")),
                   color = factor(ifelse(q_value > 0.05, "n",
                                         ifelse(coefficient > 0, "pos", "neg")))
      )) +

    ggplot2::xlim(-max_value, max_value) +

    ggplot2::scale_color_manual(values =
                                  c("pos" = "brown1", "neg" = "cornflowerblue",
                                    "n" = "grey"),
                                labels = c("pos" = "Positive", "neg" = "Negative",
                                           "n" = "Not significant")) +

    ggplot2::scale_shape_manual(values = c("s" = 20, "n" = 46),
                                labels = c("n" = "Not significant", "s" = "Significant (< 0.05)")) +

    ggplot2::geom_hline(
      yintercept = -log10(0.05),
      linetype = 2,
      linewidth = 0.5,
      colour = "grey"
    ) +

    ggplot2::geom_vline(
      xintercept = 0,
      linetype = 2,
      linewidth = 0.25,
      colour = "grey"
    ) +

    ggplot2::labs(x = "coefficient", y = "-log10 p-value", color = "coefficient sign",
                  shape = "q-value") +

    ggrepel::geom_text_repel(
      data = dplyr::filter(
        da_result, q_value < 0.05),
      ggplot2::aes(label = feature),
      size = 2.5
    ) +

    ggplot2::theme_minimal()

}


#' Plot DESeq2 Volcano Plot
#'
#' This function generates a volcano plot for DESeq2 differential abundance
#' analysis results.
#'
#' @param da_result A tibble containing the results of the DESeq2 DA analysis.
#'
#' @return A ggplot2 object representing the volcano plot.
.deseq_volcano <- function(da_result){

  da_result %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = log2FC, y = -log10(p_value))) +

    ggplot2::geom_point(
      ggplot2::aes(shape = factor(ifelse(padj_value < 0.05, "s", "n")),
                   color = factor(ifelse(padj_value > 0.05, "n",
                                         ifelse(log2FC > 0, "pos", "neg")))
      )) +

    ggplot2::scale_color_manual(values =
                                  c("pos" = "brown1", "neg" = "cornflowerblue",
                                    "n" = "grey"),
                                labels = c("pos" = "Positive", "neg" = "Negative",
                                           "n" = "Not significant")) +

    ggplot2::scale_shape_manual(values = c("s" = 20, "n" = 46),
                                labels = c("n" = "Not significant", "s" = "Significant (< 0.05)")) +

    ggplot2::geom_hline(
      yintercept = -log10(0.05),
      linetype = 2,
      linewidth = 0.5,
      colour = "grey"
    ) +

    ggplot2::geom_vline(
      xintercept = 0,
      linetype = 2,
      linewidth = 0.25,
      colour = "grey"
    ) +

    ggplot2::labs(x = "log2FoldChange", y = "-log10 p-value", color = "log2FC sign",
                  shape = "p-adjusted value") +

    ggrepel::geom_text_repel(
      data = dplyr::filter(
        da_result, padj_value < 0.05),
      ggplot2::aes(label = feature),
      size = 2.5
    ) +

    ggplot2::theme_minimal()


}

