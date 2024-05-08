#' Generate a Volcano Plot for MaAsLin2 Analysis Results
#'
#' This function generates a volcano plot based on the results of a MaAsLin2
#' analysis in order to visualize the differential abundance of microbial gene
#' families or pathways (features) associated with specific fixed effects.
#'
#' The x-axis represents the coefficient or effect size of each feature,
#' indicating the magnitude and direction of change associated with the chosen
#' fixed effect. The y-axis represents the statistical significance of the
#' differences (negative logarithm of the p-values).
#'
#' Each point on the plot corresponds to a gene family or pathway, with its
#' position determined by both its effect size and statistical significance.
#' Only features with a q-value < 0.05 are deemed significant and therefore
#' colored, in red if its abundance increases or blue if its abundance
#' decreases.
#'
#' @param da_result The result of the MaAsLin2 differential analysis.
#'
#' @return A volcano plot visualizing differential abundance of microbial
#'   features.
#' @export
#' @autoglobal
#' @tests
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
#' plt <- plt_volcano(da_result)
#'
#' # Check ggplot object
#' testthat::expect_true("ggplot" %in% class(plt))
#'
#' # Check number of layers
#' testthat::expect_equal(length(plt$layers), 4)
#'
#' # Check labels
#' labels <- c("coefficient", "-log10 p-value", "coefficient sign", "q-value",
#'             "yintercept", "xintercept", "feature")
#'
#' testthat::expect_equal(labels, unname(unlist(plt$labels)))
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

  plt_df <- da_result$results %>% tibble::as_tibble()

  lower <- stats::quantile(plt_df$coef, 0.01)
  upper <- stats::quantile(plt_df$coef, 0.99)

  max_value <- max(abs(lower), abs(upper))

  plt_df %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = coef, y = -log10(pval))) +

    ggplot2::geom_point(
      ggplot2::aes(shape = factor(ifelse(qval < 0.05, "s", "n")),
          color = factor(ifelse(qval > 0.05, "n",
                                ifelse(coef > 0, "pos", "neg")))
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
        plt_df, qval < 0.05),
      ggplot2::aes(label = feature),
      size = 2.5
    ) +

    ggplot2::theme_minimal()
}
