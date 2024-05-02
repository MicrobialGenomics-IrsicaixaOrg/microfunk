#' Generate a Lollipop Plot for MaAsLin2 Analysis Results
#'
#' This function generates a lollipop plot based on the results of a MaAsLin2
#' analysis to visualize the differential abundance of microbial gene families
#' or pathways (features) associated with specific fixed effects. It only
#' includes significant features with a q-value < 0.05.
#'
#' The x-axis represents the coefficient or effect size of each feature,
#' indicating the magnitude and direction of change associated with the chosen
#' fixed effect. The y-axis represents the features ordered by their coefficient
#' values.
#'
#' Each feature is depicted as a line segment originating from the origin (0 on
#' the x-axis) and extending to the corresponding coefficient value on the
#' x-axis. Negative coefficient values indicate a decrease in abundance, while
#' positive ones indicate an increase in abundance. The color of each feature
#' segment is determined by the statistical significance of the differences,
#' represented by a color gradient scale using the negative logarithm of
#' p-values.
#'
#' @param da_result The result of the MaAsLin2 differential analysis.
#'
#' @return A lollipop plot visualizing differential abundance of microbial
#'   features.
#' @export
#' @autoglobal
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
#' # Lollipop Plot
#' plt_lollipop(da_result)
plt_lollipop <- function(da_result){

  plt_df <- da_result$results %>% tibble::as_tibble()

  # Calculate limits for color gradient according to p-value
  lim <- plt_df %>%
    dplyr::filter(qval < 0.05) %>%
    dplyr::pull(pval) %>%
    -log10(.) %>%
    abs() %>%
    max() %>%
    c(-.,.)


  # Lollipop plot
  plt_df %>%

    # Include only significant results
    dplyr::filter(qval < 0.05) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = coef, y = stats::reorder(feature, coef),
                                           fill = -log10(pval),
                                           color = -log10(pval))) +

    # Segment connecting origin to corresponding feature
    ggplot2::geom_segment(ggplot2::aes(xend = 0, yend = feature),
                          linewidth = 0.5) +

    # Points at the end of each segment (color = significance)
    ggplot2::geom_point(shape = 21, size = 3) +

    # Set color gradient with calculated limits
    ggplot2::scale_color_gradientn(colors =
                                     RColorBrewer::brewer.pal(11, "Spectral"),
                                   limits = lim) +
    ggplot2::scale_fill_gradientn(colors =
                                    RColorBrewer::brewer.pal(11, "Spectral"),
                                  limits = lim) +

    # Add labels
    ggplot2::labs(x = "coefficient", y = "feature") +

    # Change y-axis feature label size
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7.5)) +

    ggplot2::theme_minimal()
}
