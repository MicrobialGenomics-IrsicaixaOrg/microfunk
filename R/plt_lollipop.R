#' Generate a Lollipop Plot for MaAsLin2 Analysis Results
#'
#' This function generates a lollipop plot based on the results of a MaAsLin2 or
#' DESeq2 analysis to visualize the differential abundance of microbial gene
#' families or pathways (features) associated with specific categorical
#' variables. It only includes significant features with a q-value < 0.05 (for
#' MaAsLin2) or a p-adjusted value < 0.05 (for DESeq2).
#'
#' The x-axis represents the coefficient/effect size (for MaAsLin2) or
#' log2FoldChange (for DESeq2) of each feature, indicating the magnitude and
#' direction of change associated with the chosen categorical variable. The
#' y-axis represents the features ordered by their coefficient values (for
#' MaAsLin2) or log2FoldChange (for DESeq2).
#'
#' Each feature is depicted as a line segment originating from the origin (0 on
#' the x-axis) and extending to the corresponding coefficient/log2FC value on
#' the x-axis. Negative coefficient/log2FC values indicate a decrease in
#' abundance, while positive ones indicate an increase in abundance. The color
#' of each feature segment is determined by the statistical significance of the
#' differences, represented by a color gradient scale using the negative
#' logarithm of p-values.
#'
#' @param da_result A tibble containing the result of the differential analysis.
#' @param n Maximum number of significant associations to be displayed on the
#'   plot in ascending order of q-value (for MaAsLin2) or p-adjusted value (for
#'   DESeq2).
#'
#' @return A lollipop plot visualizing differential abundance of microbial
#'   features.
#' @export
#' @autoglobal
#' @import ggplot2
#' @tests
#' # Def data paths
#' metadata <- system.file("extdata", "reduced_meta.csv", package = "microfunk")
#' file_path <- system.file("extdata", "reduced_genefam_cpm_kegg.tsv", package = "microfunk")
#'
#' # Read HUMAnN3 & MaAsLin2 Analysis
#' da_maaslin <- read_humann(file_path, metadata) %>%
#'  run_maaslin2(fixed_effects = "ARM")
#'
#' # Read HUMAnN3 & DESeq2 Analysis
#' da_deseq <- read_humann(file_path, metadata) %>%
#'  run_deseq2(factor = "ARM")
#'
#' # Lollipop Plot
#' plt1 <- plt_lollipop(da_maaslin, n = 20)
#' plt2 <- plt_lollipop(da_deseq, n = 20)
#'
#' # Check ggplot object
#' testthat::expect_true("ggplot" %in% class(plt1))
#' testthat::expect_true("ggplot" %in% class(plt2))
#'
#' # Check number of layers
#' testthat::expect_equal(length(plt1$layers), 2)
#' testthat::expect_equal(length(plt2$layers), 2)
#'
#' # Check labels
#' labels1 <- c("coefficient", "feature", "-log10(p_value)", "-log10(p_value)",
#'              "xend", "feature")
#' labels2 <- c("log2FoldChange", "feature", "-log10(p_value)", "-log10(p_value)",
#'             "xend", "feature")
#' testthat::expect_equal(labels1, unname(unlist(plt1$labels)))
#' testthat::expect_equal(labels2, unname(unlist(plt2$labels)))
#'
#' # Check unknown input type
#' da_unknown <- da_maaslin %>%
#'  dplyr::mutate(da_method = "unkown")
#' testthat::expect_error(plt_lollipop(da_unknown, n = 20))
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
#' plt_lollipop(da_result, n = 20)
plt_lollipop <- function(da_result, n = 20){

  input_type <- da_result %>%
    dplyr::pull(da_method) %>%
    unique()

  if(input_type == "maaslin2"){
    .maaslin_lollipop(da_result = da_result, n = n)
  } else if(input_type == "deseq2"){
    .deseq_lollipop(da_result = da_result, n = n)
  } else {
    cli::cli_abort(c(
      "x" = "Unknown input type.",
      "i" = "Please provide DA results of DESeq2 or MaAsLin2 analysis."))
  }
}


#' Plot MaAsLin2 Lollipop Plot
#'
#' This function generates a lollipop plot for MaAsLin2 differential abundance
#' analysis results.
#'
#' @param da_result A tibble containing the results of the MaAsLin2 DA analysis.
#' @param n Maximum number of significant associations to be displayed on the
#'   plot in ascending order of q-value.
#'
#' @autoglobal
#' @return A ggplot2 object representing the lollipop plot.
.maaslin_lollipop <- function(da_result, n) {

  da_result <- da_result %>%
    dplyr::filter(signif == TRUE) %>%
    dplyr::arrange(q_value)

  if (nrow(da_result) > n) {
    cli::cli_inform(c(
      "!" ="Features present in all methods are greater than the cutoff n = {n}",
      "i" = "The top {n} significant features will be used"
    ))
  }

  da_result <- dplyr::slice(da_result, 1:n)

  lim <- da_result %>%
    dplyr::filter(q_value < 0.05) %>%
    dplyr::pull(p_value) %>%
    -log10(.) %>%
    abs() %>%
    max() %>%
    c(-.,.)

  da_result %>%
    dplyr::filter(q_value < 0.05) %>%
    ggplot(aes(
      x = coefficient,
      y = stats::reorder(feature, coefficient),
      fill = -log10(p_value),
      color = -log10(p_value)
    )) +
    geom_segment(aes(xend = 0, yend = feature), linewidth = 0.5) +
    geom_point(shape = 21, size = 3) +
    scale_color_gradientn(
      colors = RColorBrewer::brewer.pal(11, "Spectral"),
      limits = lim
    ) +
    scale_fill_gradientn(
      colors = RColorBrewer::brewer.pal(11, "Spectral"),
      limits = lim
    ) +
    labs(x = "coefficient", y = "feature") +
    theme(axis.text.y = element_text(size = 7.5)) +
    theme_minimal()
}

#' Plot DESeq2 Lollipop Plot
#'
#' This function generates a lollipop plot for DESeq2 differential abundance
#' analysis results.
#'
#' @param da_result A tibble containing the results of the DESeq2 DA analysis.
#' @param n Maximum number of significant associations to be displayed on the
#'   plot in ascending order of p-adjusted value.
#'
#' @autoglobal
#' @return A ggplot2 object representing the lollipop plot.
.deseq_lollipop <- function(da_result, n){

  da_result <- da_result %>%
    dplyr::filter(signif == TRUE) %>%
    dplyr::arrange(padj_value)

  if (nrow(da_result) > n) {
    cli::cli_inform(c(
      "!" = "Features present in all methods are greater than the cutoff n = {n}",
      "i" = "The top {n} significant features will be used"
    ))
  }

  da_result <- dplyr::slice(da_result, 1:n)

  lim <- da_result %>%
    dplyr::filter(padj_value < 0.05) %>%
    dplyr::pull(p_value) %>%
    -log10(.) %>%
    abs() %>%
    max() %>%
    c(-.,.)

  da_result %>%
    dplyr::filter(padj_value < 0.05) %>%
    ggplot(mapping = aes(
      x = log2FC,
      y = stats::reorder(feature, log2FC),
      fill = -log10(p_value),
      color = -log10(p_value)
    )) +
    geom_segment(aes(xend = 0, yend = feature), linewidth = 0.5) +
    geom_point(shape = 21, size = 3) +
    scale_color_gradientn(
      colors = RColorBrewer::brewer.pal(11, "Spectral"),
      limits = lim
    ) +
    scale_fill_gradientn(
      colors = RColorBrewer::brewer.pal(11, "Spectral"),
      limits = lim
    ) +
    labs(x = "log2FoldChange", y = "feature") +
    theme(axis.text.y = element_text(size = 7.5)) +
    theme_minimal()
}

