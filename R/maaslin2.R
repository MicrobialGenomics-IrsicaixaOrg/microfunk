#' Run MaAsLin2 Analysis Quietly
#'
#' This function performs a MaAsLin2 analysis
#' (https://huttenhower.sph.harvard.edu/maaslin/), which finds associations
#' between microbiome meta-omics features and metadata. Features consist on the
#' abundance of gene families or metabolic pathways of multiple samples,
#' obtained using HUMAnN3 (https://huttenhower.sph.harvard.edu/humann).
#'
#' The function takes a SummarizedExperiment object containing normalized gene/
#' pathway abundance values as the assay and metadata as the column data, as
#' well as several parameter options to run the analysis. It creates a temporary
#' output folder to store created files and returns the result of the analysis
#' on the console.
#'
#' @param se A SummarizedExperiment object containing normalized abundance
#'   values as the assay and metadata as the column data.
#' @param output Character string specifying the output directory path. Default
#'   is a tmp folder.
#' @param min_abundance Minimum abundance for each feature. Default is 0.0.
#' @param min_prevalence Minimum percent of samples for which a feature is
#'   detected at minimum abundance. Default is 0.1.
#' @param min_variance Minimum variance threshold for features to be included.
#'   Default is 0.0.
#' @param normalization Normalization method to apply. Default is 'NONE'.
#' @param transform Transformation method to be used. Default is 'NONE'.
#' @param analysis_method Analysis method to apply. Default is 'LM'.
#' @param max_significance Q-value threshold for significance. Default is 0.25.
#' @param random_effects Random effects to be included in the analysis
#'   (comma-delimited for multiple effects). Default is NULL.
#' @param fixed_effects Fixed effects to be included in the analysis
#'   (comma-delimited for multiple effects). Default is NULL.
#' @param correction Correction method for computing the Q-value. Default is
#'   'BH'.
#' @param standardize Apply Z-score so continuous metadata are on the same
#'   scale. Default is FALSE.
#' @param plot_heatmap Generate a heatmap for significant associations. Default
#'   is TRUE.
#' @param heatmap_first_n In the heatmap, plot only top N features with
#'   significant associations. Default is 50.
#' @param plot_scatter Generate scatter plots for the significant associations.
#'   Default is TRUE.
#' @param max_pngs Set the maximum number of scatter plots for significant
#'   associations to save as png files. Default is 10.
#' @param save_scatter Save all scatter plot ggplot objects to an RData file.
#'   Default is FALSE.
#' @param save_models Return the full model outputs and save to an RData file.
#'   Default is FALSE.
#' @param reference Factor to use as reference for categorical variables,
#'   provided as a string of 'variable,reference' (semi-colon delimited for
#'   multiple variables).
#'
#' @return Results of the MaAsLin2 analysis on the console.
#' @export
#' @autoglobal
#' @tests
#' # Def data paths
#' metadata <- system.file("extdata", "reduced_meta.csv", package = "microfunk")
#' file_path <-
#'   system.file("extdata", "reduced_genefam_cpm_kegg.tsv", package = "microfunk")
#'
#' # Read HUMAnN3 & MaAsLin2 Analysis
#' da_result <-
#'   read_humann(file_path, metadata) %>%
#'   run_maaslin2(fixed_effects = "ARM")
#'
#' # Test number of significant associations
#' da_result$results %>%
#'   dplyr::filter(qval < 0.25) %>%
#'   nrow() %>%
#'   testthat::expect_equal(63)
#'
#' # Test P-values
#' da_result$results %>%
#'   dplyr::pull(pval) %>%
#'   mean() %>%
#'   round(3) %>%
#'   testthat::expect_equal(0.538)
#'
#' # Test Q-values
#' da_result$results %>%
#'   dplyr::pull(qval) %>%
#'   mean() %>%
#'   round(3) %>%
#'   testthat::expect_equal(0.802)
#'
#' @examples
#' # Read HUMAnN3 & MaAsLin2 Analysis
#' da_result <- read_humann(
#'   file_path = system.file("extdata", "reduced_genefam_cpm_kegg.tsv", package = "microfunk"),
#'   metadata = system.file("extdata", "reduced_meta.csv", package = "microfunk")
#' ) %>% run_maaslin2(fixed_effects = "ARM")
#'
#' da_result
run_maaslin2 <- function(se,
                         output = paste0(tempdir(), "/output_folder"),
                         min_abundance = 0.0,
                         min_prevalence = 0.1,
                         min_variance = 0.0,
                         normalization = "NONE",
                         transform = "NONE",
                         analysis_method = "LM",
                         max_significance = 0.25,
                         random_effects = NULL,
                         fixed_effects = NULL,
                         correction = "BH",
                         standardize = TRUE,
                         plot_heatmap = TRUE,
                         heatmap_first_n = 50,
                         plot_scatter = TRUE,
                         max_pngs = 10,
                         save_scatter = FALSE,
                         save_models = FALSE,
                         reference = NULL) {

  # Retrieve tbls from SE object
  input_meta <- SummarizedExperiment::colData(se) %>% as.data.frame()

  # Filter input data
  input_group <-
    SummarizedExperiment::assays(se)$humann %>%
    tibble::as_tibble(rownames = "function_id") %>%
    dplyr::filter(!stringr::str_detect(function_id, "[|]")) %>%
    tibble::column_to_rownames(var = "function_id")

  # Run quietly MaAsLin2 analysis
  .maaslin2_quietly(
    input_group = input_group,
    input_meta = input_meta,
    output = output,
    min_abundance = min_abundance,
    min_prevalence = min_prevalence,
    min_variance = min_variance,
    normalization = normalization,
    transform = transform,
    analysis_method = analysis_method,
    max_significance = max_significance,
    random_effects = random_effects,
    fixed_effects = fixed_effects,
    correction = correction,
    standardize = standardize,
    plot_heatmap = plot_heatmap,
    heatmap_first_n = heatmap_firt_n,
    plot_scatter = plot_scatter,
    max_pngs = max_pngs,
    save_scatter = save_scatter,
    save_models = save_models,
    reference = reference
  )
}

#' Internal Function to Perform MaAsLin2 Analysis Quietly
#'
#' This function performs the MaAsLin2 analysis without printing any messages.
#' It uses a nested function `f_quietly()`, which takes the same arguments as
#' the parent function`.maaslin2_quietly()` and performs the actual analysis.
#' This includes data preparation, filtering and calling the `Maaslin2()`
#' function from the `Maaslin2` package. Then, the `quietly()` function from the
#' `purrr` package takes the `f_quietly()`as an argument and returns a modified
#' version that does not print any messages, warnings or errors to the console
#' (any output messages generated by the MaAsLin2 analysis are supressed).
#'
#' @param input_group  A tibble containing a filtered subset of the input data
#'   with the abundance values of main categories of gene families or metabolic
#'   pathways.
#' @param input_meta A tibble containing the metadata associated with each
#'   sample in the analysis.
#' @param output Character string specifying the output directory path. Default
#'   is a tmp folder.
#' @param min_abundance Minimum abundance for each feature. Default is 0.0.
#' @param min_prevalence Minimum percent of samples for which a feature is
#'   detected at minimum abundance. Default is 0.1.
#' @param min_variance Minimum variance threshold for features to be included.
#'   Default is 0.0.
#' @param normalization Normalization method to apply. Default is 'NONE'.
#' @param transform Transformation method to be used. Default is 'NONE'.
#' @param analysis_method Analysis method to apply. Default is 'LM'.
#' @param max_significance Q-value threshold for significance. Default is 0.25.
#' @param random_effects Random effects to be included in the analysis
#'   (comma-delimited for multiple effects). Default is NULL.
#' @param fixed_effects Fixed effects to be included in the analysis
#'   (comma-delimited for multiple effects). Default is NULL.
#' @param correction Correction method for computing the Q-value. Default is
#'   'BH'.
#' @param standardize Apply Z-score so continuous metadata are on the same
#'   scale. Default is FALSE.
#' @param plot_heatmap Generate a heatmap for significant associations. Default
#'   is TRUE.
#' @param heatmap_first_n In the heatmap, plot only top N features with
#'   significant associations. Default is 50.
#' @param plot_scatter Generate scatter plots for the significant associations.
#'   Default is TRUE.
#' @param max_pngs Set the maximum number of scatter plots for significant
#'   associations to save as png files. Default is 10.
#' @param save_scatter Save all scatter plot ggplot objects to an RData file.
#'   Default is FALSE.
#' @param save_models Return the full model outputs and save to an RData file.
#'   Default is FALSE.
#' @param reference Factor to use as reference for categorical variables,
#'   provided as a string of 'variable,reference' (semi-colon delimited for
#'   multiple variables).
#'
#' @return Results of the MaAsLin2 analysis.
#' @autoglobal
#' @keywords internal
#' @noRd
.maaslin2_quietly <- function(input_group,
                               input_meta,
                               output,
                               min_abundance,
                               min_prevalence,
                               min_variance,
                               normalization,
                               transform,
                               analysis_method,
                               max_significance,
                               random_effects,
                               fixed_effects,
                               correction,
                               standardize,
                               plot_heatmap,
                               heatmap_first_n,
                               plot_scatter,
                               max_pngs,
                               save_scatter,
                               save_models,
                               reference) {

  f_quietly <- function(input_group = input_group,
                        input_meta = input_meta,
                        output = output,
                        min_abundance = min_abundance,
                        min_prevalence = min_prevalence,
                        min_variance = min_variance,
                        normalization = normalization,
                        transform = transform,
                        analysis_method = analysis_method,
                        max_significance = max_significance,
                        random_effects = random_effects,
                        fixed_effects = fixed_effects,
                        correction = correction,
                        standardize = standardize,
                        plot_heatmap = plot_heatmap,
                        heatmap_first_n = heatmap_firt_n,
                        plot_scatter = plot_scatter,
                        max_pngs = max_pngs,
                        save_scatter = save_scatter,
                        save_models = save_models,
                        reference = reference) {

    # Run MaAsLin2 analysis
    fit_data <- Maaslin2::Maaslin2(
      input_data = input_group,
      input_metadata = input_meta,
      output = output,
      fixed_effects = fixed_effects,
      min_abundance = min_abundance,
      min_prevalence = min_prevalence,
      min_variance = min_variance,
      normalization = normalization,
      transform = transform,
      analysis_method = analysis_method,
      max_significance = max_significance,
      random_effects = random_effects,
      correction = correction,
      standardize = standardize,
      reference = reference,
      plot_heatmap = plot_heatmap,
      plot_scatter = plot_scatter)
  }

  purrr::quietly(f_quietly)(
    input_group = input_group,
    input_meta = input_meta,
    output = output,
    min_abundance = min_abundance,
    min_prevalence = min_prevalence,
    min_variance = min_variance,
    normalization = normalization,
    transform = transform,
    analysis_method = analysis_method,
    max_significance = max_significance,
    random_effects = random_effects,
    fixed_effects = fixed_effects,
    correction = correction,
    standardize = standardize,
    plot_heatmap = plot_heatmap,
    heatmap_first_n = heatmap_firt_n,
    plot_scatter = plot_scatter,
    max_pngs = max_pngs,
    save_scatter = save_scatter,
    save_models = save_models,
    reference = reference
  )$result
}
