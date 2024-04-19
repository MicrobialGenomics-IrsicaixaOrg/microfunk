#' Run MaAsLin2 Analysis Quietly
#'
#' This function performs a MaAsLin2 analysis
#' (https://huttenhower.sph.harvard.edu/maaslin/), which finds associations
#' between microbiome meta-omics features and metadata. Features consist on the
#' abundance of gene families or metabolic pathways of multiple samples,
#' obtained using HUMAnN3 (https://huttenhower.sph.harvard.edu/humann). The
#' function takes two arguments, the path to the HUMAnN3 output file and the
#' metadata file, and applies linear models to asses the relationship between
#' microbial abundance data and metadata. It creates an output folder with a
#' list of significant associations in the `significant_results.tsv` file as
#' well as plots generated for each significant association, among other files.
#'
#' @param file Path to the HUMAnN3 output (gene families or pathway) collapsed
#'   TSV file. The file should contain aggregated abundance data across multiple
#'   samples. Each row represents a unique functional category (gene family /
#'   pathway) and each column corresponds to a sample.
#' @param metadata Path to the CSV file containing sample IDs and additional
#'   information, where each row corresponds to a sample and each column
#'   represents a specific annotation associated with that sample.
#'
#' @return Results of the MaAsLin2 analysis.
#' @export
#' @autoglobal
#'
#' @examples
#' # Def data paths
#' metadata <-
#'   system.file("extdata", "ex_meta.csv", package = "microfunk")
#' file_path <-
#'   system.file("extdata", "All_genefam_cpm_kegg.tsv", package = "microfunk")
#'
#' # MaAsLin2 Analysis
#' maaslin2_quietly(file_path, metadata)
maaslin2_quietly <- function(file, metadata){
  run_maaslin2 <- function(file, metadata) {

    # Read files and create SE object
    result <- read_humann(file, metadata)

    # Retrieve tbls from SE object
    data <- SummarizedExperiment::assays(result)[[1]] %>% tibble::rownames_to_column(var = "function_id") %>% tibble::as_tibble()
    meta <- SummarizedExperiment::colData(result) %>% as.data.frame()

    # Filter input data
    group <- data %>%
      dplyr::filter(!stringr::str_detect(function_id, "[|]")) %>%
      tibble::column_to_rownames(var = "function_id")

    # Run MaAsLin2 analysis
    fit_data <- Maaslin2::Maaslin2(
      input_data = group,
      input_metadata = meta,
      output = "output_folder",
      transform = "NONE",
      fixed_effects = colnames(meta))

  }
  purrr::quietly(run_maaslin2)(file, metadata)$result

}

