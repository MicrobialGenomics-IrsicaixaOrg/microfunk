#' Read HUMAnN3 output file and create SummarizedExperiment objects
#'
#' This function reads a HUMAnN3 output file
#' (https://huttenhower.sph.harvard.edu/humann), which contains information
#' about the abundance of gene families or metabolic pathways of a given sample
#' in a TSV format. The function performs several checks on the structure and
#' content of the file, as well as reading and validating the metadata CSV file,
#' ensuring it contains the required information and matching sample IDs to the
#' provided TSV data file. Finally, it creates a SummarizedExperiment object
#' (https://bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html)
#' containing gene/pathway abundance values as the assay and metadata as the
#' column data.
#'
#' @param file_path Path to the HUMAnN3 output (gene families or pathway)
#'   collapsed TSV file. The file should contain aggregated abundance data
#'   across multiple samples. Each row represents a unique functional category
#'   (gene family / pathway) and each column corresponds to a sample.
#' @param metadata Path to the CSV file containing sample IDs and additional
#'   information, where each row corresponds to a sample and each column
#'   represents a specific annotation associated with that sample.
#'
#' @return A SummarizedExperiment object containing gene/pathway abundance
#'   values as the assay and metadata as the column data.
#' @export
#' @autoglobal
#' @tests
#' file_c <-
#'   system.file("extdata", "All_genefam_cpm_kegg.tsv", package =
#'   "microfunk")
#' metadata <-
#'   system.file("extdata", "ex_meta.csv", package =
#'   "microfunk")
#' file_r <-
#'   system.file("extdata", "All_genefam_rpk_kegg.tsv",
#'   package = "microfunk")
#'
#' # Test metadata errors
#' system.file("extdata", package = "microfunk") %>%
#'   list.files(pattern = "ex_meta_error", full.names = TRUE) %>%
#'   purrr::walk(~testthat::expect_error(read_humann(file_c, .x)))
#'
#' # Test HUMAnN3 file errors
#' system.file("extdata", package = "microfunk") %>%
#'   list.files(pattern = "All_genefam_cpm_kegg_error", full.names = TRUE) %>%
#'   purrr::walk(~ testthat::expect_error(read_humann(.x, metadata)))
#'
#' # Test output
#' testthat::expect_equal(class(read_humann(file_r, metadata))[1],
#'   "SummarizedExperiment")
#'
#' @examples
#' # Def data paths
#' metadata <-
#'   system.file("extdata", "ex_meta.csv", package = "microfunk")
#' file_path <-
#'   system.file("extdata", "All_genefam_cpm_kegg.tsv", package = "microfunk")
#'
#' # Read HUMAnN3
#' read_humann(file_path, metadata)
read_humann <- function(file_path, metadata){

  # check input paths
  .humann_path_checks(file_path, metadata)

  # read file
  tbl <-
    data.table::fread(file_path, sep = "\t") %>%
    tibble::as_tibble()

  # check humann tbl
  .humann_tbl_checks(tbl)

  # change column names
  tbl <-
    dplyr::rename_with(tbl, ~ stringr::str_remove(., "_Abundance.*$")) %>%
    dplyr::rename("function_id" = 1)

  # read and check metadata
  m_data <-
    data.table::fread(metadata) %>%
    tibble::as_tibble() %>%
    .humann_metadata_check(tbl) %>%
    data.frame(row.names = 1)

  # reorder by id
  tbl <- tbl %>%
    dplyr::arrange(function_id) %>%
    dplyr::select(function_id, rownames(m_data)) %>%
    data.frame(row.names = 1)

  # create SE object
  SummarizedExperiment::SummarizedExperiment(
    assays = list(humann = tbl),
    colData = m_data
  )
}


#' Check HUMAnN3 output file and metadata input paths
#'
#' This function checks that input paths of the HUMAnN3 output file and metadata
#' file exist and that metadata is contained in a CSV file.
#'
#' @param file_path Path to the HUMAnN3 output TSV file.
#' @param metadata Path to the metadata CSV file.
#'
#' @return 'NULL' if input paths pass the checks. If any issue is found, an
#'   error is raised.
#' @keywords internal
#' @autoglobal
#' @tests
#' nfile <-
#'   system.file("extdata", "nfile", package = "microfunk")
#' meta <-
#'   system.file("extdata", "ex_meta.csv", package = "microfunk")
#' file <-
#'   system.file("extdata", "All_genefam_cpm_kegg.tsv", package = "microfunk")
#'
#' testthat::expect_error(.humann_path_checks(metadata = meta))
#' testthat::expect_error(.humann_path_checks(file))
#' testthat::expect_error(.humann_path_checks(nfile, meta))
.humann_path_checks <- function(file_path, metadata) {
  # check input
  if (missing(file_path)) { cli::cli_abort(c("x" = "File path is missing")) }
  if (missing(metadata)) { cli::cli_abort(c("x" = "Metadata is missing")) }

  # check file path exists
  if (!file.exists(file_path)) { cli::cli_abort(c("x" = "File does not exist")) }

  # check metadata input
  if (!stringr::str_detect(metadata, ".csv$")) {
    cli::cli_abort(c(
      "x" = "Unsupported metadata file extension.",
      "i" = "Please provide a CSV file."
    ))
  }
}

#' Check validity of the tibble derived from the HUMAnN3 output file
#'
#' This function performs several checks on the tibble: its abundance values
#' correspond to either '# Gene Family' or '# Pathway', it is not empty, it
#' contains at least two samples and column classes (character/numeric) are
#' correct.
#'
#' @param humann_tbl A tibble containing the data from the HUMAnN3 output file.
#'
#' @return 'NULL' if all checks pass. If any issue is found, an error is raised.
#' @keywords internal
#' @autoglobal
.humann_tbl_checks <- function(humann_tbl) {
  # check header
  header <- colnames(humann_tbl)[1]
  if (!(
    stringr::str_detect(header, "^# Gene Family$") ||
    stringr::str_detect(header, "^# Pathway$")
  )) {
    cli::cli_abort(c("x" = "The header of the file should be '# Gene Family' or '# Pathway'."))
  }

  # check number of columns
  if (ncol(humann_tbl) < 2) {
    cli::cli_abort(c(
      "x" = "ncol(input) < 2.",
      "i" = "The file should contain at least two columns."
    ))
  }

  # check number of rows
  if (nrow(humann_tbl) <= 1) {
    cli::cli_abort(c(
      "x"="The file does not contain any data",
      "i"="Ensure the file is not empty."
    ))
  }

  # check column classes
  class_1 <- class(humann_tbl[[1]]) == "character"
  class_2 <- purrr::map_lgl(humann_tbl[-1], ~ class(.x) == "numeric") %>% unname()

  if (!class_1 || !all(class_2)) {
    cli::cli_abort(c(
      "x" = "Column classes are not correct.",
      "i" = "First column should be 'character' and the rest should be 'numeric'."
    ))
  }
}

#' Check metadata consistency with HUMAnN3 output file
#'
#' This function verifies the consistency of the provided metadata with the
#' HUMAnN3 file. It ensures metadata contains sufficient annotations and that
#' sample IDs in the metadata match the column names of the HUMAnN3 file.
#'
#' @param meta_df A tibble containing sample metadata.
#' @param humann_tbl A tibble containing the data from the HUMAnN3 output file.
#'
#' @return The input metadata tibble if all checks pass.
#' @keywords internal
#' @autoglobal
.humann_metadata_check <- function(meta_df, humann_tbl) {
  # check metadata content
  if (ncol(meta_df) < 2) {
    cli::cli_abort(c(
      "x" = "ncol(metadata) < 2.",
      "i" = "Metadata should contain at least one annotation."
    ))
  }

  # check sample ids
  if (length(meta_df[[1]]) != length(names(humann_tbl)[-1])) {
    cli::cli_abort(c("x" = "Different number of samples in assay and metadata."))
  }
  if (!all(meta_df[[1]] %in% names(humann_tbl))) {
    cli::cli_abort(c("x" = "Sample names in assay and metadata do not match."))
  }

  meta_df
}

