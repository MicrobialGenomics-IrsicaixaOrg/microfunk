#' Read HUMAnN3 output file and create SummarizedExperiment objects
#'
#' This function reads an output file from HUMAnN3 and creates
#' SummarizedExperiment objects containing normalized abundance values for
#' groups (gene families or pathway categories) and stratifications (breakdown
#' of groups by species).
#'
#' @param file_path Path to the HUMAnN3 output TSV file.
#' @param metadata A CSV file containing sample information.
#'
#' @return A list containing two SummarizedExperiment objects: one for groups
#'   and one for stratifications.
#'
#' @export
#'
#' @examples
#' # Def data paths
#' metadata <- system.file("extdata", "ex_meta.csv", package = "microfunk")
#' file_path <-
#'   system.file("extdata", "All_genefam_cpm_kegg.tsv", package = "microfunk")
#'
#' # Read HUMAnN3
#' draft_read_humann(file_path, metadata)
draft_read_humann <- function(file_path, metadata){

  # check input paths
  .humann_path_checks(file_path, metadata)

  # read file
  tbl <- data.table::fread(file_path, sep = "\t") %>% tibble::as_tibble()

  # check humann tbl
  .humann_tbl_checks(tbl)

  # change column names
  tbl <- dplyr::rename_with(tbl, ~ stringr::str_remove(., "_Abundance.*$"))

  # check data normalization
  norm_method <- .extract_norm_method(tbl)

  if (norm_method == "rpk") {
    tbl <- .rpk2cpm(tbl)
  }

  # read and check metadata
  m_data <-
    data.table::fread(metadata) %>%
    tibble::as_tibble() %>%
    .humann_metadata_check(tbl) %>%
    data.frame(row.names = 1)

  # reorder by id
  tbl <- tbl %>%
    dplyr::select(1, rownames(m_data)) %>%
    data.frame(row.names = 1)

  # create SE object
  SummarizedExperiment::SummarizedExperiment(
    assays = tbl,
    colData = m_data
  )
}


.humann_path_checks <- function(file_path, metadata) {
  # check input
  if (missing(file_path)) { cli::cli_abort(c("x" = "File path is missing")) }
  if (missing(metadata)) { cli::cli_abort(c("x" = "Metadata is missing")) }

  # check file path exists
  if (!file.exists(file_path)) { cli::cli_abort(c("x" = "File does not exist")) }

  # check metadata input
  if (tools::file_ext(metadata) != "csv"){
    cli::cli_abort(c(
      "x" = "Unsupported metadata file extension.",
      "i" = "Please provide a CSV file."
    ))
  }
}

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
  class_2 <- purrr::map_lgl(tbl[-1], ~ class(.x) == "numeric") %>% unname()

  if (!class_1 || !all(class_2)) {
    cli::cli_abort(c(
      "x" = "Column classes are not correct.",
      "i" = "First column should be 'character' and the rest should be 'numeric'."
    ))
  }
}

.extract_norm_method <- function(tbl) {
  # check data normalization
  sum_counts <-
    dplyr::filter(tbl, !stringr::str_detect(.data[[header]],"[|]")) %>%
    dplyr::select(-1) %>%
    colSums()

  classification <-
    ifelse(dplyr::between(sum_counts, 99900, 1001000), "cpm", "rpk") %>%
    unique()

  # Check if all are same normalization
  if (length(classification) != 1) {
    cli::cli_abort(c(
      "x" = "Columns have mixed units.",
      "i" = "Please ensure all samples have either 'rpk' or 'cpm' units."
    ))
  }

  classification
}

.rpk2cpm <- function(tbl) {
  # scale by cat
  tbl_group <-
    dplyr::filter(tbl, !stringr::str_detect(.data[[header]],"[|]")) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ . / sum(.) * 1e3))

  # scale by strat
  tbl_strat <-
    dplyr::filter(tbl, stringr::str_detect(.data[[header]], "[|]|UNMAPPED")) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ . / sum(.) * 1e3))

  dplyr::bind_rows(tbl_group, tbl_strat) %>%
    dplyr::arrange(stringr::str_remove_all(`# Gene Family`, "[|].*"))
}

.humann_metadata_check <- function(meta_df, humann_tbl) {
  # check metadata content
  if (ncol(meta_df) < 2) {
    cli::cli_abort(c(
      "x" = "ncol(metadata) < 2.",
      "i" = "Metadata should contain at least one annotation."
    ))
  }

  # check sample ids
  if (!all(meta_df[[1]] %in% names(humann_tbl))) {
    cli::cli_abort(c("x" = "Sample names in assay and metadata do not match."))
  }


  meta_df
}
