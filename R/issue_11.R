#' Read HUMAnN3 output file and create SummarizedExperiment objects
#'
#' This function reads an output file from HUMAnN3 and creates SummarizedExperiment objects
#' containing normalized abundance values for groups (gene families or pathway categories)
#' and stratifications (breakdown of groups by species).
#'
#' @param file_path Path to the HUMAnN3 output TSV file.
#' @param metadata A CSV file containing sample information.
#'
#' @return A list containing two SummarizedExperiment objects: one for groups and one for stratifications.
#'
#' @export
#'
#' @examples
#' file_path <-
#'   system.file("extdata", "All_genefam_cpm_kegg.tsv", package = "microfunk")
#' metadata <-
#'   system.file("extdata", "ex_meta.csv", package = "microfunk")
#' draft_read_humann(file_path, metadata)
draft_read_humann <- function(file_path, metadata){

  # check input
  if (missing(file_path)){cli::cli_abort(c("x"="File path is missing"))}
  if (missing(metadata)){cli::cli_abort(c("x"="Metadata is missing"))}

  # check file path exists
  if (!file.exists(file_path)) {
    cli::cli_abort(c("x"="File does not exist"))
  }

  # check metadata input
  if (tools::file_ext(metadata) != "csv"){
    cli::cli_abort(c("x"="Unsupported metadata file extension.",
                     "i"="Please provide a CSV file."))
  }

  # read file
  tbl <-
    data.table::fread(file_path, sep="\t") %>%
    tibble::as_tibble()

  # check header
  header <- colnames(tbl)[1]
  if (!(stringr::str_detect(header, "^# Gene Family$") ||
        stringr::str_detect(header, "^# Pathway$"))) {
    cli::cli_abort(c("x"="The header of the file should be '# Gene Family' or '# Pathway'."))
  }

  # check number of columns
  if (ncol(tbl) < 2) { cli::cli_abort(c("x"="ncol(input) < 2.",
                        "i"="The file should contain at least two columns.")) }

  # check number of rows
  if (nrow(tbl) <= 1) {
    cli::cli_abort(c("x"="The file does not contain any data",
                    "i"="Ensure the file is not empty."))
  }

  # check column classes
  class_1 <- class(tbl[[1]]) == "character"
  class_2 <- purrr::map_lgl(tbl[-1], ~class(.x) == "numeric") %>%
    unname()

  if (!class_1 || !all(class_2)) {
    cli::cli_abort(c("x" = "Column classes are not correct.",
               "i" = "First column should be 'character' and the rest should be 'numeric'."))
  }

  # change column names
  pattern <- "_Abundance.*$"
  names <- colnames(tbl) %>%
    stringr::str_remove(pattern)
  colnames(tbl) <- names

  # new tibble without stratification
  tbl_group <- tbl %>%
    dplyr::filter(!stringr::str_detect(.data[[header]],"[|]")) %>%
    dplyr::select(gtools::mixedsort(names(.))) %>%
    tibble::as_tibble()

  # new tibble with stratified values
  tbl_strat <- tbl %>%
    dplyr::filter(stringr::str_detect(.data[[header]], "[|]")) %>%
    dplyr::select(gtools::mixedsort(names(.))) %>%
    tibble::as_tibble()

  # check data normalization
  sum_counts <- tbl_group %>%
    dplyr::select(-1) %>%
    colSums() %>%
    unname()

  # classify columns in rpk/cpm
  classification <- ifelse(dplyr::between(sum_counts, 99900, 1001000), "cpm", "rpk")

  # check for different units across columns
  mixed <- length(unique(classification)) == 1
  if (!mixed){cli::cli_abort(c("x"="Columns have mixed units.",
              "i"="Please ensure all samples have either 'rpk' or 'cpm' units."))}

  # normalize rpk to cpm if needed
  if(unique(classification) == "rpk"){

    # scale groups
    norm_group <- tbl_group %>%
      dplyr::select(-1) %>%
      purrr::map2_dfc(., sum_counts, ~ .x / .y * 1e6) %>%
      dplyr::bind_cols(tbl_group[1],.) %>%
      tibble::as_tibble()

    # scale stratification
    norm_strat <- tbl_strat %>%
      dplyr::select(-1) %>%
      purrr::map2_dfc(., sum_counts, ~ .x / .y * 1e6) %>%
      dplyr::bind_cols(tbl_strat[1],.) %>%
      tibble::as_tibble()
  }
  norm_group <- tbl_group
  norm_strat <- tbl_strat

  # read metadata
  m_data <- data.table::fread(metadata) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(readr::parse_number(.[[1]]))

  # check metadata content
  if (ncol(m_data) < 2) { cli::cli_abort(c("x"="ncol(metadata) < 2.",
                       "i"="Metadata should contain at least one annotation.")) }

  names_match <- identical(m_data[[1]], names(tbl_group)[-1])
  if (!names_match){cli::cli_abort(c(
    "x"="Sample names in assay and metadata do not match."))}


  # store data into matrix
  m_groups <- as.matrix(unname(norm_group[, -1]))
  m_strat <- as.matrix(unname(norm_strat[,-1]))

  # create SE objects
  se_groups <- SummarizedExperiment::SummarizedExperiment(assays=m_groups, colData=m_data)
  se_strat <- SummarizedExperiment::SummarizedExperiment(assays=m_strat, colData=m_data)

  return(list(se_groups=se_groups, se_strat=se_strat))
}



