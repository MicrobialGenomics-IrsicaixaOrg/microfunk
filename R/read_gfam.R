#' Read HUMAnN Gene Families
#'
#' This function reads in a HUMAnN gene families file and performs several
#' checks on the data. It also normalizes the data and renames the columns for
#' easier downstream analysis.
#'
#' @param file_path A string. The path to the HUMAnN gene families file.
#'
#' @return A tibble. The gene families data with columns: 'class', 'annotation',
#'   and 'norm' (normalized counts).
#' @export
#' @tests
#' file_path <-
#'   system.file("extdata", "demo_genefamilies.tsv", package = "microfunk")
#'
#' # Test output has tree columns
#' testthat::expect_equal(ncol(read_gfam(file_path)), 3)
#'
#' @examples
#' file_path <-
#'   system.file("extdata", "demo_genefamilies.tsv", package = "microfunk")
#'
#' gfam <- read_gfam(file_path)
read_gfam <- function(file_path) {
  # read in the file
  tbl <-
    data.table::fread(file_path) %>%
    tibble::as_tibble()

  # check header
  header <- colnames(tbl)
  if (!stringr::str_detect(header[1], "^# Gene Family$")) {
    cli::cli_abort("The header of the file is not correct.")
  }

  # check number of rows
  if (nrow(tbl) <= 1) {
    cli::cli_abort("The file does not contain any data.")
  }

  # check number of columns
  if (ncol(tbl) != 2) { cli::cli_abort("ncol(input) != 2") }

  # check column classes
  class_1 <- class(tbl[[1]]) == "character"
  class_2 <- class(tbl[[2]]) == "numeric"
  if (!class_1 || !class_2) {
    cli::cli_abort("The column classes are not correct.")
  }

  # check data normalitzation
  sum_counts <-
    tbl %>%
    dplyr::filter(!stringr::str_detect(`# Gene Family`, "[|]")) %>%
    dplyr::pull(2) %>%
    sum()

  norm <- "rpk"
  if (dplyr::between(sum_counts, 999900, 1000100)) { norm <- "cpm" }

  # Prepare tbl
  tbl <-
    tbl %>%
    stats::setNames(c("annotation", norm)) %>%
    dplyr::mutate(class = "gene_family", .before = 1)

  tbl
}
