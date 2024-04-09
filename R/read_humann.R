#' Read HUMAnN Gene Families
#'
#' This function reads a HUMAnN gene family file (https://huttenhower.sph.harvard.edu/humann),
#' which contains information about the abundance of gene families of a given sample in a .tsv format.
#' The function performs several data validation checks, distinguishing between normalized
#' (CPM, copies per million)and non-normalized abundance data (RPKs, reads per kilobase).
#' It returns a modified tibble with renamed columns for easier downstream analysis.
#'
#' @param file_path A string. The path to the HUMAnN gene families file.
#'
#' @return A tibble of the gene families data with  three columns: 'class' (file type), 'annotation' (gene family)
#'   and 'norm' (abundance units).
#' @export
#' @tests
#' files <- system.file("extdata", package = "microfunk") %>%
#' list.files(pattern = "demo_genefamilies", full.names = TRUE)
#'
#' error_files <-
#' which(stringr::str_detect(files, "error")) %>%
#' files[.]
#'
#' norm_file <-
#' which(stringr::str_detect(files, "cpm")) %>%
#' files[.]
#'
#' purrr::map(error_files, ~ testthat::expect_error(read_genefamily(.)))
#'
#' testthat::expect_equal(colnames(read_genefamily(norm_file))[3], "cpm")
#'
#' @examples
#' file_path <-
#'   system.file("extdata", "demo_genefamilies.tsv", package = "microfunk")
#'
#' gfam <- read_genefamily(file_path)

read_genefamily <- function(file_path) {
  # read in the file
  tbl <-
    data.table::fread(file_path, sep="\t") %>%
    tibble::as_tibble()

  # check header
  header <- colnames(tbl)
  if (!stringr::str_detect(header[1], "^# Gene Family$")) {
    cli::cli_abort("The header of the file should be '# Gene Family'.")
  }

  # check number of rows
  if (nrow(tbl) <= 1) {
    cli::cli_abort("The file does not contain any data. Ensure the file is not empty.")
  }

  # check number of columns
  if (ncol(tbl) != 2) { cli::cli_abort("ncol(input) != 2. The file should contain two columns.") }

  # check column classes
  class_1 <- class(tbl[[1]]) == "character"
  class_2 <- class(tbl[[2]]) == "numeric"
  if (!class_1 || !class_2) {
    cli::cli_abort(c("x" = "Column classes are not correct.",
                   "i" = "First column should be 'character' and second one should be 'numeric'."))
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



#' Read HUMAnN Pathway Abundances
#'
#' This function reads a HUMAnN pathway abundance file (https://huttenhower.sph.harvard.edu/humann),
#' which contains information about the abundance of metabolic pathways present in a given sample
#' in a .tsv format. The function performs several data validation checks, distinguishing between normalized
#' (CPM, copies per million)and non-normalized abundance data (RPKs, reads per kilobase).
#' It returns a modified tibble with renamed columns for easier downstream analysis.
#'
#' @param file_path A string. The path to the HUMAnN pathway abundance file.
#'
#' @return A tibble of the pathway abundance data with  three columns: 'class' (file type), 'annotation' (pathway)
#'   and 'norm' (abundance units).
#' @export
#' @tests
#' files <- system.file("extdata", package = "microfunk") %>%
#' list.files(pattern = "demo_pathabundance", full.names = TRUE)
#'
#' error_files <-
#' which(stringr::str_detect(files, "error")) %>%
#' files[.]
#'
#' norm_file <-
#' which(stringr::str_detect(files, "cpm")) %>%
#' files[.]
#'
#' purrr::map(error_files, ~ testthat::expect_error(read_pathabundance(.)))
#'
#' testthat::expect_equal(colnames(read_pathabundance(norm_file))[3], "cpm")
#'
#' @examples
#' file_path <-
#'   system.file("extdata", "demo_pathabundance.tsv", package = "microfunk")
#'
#' pthab <- read_pathabundance(file_path)
read_pathabundance <- function(file_path) {
  # read in the file
  tbl <-
    data.table::fread(file_path, sep="\t") %>%
    tibble::as_tibble()

  # check header
  header <- colnames(tbl)
  if (!stringr::str_detect(header[1], "^# Pathway$")) {
    cli::cli_abort("The header of the file should be '# Pathway'.")
  }

  # check number of rows
  if (nrow(tbl) <= 1) {
    cli::cli_abort("The file does not contain any data. Ensure the file is not empty.")
  }

  # check number of columns
  if (ncol(tbl) != 2) { cli::cli_abort("ncol(input) != 2. The file should contain two columns.") }

  # check column classes
  class_1 <- class(tbl[[1]]) == "character"
  class_2 <- class(tbl[[2]]) == "numeric"
  if (!class_1 || !class_2) {
    cli::cli_abort(c("x" = "Column classes are not correct.",
                     "i" = "First column should be 'character' and second one should be 'numeric'."))
  }

  # check data normalitzation
  sum_counts <-
    tbl %>%
    dplyr::filter(!stringr::str_detect(`# Pathway`, "[|]")) %>%
    dplyr::pull(2) %>%
    sum()

  norm <- "rpk"
  if (dplyr::between(sum_counts, 999900, 1000100)) { norm <- "cpm" }

  # Prepare tbl
  tbl <-
    tbl %>%
    stats::setNames(c("annotation", norm)) %>%
    dplyr::mutate(class = "pathway", .before = 1)

  tbl
}
