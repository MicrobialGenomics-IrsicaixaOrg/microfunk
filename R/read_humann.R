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
#' file_path1 <-
#'   system.file("extdata", "demo_genefamilies.tsv", package = "microfunk")
#'
#' file_path2 <-
#'   system.file("extdata", "demo_genefamilies_error_header.tsv", package = "microfunk")
#'
#' file_path3 <-
#'   system.file("extdata", "demo_genefamilies_error_nrow.tsv", package = "microfunk")
#'
#' file_path4 <-
#'   system.file("extdata", "demo_genefamilies_error_ncol.tsv", package = "microfunk")
#'
#' file_path5 <-
#'   system.file("extdata", "demo_genefamilies_error_colclass.tsv", package = "microfunk")
#'
#' file_path6 <-
#'   system.file("extdata", "demo_genefamilies-cpm.tsv", package = "microfunk")
#'
#'
#' # Test output has three columns
#' testthat::expect_equal(ncol(read_genefamily(file_path1)), 3)
#'
#' # Test header
#' testthat::expect_error(read_genefamily(file_path2))
#'
#' # Test number of rows
#' testthat::expect_error(read_genefamily(file_path3))
#'
#' # Test number of columns
#' testthat::expect_error(read_genefamily(file_path4))
#'
#' # Test column class
#' testthat::expect_error(read_genefamily(file_path5))
#'
#' # Test normalization
#' testthat::expect_equal(colnames(read_genefamily(file_path6))[3], "cpm")
#'
#' @examples
#' file_path <-
#'   system.file("extdata", "demo_genefamilies.tsv", package = "microfunk")
#'
#' gfam <- read_genefamily(file_path)

read_genefamily <- function(file_path) {
  # read in the file
  tbl <-
    data.table::fread(file_path) %>%
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
    cli::cli_abort("Column classes are not correct. The first column should contain strings
                   ('character') and the second column should contain numbers ('numeric').")
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
#' file_path1 <-
#'   system.file("extdata", "demo_pathabundance.tsv", package = "microfunk")
#'
#' file_path2 <-
#'   system.file("extdata", "demo_pathabundance_error_header.tsv", package = "microfunk")
#'
#' file_path3 <-
#'   system.file("extdata", "demo_pathabundance_error_nrow.tsv", package = "microfunk")
#'
#' file_path4 <-
#'   system.file("extdata", "demo_pathabundance_error_ncol.tsv", package = "microfunk")
#'
#' file_path5 <-
#'   system.file("extdata", "demo_pathabundance_error_colclass.tsv", package = "microfunk")
#'
#' file_path6 <-
#'   system.file("extdata", "demo_pathabundance-cpm.tsv", package = "microfunk")
#'
#'
#' # Test output has three columns
#' testthat::expect_equal(ncol(read_pathabundance(file_path1)), 3)
#'
#' # Test header
#' testthat::expect_error(read_pathabundance(file_path2))
#'
#' # Test number of rows
#' testthat::expect_error(read_pathabundance(file_path3))
#'
#' # Test number of columns
#' testthat::expect_error(read_pathabundance(file_path4))
#'
#' # Test column class
#' testthat::expect_error(read_pathabundance(file_path5))
#'
#' # Test normalization
#' testthat::expect_equal(colnames(read_pathabundance(file_path6))[3], "cpm")
#'
#' @examples
#' file_path <-
#'   system.file("extdata", "demo_pathabundance.tsv", package = "microfunk")
#'
#' pthab <- read_pathabundance(file_path)
read_pathabundance <- function(file_path) {
  # read in the file
  tbl <-
    data.table::fread(file_path) %>%
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
    cli::cli_abort("Column classes are not correct. The first column should contain strings
                   ('character') and the second column should contain numbers ('numeric').")
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
