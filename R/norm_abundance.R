#' Normalize gene family or pathway abundance data in a SummarizedExperiment
#' object
#'
#' This function normalizes gene family or pathway abundance data stored in a
#' SummarizedExperiment object, containing gene/pathway abundance values as the
#' assay and metadata as the column data
#' (https://bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html),
#' to either CPM units (copies per million) or relative abundance.
#'
#' It first checks whether the requested normalization method is supported. If
#' so, then retrieves the abundance data from the assays slot in the SE object,
#' determines the current normalization method applied to the data (RPK, CPM or
#' relative abundance) and applies the requested normalization if needed.
#'
#' @param se A SummarizedExperiment object containing gene/pathway abundance
#'   values as the assay and metadata as the column data.
#' @param norm A character string specifying the normalization method to apply.
#'   Supported methods are 'cpm' for copies per million and 'relab' for relative
#' @param .batch_size An integer specifying the number of columns to process in each
#'   batch. Default is 10.
#'
#' @return An updated SE object with normalized abundance values in the assay
#'   slot.
#' @export
#' @autoglobal
#' @tests
#' # Sample RPK testing data
#' column1 <- c("FIRST", "FIRST|Species1", "UNGROUPED","UNGROUPED|Species1", "UNGROUPED|Species2", "UNMAPPED")
#' column2 <- c(5, 5, 100, 90, 10, 20)
#' rpk_data <- tibble::tibble(
#'  function_id = column1,
#'  sample1 = column2,
#'  sample2 = column2 ) %>%
#'  data.frame(row.names = 1)
#'
#' # Sample RPK SE object
#' rpk_se <-
#'  SummarizedExperiment::SummarizedExperiment(assays = list(humann = rpk_data))
#'
#' # Sample relative abundance testing data
#' relab_data <- rpk_data %>%
#'   dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ .x / 125))
#'
#' # Sample relative abundance SE object
#' relab_se <-
#'  SummarizedExperiment::SummarizedExperiment(assays = list(humann = relab_data))
#'
#' # Sample CPM testing data
#' cpm_data <- relab_data %>%
#'   dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ .x * 1e6))
#'
#' # Sample CPM SE object
#' cpm_se <-
#'  SummarizedExperiment::SummarizedExperiment(assays = list(humann = cpm_data))
#'
#' # Test errors
#' testthat::expect_error(norm_abundance(rpk_se, norm = "unsupported"))
#' testthat::expect_error(norm_abundance(cpm_se, norm = "cpm"))
#' testthat::expect_error(norm_abundance(relab_se, norm = "relab"))
#'
#' # Test outputs
#' inputs <- list(
#'   list(rpk_se, "cpm"),
#'   list(rpk_se, "relab"),
#'   list(cpm_se, "relab"),
#'   list(relab_se, "cpm")
#' )
#' expected_outputs <- list(
#'   cpm_se,
#'   relab_se,
#'   relab_se,
#'   cpm_se
#' )
#'
#' purrr::walk2(inputs, expected_outputs, ~ {
#'   testthat::expect_equal(norm_abundance(.x[[1]], .x[[2]]), .y)
#' })
#'
#' @examples
#' # Def data paths
#' metadata <- system.file("extdata", "ex_meta.csv", package = "microfunk")
#' file_path <-
#'   system.file("extdata", "All_genefam_cpm_kegg.tsv", package = "microfunk")
#'
#' # Read HUMAnN3 & Normalize to Relative Abundance
#' transformed_data <-
#'   read_humann(file_path, metadata) %>%
#'   norm_abundance(norm = "relab")
norm_abundance <- function(se, norm, .batch_size = 10){

  if (norm != "cpm" && norm != "relab"){
    cli::cli_abort(c(
      "x" = "The requested normalization method is not supported.",
      "i" = "Supported methods are copies per million ('cpm') and relative abundance ('relab')."
    ))
  }

  # get abundance data
  data <-
    SummarizedExperiment::assays(se)$humann %>%
    tibble::as_tibble(rownames = "function_id")

  # check normalization method
  classification <- .extract_norm_method(data)

  if (norm == classification) {
    cli::cli_abort(c(
      "x" = "The requested normalization method is already applied to the data."
    ))
  }

  # Batch Normalization
  batch_size <- 25
  batches <- split(seq(2, ncol(data)), ceiling(seq_along(seq(2, ncol(data)))/batch_size))
  norm_df <-
    seq_along(batches) %>%
    purrr::map_dfc( ~ {
      batch <- c("function_id", colnames(data)[batches[[.x]]])
      temp_df <- data[, batch, drop = FALSE]
      temp_res <- .normalize(temp_df, norm, classification)
      if (.x != 1) { temp_res <- temp_res[, -1] }
      temp_res
    }) %>% data.frame(row.names = 1)

  SummarizedExperiment::assays(se) <- list(humann = norm_df)
  se
}

#' Extract abundance data normalization method
#'
#' This function determines if abundance data for each sample is in RPK units
#' (reads per kilobase) or is normalized to CPM units (copies per million) or
#' relative abundance. It also ensures consistency of units across all samples
#' (columns).
#'
#' @param tbl A tibble containing gene/pathway abundance data.
#'
#' @return A character string indicating the normalization method ('cpm',
#'   'rpk' or 'relab').
#' @keywords internal
#' @autoglobal
#' @noRd
#' @tests
#' se <- read_humann(
#'  file_path = system.file("extdata", "cpm_kegg_error_units.tsv", package = "microfunk"),
#'  metadata = system.file("extdata", "meta_error_units.csv", package = "microfunk")
#' )
#'  data <- SummarizedExperiment::assays(se)$humann %>%
#'    tibble::as_tibble(rownames = "function_id")
#'  testthat::expect_error(.extract_norm_method(data))
.extract_norm_method <- function(tbl) {

  # check data normalization
  sum_counts <-
    dplyr::filter(tbl, !stringr::str_detect(function_id, "[|]")) %>%
    dplyr::select(-1) %>%
    colSums()

  classification <-
    dplyr::case_when(
      dplyr::between(sum_counts, 99900, 1001000) ~ "cpm",
      dplyr::between(sum_counts, 0.9, 1.1) ~ "relab",
      TRUE ~ "rpk"
    ) %>% unique()

  # Check if all are same normalization
  if (length(classification) != 1) {
    cli::cli_abort(c(
      "x" = "Columns have mixed units.",
      "i" = "Please ensure all samples have either 'rpk', 'cpm' or
      'relative abundance' units."
    ))
  }

  classification
}

#' Convert RPK units to CPM
#'
#' This function divides the data into two separate tibbles: one containing main
#' categories (gene families or pathways) and another containing the stratified
#' abundance values by species. After separating the data, it scales the values
#' within each tibble from RPK to CPM and merges normalized data into a single
#' tibble.
#'
#' @param tbl A tibble containing the gene/pathway abundance data.
#'
#' @return A tibble with the abundance data converted from RPK to CPM.
#' @keywords internal
#' @autoglobal
#' @noRd
#' @tests
#' # Input tibble
#' column1 <- c("UNMAPPED", "UNGROUPED", "UNGROUPED|Species1", "UNGROUPED|Species2", "FIRST", "FIRST|Species1")
#' column2 <- c(20, 100, 90, 10, 5, 5)
#' input_tbl <- tibble::tibble(
#'   function_id = column1,
#'   sample1 = column2,
#'   sample2 = column2
#' )
#'
#' # Expected output tibble
#' column1 <- c("FIRST", "FIRST|Species1", "UNGROUPED","UNGROUPED|Species1", "UNGROUPED|Species2", "UNMAPPED")
#' column2 <- c(40000, 40000, 800000, 720000, 80000, 160000)
#' output_tbl <- tibble::tibble(
#'   function_id = column1,
#'   sample1 = column2,
#'   sample2 = column2
#' ) %>% data.table::as.data.table()
#'
#' testthat::expect_equal(.rpk2cpm(input_tbl), output_tbl)
.rpk2cpm <- function(tbl) {
  dt <- data.table::as.data.table(tbl)
  c_num <- names(dt)[sapply(dt, is.numeric)]
  cat <- dt[!grepl("\\|", function_id)]
  strat <- dt[grepl("\\|", function_id) | function_id == "UNMAPPED"]
  cat[, (c_num) := lapply(.SD, function(x) x / sum(x) * 1e6), .SDcols = c_num]
  strat[, (c_num) := lapply(.SD, function(x) x / sum(x) * 1e6), .SDcols = c_num]
  norm_tbl <- data.table::rbindlist(list(cat, strat))
  data.table::setorder(norm_tbl, function_id)
  unique(norm_tbl, by = "function_id")
}

#' Convert relative abundance units to CPM
#'
#' This function divides the data into two separate tibbles: one containing main
#' categories (gene families or pathways) and another containing the stratified
#' abundance values by species. After separating the data, it scales the values
#' within each tibble from relative abundance to CPM and merges normalized data
#' into a single tibble.
#'
#' @param tbl A tibble containing the gene/pathway abundance data.
#'
#' @return A tibble with the abundance data converted from relative abundance to
#'   CPM.
#' @keywords internal
#' @autoglobal
#' @noRd
#' @tests
#' # Input tibble
#' column1 <- c("UNMAPPED", "UNGROUPED", "UNGROUPED|Species1", "UNGROUPED|Species2", "FIRST", "FIRST|Species1")
#' column2 <- c(0.16, 0.8, 0.72, 0.08, 0.04, 0.04)
#' input_tbl <- tibble::tibble(
#'   function_id = column1,
#'   sample1 = column2,
#'   sample2 = column2
#' )
#'
#' # Expected output tibble
#' column1 <- c("FIRST", "FIRST|Species1", "UNGROUPED","UNGROUPED|Species1", "UNGROUPED|Species2", "UNMAPPED")
#' column2 <- c(40000, 40000, 800000, 720000, 80000, 160000)
#' output_tbl <- tibble::tibble(
#'   function_id = column1,
#'   sample1 = column2,
#'   sample2 = column2
#' ) %>% data.table::as.data.table()
#'
#' testthat::expect_equal(.relab2cpm(input_tbl), output_tbl)
.relab2cpm <- function(tbl){
  dt <- data.table::as.data.table(tbl)
  c_num <- names(dt)[sapply(dt, is.numeric)]
  cat <- dt[!grepl("\\|", function_id)]
  strat <- dt[grepl("\\|", function_id) | function_id == "UNMAPPED"]
  cat[, (c_num) := lapply(.SD, function(x) x * 1e6), .SDcols = c_num]
  strat[, (c_num) := lapply(.SD, function(x) x * 1e6), .SDcols = c_num]
  norm_tbl <- data.table::rbindlist(list(cat, strat))
  data.table::setorder(norm_tbl, function_id)
  unique(norm_tbl, by = "function_id")
}


#' Convert RPK units to relative abundance
#'
#' This function divides the data into two separate tibbles: one containing main
#' categories (gene families or pathways) and another containing the stratified
#' abundance values by species. After separating the data, it scales the values
#' within each tibble from RPK to relative abundance and merges normalized data
#' into a single tibble.
#'
#' @param tbl A tibble containing the gene/pathway abundance data.
#'
#' @return A tibble with the abundance data converted from RPK to relative
#'   abundance.
#' @keywords internal
#' @autoglobal
#' @noRd
#' @tests
#' # Input tibble
#' column1 <- c("UNMAPPED", "UNGROUPED", "UNGROUPED|Species1", "UNGROUPED|Species2", "FIRST", "FIRST|Species1")
#' column2 <- c(20, 100, 90, 10, 5, 5)
#' input_tbl <- tibble::tibble(
#'   function_id = column1,
#'   sample1 = column2,
#'   sample2 = column2
#' )
#'
#' # Expected output tibble
#' column1 <- c("FIRST", "FIRST|Species1", "UNGROUPED","UNGROUPED|Species1", "UNGROUPED|Species2", "UNMAPPED")
#' column2 <- c(0.04, 0.04, 0.8, 0.72, 0.08, 0.16)
#' output_tbl <- tibble::tibble(
#'   function_id = column1,
#'   sample1 = column2,
#'   sample2 = column2
#' ) %>% data.table::as.data.table()
#'
#' testthat::expect_equal(.rpk2relab(input_tbl), output_tbl)
.rpk2relab <- function(tbl){
  dt <- data.table::as.data.table(tbl)
  c_num <- names(dt)[sapply(dt, is.numeric)]
  cat <- dt[!grepl("\\|", function_id)]
  strat <- dt[grepl("\\|", function_id) | function_id == "UNMAPPED"]
  cat[, (c_num) := lapply(.SD, function(x) x / sum(x)), .SDcols = c_num]
  strat[, (c_num) := lapply(.SD, function(x) x / sum(x)), .SDcols = c_num]
  norm_tbl <- data.table::rbindlist(list(cat, strat))
  data.table::setorder(norm_tbl, function_id)
  unique(norm_tbl, by = "function_id")
}

#' Convert CPM units to relative abundance
#'
#' This function divides the data into two separate tibbles: one containing main
#' categories (gene families or pathways) and another containing the stratified
#' abundance values by species. After separating the data, it scales the values
#' within each tibble from CPM to relative abundance and merges normalized data
#' into a single tibble.
#'
#' @param tbl A tibble containing the gene/pathway abundance data.
#'
#' @return A tibble with the abundance data converted from CPM to relative
#'   abundance.
#' @keywords internal
#' @autoglobal
#' @noRd
#' @tests
#' # Input tibble
#' column1 <- c("UNMAPPED", "UNGROUPED", "UNGROUPED|Species1", "UNGROUPED|Species2", "FIRST", "FIRST|Species1")
#' column2 <- c(160000, 800000, 720000, 80000, 40000, 40000)
#' input_tbl <- tibble::tibble(
#'   function_id = column1,
#'   sample1 = column2,
#'   sample2 = column2
#' )
#' # Expected output tibble
#' column1 <- c("FIRST", "FIRST|Species1", "UNGROUPED","UNGROUPED|Species1", "UNGROUPED|Species2", "UNMAPPED")
#' column2 <- c(0.04, 0.04, 0.8, 0.72, 0.08, 0.16)
#' output_tbl <- tibble::tibble(
#'   function_id = column1,
#'   sample1 = column2,
#'   sample2 = column2
#' ) %>% data.table::as.data.table()
#'
#' testthat::expect_equal(.cpm2relab(input_tbl), output_tbl)
.cpm2relab <- function(tbl){
  dt <- data.table::as.data.table(tbl)
  c_num <- names(dt)[sapply(dt, is.numeric)]
  cat <- dt[!grepl("\\|", function_id)]
  strat <- dt[grepl("\\|", function_id) | function_id == "UNMAPPED"]
  cat[, (c_num) := lapply(.SD, function(x) x / 1e6), .SDcols = c_num]
  strat[, (c_num) := lapply(.SD, function(x) x / 1e6), .SDcols = c_num]
  norm_tbl <- data.table::rbindlist(list(cat, strat))
  data.table::setorder(norm_tbl, function_id)
  unique(norm_tbl, by = "function_id")
}

#' Normalize data based on specified method and classification
#'
#' This function normalizes a given data table based on the specified
#' normalization method and classification. It supports converting between CPM
#' and RPK, and between CPM and relative abundance (relab).
#'
#' @param tbl A data frame containing the data to be normalized.
#' @param norm A character string specifying the normalization method. Accepts
#'   "cpm" or "relab".
#' @param classification A character string specifying the classification of the
#'   data. Accepts "rpk" or "relab".
#' @return A data frame with normalized data.
#' @keywords internal
#' @autoglobal
#' @noRd
.normalize <- function(tbl, norm, classification) {
  if (norm == "cpm" & classification == "rpk") { res <- .rpk2cpm(tbl) }
  if (norm == "cpm" & classification == "relab") { res <- .relab2cpm(tbl) }
  if (norm == "relab" & classification == "rpk") { res <- .rpk2relab(tbl) }
  if (norm == "relab" & classification == "cpm") { res <- .cpm2relab(tbl) }
  res
}
