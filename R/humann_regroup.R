#' Regroup HUMAnN3 Abundance Data in a SummarizedExperiment Object
#'
#' This function regroups gene family or metabolic pathway abundance data stored
#' in a SummarizedExperiment object and annotated with UniRef90 IDs
#' (https://www.uniprot.org/help/uniref).
#'
#' It takes as inputs a SummarizedExperiment object, the new annotation method
#' to regroup to and the version of the database to use. First, the validity of
#' inputs is checked using the `.check_regroup()` function. Then, it processes
#' the abundance data, extracting its relevant UniRef90 IDs, and ensures that
#' the necessary database file for the annotation conversion are available,
#' fetching them if required. The function reads those files and matches them
#' with abundance data UniRef90 annotation. Finally, it merges and aggregates
#' the data based on the new annotations and constructs a new
#' SummarizedExperiment object.
#'
#' @param se A SummarizedExperiment object containing gene/pathway abundance
#'   values as the assay and metadata as the column data.
#' @param to A character string specifying the target annotation type for
#'   regrouping. Valid values include "ec" (Enzyme Comission), "eggnog"
#'   (evolutionary genealogy of genes: Non-supervised Orthologous Groups), "go"
#'   (Gene Ontology), "ko" (KEGG Orthology) and "pfam" (Protein families
#'   database of alignments and HMMs).

#' @param version A character string specifying the version of the database to
#'   use. Default is "v201901b".
#'
#' @return A SummarizedExperiment object containing gene/pathway abundance
#'   values regrouped according to the specified annotation as the assay,
#'   metadata as the column data and information about the regrouped functions
#'   as the row data.
#' @export
#' @autoglobal
#' @tests
#' se_regrouped <- read_humann(
#'  file_path = system.file("extdata", "reduced_genefam_rpk_uniref.tsv", package = "microfunk"),
#'  metadata = system.file("extdata", "reduced_meta.csv", package = "microfunk")
#'  ) %>%
#'  humann_regroup(to = "test")
#'
#' # Test output
#' testthat::expect_equal(class(se_regrouped)[1],
#'   "SummarizedExperiment")
#'
#' @examples
#' # Read HUMAnN3 & Regroup to PFAM annotation
#' se_regrouped <- read_humann(
#'   file_path = system.file("extdata", "reduced_genefam_rpk_uniref.tsv", package = "microfunk"),
#'   metadata = system.file("extdata", "reduced_meta.csv", package = "microfunk")
#'   ) %>%
#'   humann_regroup(to = "pfam")
humann_regroup <- function(se, to, version = "v201901b") {

  # Checks
  .check_regroup(se, to)

  dt <-
    SummarizedExperiment::assays(se)$humann %>%
    data.table::setDT(keep.rownames = "function_id")


  dt[, uniref_90 := stringr::str_remove_all(function_id, "[|].*|UniRef90_")]
  data.table::setcolorder(dt, "uniref_90", after = "function_id")
  ids_retain <- unique(dt$uniref_90) %>% paste0(collapse = "|")

  db_exists <- .make_humann_name(to) %>% .file_db_exists()
  if (!db_exists) { fetch_humann_db(annot = to) }

  ids <- data.table::fread(file.path(.cache_path(), .make_humann_name(to))) # Optimizar
  ids <- ids[uniref_90 %in% dt$uniref_90,]

  dt <- dt[ids, on = "uniref_90"]
  data.table::setcolorder(dt, c("id", "id_name"), after = "function_id")
  dt[, function_id := stringr::str_replace_all(function_id, paste0("UniRef90_", uniref_90), id)]
  dt <- dplyr::select(dt, -uniref_90, -id)
  dt <- dt[, lapply(.SD, sum), .SDcols = 3:ncol(dt), by = .(function_id, id_name)]

  se_dt <- dplyr::select(dt, -id_name) %>% data.frame(row.names = 1)
  SummarizedExperiment::SummarizedExperiment(
    assays = list(humann = se_dt),
    rowData = dt[, .(function_id, id_name)] %>% data.frame(row.names = 1),
    colData = SummarizedExperiment::colData(se)
  )
}

#' Check annotation regrouping validity for abundance data in a SE object
#'
#' This function checks the validity of regrouping the gene family/pathway
#' abundance data to a requested target annotation. It ensures that the IDs in
#' the abundance assay data from the SE object are known and belong to a valid
#' annotation type. Then, it checks that only UniRef90 IDs are present for
#' regrouping. Finally, it validates that the target annotation type is among
#' the supported annotations.
#'
#' If any of these checks fail, an appropriate error message is displayed, and
#' the function aborts the regrouping process.
#'
#' @param se A SummarizedExperiment object containing gene/pathway abundance
#'   values as the assay and metadata as the column data.
#' @param to A character string specifying the target annotation type for
#'   regrouping. Valid values include "ec" (Enzyme Comission), "eggnog"
#'   (evolutionary genealogy of genes: Non-supervised Orthologous Groups), "go"
#'   (Gene Ontology), "ko" (KEGG Orthology) and "pfam" (Protein families
#'   database of alignments and HMMs).
#'
#' @return 'NULL' if all checks pass. If any issue is found, an error is raised.
#' @keywords internal
#' @autoglobal
#' @tests
#' # Test SE annotation errors
#' meta <-  system.file("extdata", "reduced_meta.csv", package = "microfunk")
#' se_list <- system.file("extdata", package = "microfunk") %>%
#'  list.files(pattern = "reduced_genefam_rpk_uniref_error", full.names = TRUE) %>%
#'    purrr::map(~ read_humann(.x, meta))
#'  purrr::walk(se_list, ~ testthat::expect_error(.check_regroup(.x, to = "test")))
#'
#' # Test invalid annotation requested
#' se <- read_humann(
#'  file_path = system.file("extdata", "reduced_genefam_rpk_uniref.tsv", package = "microfunk"),
#'  metadata = meta
#'  )
#'  testthat::expect_error(.check_regroup(se, to = "abc"))
.check_regroup <- function(se, to) {
  ids <- SummarizedExperiment::assays(se)$humann %>% rownames() %>% .[. != "UNMAPPED"]
  annot <- dplyr::case_when(
    stringr::str_detect(ids, "UniRef90_") ~ "uniref90",
    stringr::str_detect(ids, "^K") ~ "kegg",
    stringr::str_detect(ids, "^PF") ~ "pfam",
    stringr::str_detect(ids, "^COG|^ENO") ~ "eggnog",
    stringr::str_detect(ids, "^[0-9][.][0-9]") ~ "ec",
    stringr::str_detect(ids, "^GO") ~ "go",
    TRUE ~ "unknown"
  ) %>% unique()

  if (length(annot) > 1) {
    cli::cli_alert_warning(
      "Multiple annotation types detected: {cli::col_blue(annot)}"
    )
    cli::cli_li("Skipping {.code humann_regroup}")
  }

  if (annot == "unknown") {
    cli::cli_abort(c("x" = "No known annotation types detected"))
  }

  if (annot != "uniref90") {
    cli::cli_abort(c("x" = "Only UniRef90 IDs are supported for regrouping"))
  }

  available <- c("ec", "eggnog", "go", "ko", "pfam", "test")
  if (!all(to %in% available)) {
    cli::cli_abort(c(
      "x" = "Invalid regrouping annotation requested",
      "i" = "Available annotations: {cli::col_blue(available)}"
    ))
  }

}





