
a <- function() {
  se <- read_humann(
    file_path = "~/Downloads/All_genefam.tsv",
    metadata = system.file("extdata", "ex_meta.csv", package = "microfunk")
  ) %>% norm_abundance("cpm")

  to = "pfam"
}



humann_regroup <- function(se, to = "pfam") {

  .check_regroup()

  dt <-
    SummarizedExperiment::assays(se)$humann %>%
    data.table::setDT(keep.rownames = "function_id")


  dt[, uniref_90 := stringr::str_remove_all(function_id, "[|].*|UniRef90_")]
  data.table::setcolorder(dt, "uniref_90", after = "function_id")
  ids_retain <- unique(dt$uniref_90) %>% paste0(collapse = "|")

  db_exists <- .make_humann_name(to) %>% .file_db_exists()
  if (!db_exists) { fetch_humann_db(annot = to) }

  ids <- data.table::fread(file.path(.cache_path(), .make_humann_name(to)))
  ids <- ids[uniref_90 %in% dt$uniref_90,]

  dt <- dt[ids, on = "uniref_90"]
  data.table::setcolorder(dt, c("id", "id_name"), after = "function_id")
  dt[, function_id := stringr::str_replace_all(function_id, paste0("UniRef90_", uniref_90), id)]
  dt <- dplyr::select(dt, -uniref_90, -id)
  dt <- dt[, lapply(.SD, sum), .SDcols = 4:ncol(dt), by = .(function_id, id_name)]

  se_dt <- dplyr::select(dt, -id_name) %>% data.frame(row.names = 1)
  SummarizedExperiment::SummarizedExperiment(
    assays = list(humann = se_dt),
    rowData = dt[, .(function_id, id_name)] %>% data.frame(row.names = 1),
    colData = SummarizedExperiment::colData(se)
  )
}

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

}





