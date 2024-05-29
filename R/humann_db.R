#' Get the path to the cache directory
#'
#' This function returns the path to the cache directory used by microfunk.
#'
#' @return A character string representing the path to the cache directory.
#'
#' @keywords internal
#' @noRd
.cache_path <- function() {
  file.path("~", ".microfunk", "humann_db")
}

#' Create a human-readable name for a database file
#'
#' This function generates a human-readable name for a database file based on
#' its annotation and version.
#'
#' @param annot A character string representing the annotation.
#' @param version A character string representing the version.
#'
#' @return A character string representing the human-readable name of the
#'   database file.
#'
#' @keywords internal
#' @noRd
.make_humann_name <- function(annot, version) {
  paste0(annot, "_", version, ".tsv.gz")
}

#' Check if a given version is valid
#'
#' This function checks if a provided version is valid.
#'
#' @param version A character string representing the version to be checked.
#'
#' @keywords internal
#' @noRd
.check_version <- function(version) {
  if (!version %in% c("v201901b")) {
    cli::cli_abort(c(
      "x" = "Invalid db version",
      "i" = "Available versions: {cli::col_blue('v201901b')}"
    ))
  }
}

#' Create the microfunk cache directory
#'
#' This function creates the cache directory for microfunk if it does not
#' already exist.
#'
#' @keywords internal
#' @noRd
.creat_microfunk_chache <- function() {
  if (dir.exists(.cache_path())) {
    cli::cli_alert_info("Cache directory {cli::col_blue('.microfunk')} already exists")
    return(invisible(TRUE))
  }

  dir.create(.cache_path(), recursive = TRUE)
  cli::cli_alert_success("Created {cli::col_blue('.microfunk')} directory")
}

#' Check if a file exists in the database directory
#'
#' This function checks if a file exists in the database directory within the
#' cache directory.
#'
#' @param file_name A character string representing the name of the file.
#'
#' @return TRUE if the file exists, FALSE otherwise.
#'
#' @keywords internal
#' @noRd
.file_db_exists <- function(file_name) {
  exists <- file.exists(file.path(.cache_path(), file_name))
  if (exists) {
    cli::cli_alert_success("File {cli::col_blue(file_name)} cached")
  } else {
    cli::cli_alert_info("File {cli::col_blue(file_name)} not cached")
  }
  !exists
}

#' Download a file from AWS
#'
#' This function downloads a file from AWS and saves it to the cache directory.
#'
#' @param file_name A character string representing the name of the file to be
#'   downloaded.
#'
#' @keywords internal
#' @noRd
.get_from_aws <- function(file_name, mode = "wb", quiet = TRUE) {
  dest_file <- file.path(.cache_path(), file_name)
  b_url <- "https://bioconductor-packages.s3.amazonaws.com/microfunk/humann_db/"
  url <- paste0(b_url, file_name)

  cli::cli_alert_info("Fetching {cli::col_blue(file_name)} from AWS")
  h <- curl::new_handle()
  curl::handle_setheaders(h, `accept-encoding` = "gz")
  timeout_seconds = 1800
  curl::handle_setopt(h, timeout_ms = timeout_seconds * 1000)
  result = tryCatch({
    curl::curl_download(url, dest_file, mode = mode, quiet = quiet, handle = h)
    return(TRUE)
  }, error = function(e) {
    cli::cli_abort(c(
      "x" = "Error downloading {cli::col_blue(file_name)}", "!" = e$message
    ))
    if (file.exists(dest_file)) { file.remove(dest_file) }
    return(FALSE)
  })
  cli::cli_alert_success(
    "File {cli::col_blue(file_name)} downloaded to {cli::col_blue(.cache_path())}"
  )
}

#' Fetch HUMAnN3 databases
#'
#' This function fetches HUMAnN3 databases from AWS based on the provided
#' annotation and version. If the databases do not exist in the cache directory,
#' or if the overwrite option is enabled, the function will download the
#' databases from AWS.
#'
#' @param annot A character vector specifying the annotations to fetch. Default
#'   is c("ec", "eggnog", "go", "ko", "pfam").
#' @param version A character string representing the version of the databases
#'   to fetch. Default is "v201901b".
#' @param overwrite A logical indicating whether to overwrite existing files in
#'   the cache directory. Default is FALSE.
#'
#' @details This function first checks if the provided annotation is valid. If
#'   not, it aborts with a message indicating the available annotations. Then,
#'   it creates the microfunk cache directory if it does not exist, and checks
#'   the version compatibility. Finally, it fetches the databases from AWS and
#'   saves them to the cache directory. If the overwrite option is TRUE,
#'   existing files in the cache directory will be overwritten.
#'
#' @return invisible TRUE
#' @export
#' @tests
#' # Remove existen file
#' file.remove(file.path("~/.microfunk/humann_db", "test_v201901b.tsv.gz"))
#'
#' # Normal tests
#' testthat::expect_invisible(fetch_humann_db("test"))
#' testthat::expect_invisible(fetch_humann_db("test"))
#' testthat::expect_error(fetch_humann_db("no_db"))
#' testthat::expect_error(fetch_humann_db(version = "v201901a"))
#' testthat::expect_invisible(fetch_humann_db("test", overwrite = TRUE))
#'
#' # Snapshot tests
#' file.remove(file.path("~/.microfunk/humann_db", "test_v201901b.tsv.gz"))
#' testthat::expect_snapshot(fetch_humann_db("test"))
#' testthat::expect_snapshot(fetch_humann_db("test"))
#' testthat::expect_snapshot(fetch_humann_db("test", overwrite = TRUE))
#'
#' @examples
#' # Fetch HUMAnN3 databases with default parameters
#' fetch_humann_db(annot = "test")
fetch_humann_db <- function(annot = c("ec", "eggnog", "go", "ko", "pfam"),
                            version = "v201901b",
                            overwrite = FALSE) {

  available <- c("ec", "eggnog", "go", "ko", "pfam", "test")
  if (!all(annot %in% available)) {
    cli::cli_abort(c(
      "x" = "Invalid db annotation",
      "i" = "Available annotations: {cli::col_blue(available)}"
    ))
  }

  cli::cli_h1("Fetching HUMAnN3 databases")
  .check_version(version)
  .creat_microfunk_chache()

  .make_humann_name(annot, version) %>%
    purrr::walk( ~ {
      need_download <- .file_db_exists(.x)
      if (need_download | overwrite) { .get_from_aws(.x) }
    })

  invisible(TRUE)
}
