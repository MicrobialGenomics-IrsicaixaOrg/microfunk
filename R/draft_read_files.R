read_humann <- function(dir, gene_pattern, path_pattern) {
  # check if directory exists
  if (!dir.exists(dir)) {
    cli::cli_abort("Directory does not exist")
  }

  # list files in the directory
  files <- list.files(dir, full.names = TRUE)

  # filter files
  gene_files <- files[stringr::str_detect(files, gene_pattern)]
  path_files <- files[stringr::str_detect(files, path_pattern)]

  if(length(gene_files) == 0 && length(path_files) == 0){
    cli::cli_abort("No files matching the provided patterns are found")
  }

  # read gene family and pathway abundance files
  gene_data <- purrr::map_dfr(gene_files, read_genefamily)
  path_data <- purrr::map_dfr(path_files, read_pathabundance)

  # merge into a single tibble
  data <- bind_rows(gene_data, path_data)

  data

}


