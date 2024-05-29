# List all files in the specified directory that match the pattern ".txt.gz"
all_files <- list.files(
  "~/humann_dbs/full_mapping_v201901/",
  pattern = ".txt.gz",
  full.names = TRUE
)

# Set environment variable to increase the connection size for vroom
Sys.setenv("VROOM_CONNECTION_SIZE" = "500000000")

# Define the list of database types to process
c("ec", "eggnog", "go", "ko", "pfam", "test") %>%
  purrr::walk(~ {
    test <- FALSE
    if (.x == "test") {
      .x = "ec"
      test <- TRUE
    }

    # Load the names associated with each database type
    names <-
      all_files[stringr::str_detect(all_files, paste0(.x, "_name.txt.gz"))] %>%
      data.table::fread(col.names = c("id", "id_name"))

    # Load the UniRef90 IDs associated with each database type
    ids <-
      all_files[stringr::str_detect(all_files, paste0(.x, "_uniref90.txt.gz"))] %>%
      vroom::vroom(delim = ",", col_names = "uniref_90", show_col_types = FALSE) %>%
      dplyr::mutate(
        id = stringr::str_remove_all(uniref_90, "\t.*"),
        uniref_90 = stringr::str_remove_all(uniref_90, paste0(id,"\t")),
        uniref_90 = stringr::str_remove_all(uniref_90, "UniRef90_"),
        .before = 1
      ) %>%
      tidyr::separate_longer_delim(cols = uniref_90, delim = "\t")

    if (test) {
      ids <- dplyr::filter(ids, id == "1.1.1.10")
      .x <- "test"
    }

    # Join the IDs with their associated names and write to S3
    dplyr::inner_join(ids, names, by = "id") %>%
      aws.s3::s3write_using(
        FUN = data.table::fwrite,
        bucket = "bioconductor-packages",
        object = paste0("microfunk/humann_db/", .x, "_v201901b.tsv.gz"),
        opts = list(multipart = TRUE)
      )
  })
