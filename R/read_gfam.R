#' Title
#'
#' @param input A tsv file of gene family abundances.
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' example.file = system.file("extdata", "demo_genefamilies.tsv",
#' package = "microfunk")
#' read_gfam(example.file)
#'

library(tidyverse)

read_gfam = function(input){

  # Check header
  header = input %>%
    select(1) %>%
    colnames()

  if (!str_detect(header, "^# Gene Family$")){
    stop("Header does not match the expected format.")
  }

  # Rename header
  data = input %>%
    rename("# gene-family" = header)

  # Check rows
  if (nrow(input) == 0){
    stop("No data rows found.")
  }

  # Check columns
  if (ncol(input) != 2){
    stop ("Incorrect number of columns.")
  } else{
    name = input %>%
      select(2) %>%
      colnames()

    if (!is.numeric(input[[2]])){
      stop ("Second column is not numeric.")
    } else{

      # Check units (RPKs/CPM)
      filtered.data = data %>%
        filter(data[[1]] == 'UNMAPPED' |
                 data[[1]] == 'UNKNOWN' |
                 grepl("^UniRef90_[^|]+$", data[[1]]))

      sum_col = sum(filtered.data[[2]], na.rm = TRUE)

      if (sum_col >= 999900 && sum_col <= 1000100){
        # CPM
        data = data %>%
          rename("CPM" = name)
      } else{
        # RPKs
        data = data %>%
          rename("RPKs" = name)
      }
    }
  }
  # Return transformed tibble
  return (data)

}












