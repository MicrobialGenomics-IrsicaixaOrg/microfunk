library(tidyverse)

example.RPK = read_delim(system.file("extdata", "demo_genefamilies.tsv", package = "microfunk"))
example.CPM = read_delim(system.file("extdata", "demo_genefamilies-cpm.tsv", package = "microfunk"))


# Output = tibble
test_that("output is tibble", {
  expect_equal(class(read_gfam(example.RPK))[2], "tbl_df" )
  expect_equal(class(read_gfam(example.CPM))[2], "tbl_df" )
})

# 2 columns
test_that("tibble has two columns", {
  expect_equal(ncol(read_gfam(example.RPK)), 2)
  expect_equal(ncol(read_gfam(example.CPM)), 2)
})

# First column is not numeric but second column is
test_that("column values type", {
  expect_equal(is.numeric(read_gfam(example.RPK)[[1]]), FALSE)
  expect_equal(is.numeric(read_gfam(example.RPK)[[2]]), TRUE)
  expect_equal(is.numeric(read_gfam(example.CPM)[[1]]), FALSE)
  expect_equal(is.numeric(read_gfam(example.CPM)[[2]]), TRUE)
})

# Column 1 header = "# gene-family" & Column 2 name = "RPKs" or "CPM"
test_that("column names", {
  expect_equal(colnames(read_gfam(example.RPK))[1], "# gene-family")
  expect_equal(colnames(read_gfam(example.CPM))[1], "# gene-family")
  expect_equal(colnames(read_gfam(example.RPK))[2], "RPKs")
  expect_equal(colnames(read_gfam(example.CPM))[2], "CPM")
})

# If column 2 name = "CPM" -> sum = 1 million
test_that("sum of 'CPM' column equals to 1 million", {
  t = read_gfam(example.CPM)
  n = colnames(t)[2]
  if (n == "CPM"){
    filtered.t = t %>%
      filter(t[[1]] == 'UNMAPPED' |
               t[[1]] == 'UNKNOWN' |
               grepl("^UniRef90_[^|]+$", t[[1]]))

    s = sum(filtered.t[[2]], na.rm = TRUE)
    expect_true(s >= 999900 && s <= 1000100)
  } else {
    skip("Units are not CPM")
  }
})
