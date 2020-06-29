file1 <- tempfile(fileext = ".csv")
file2 <- tempfile(fileext = ".csv")
x <- SDMtune:::t

test_that("Function saves object in one file correctly", {
  swd2csv(x, file1)
  f1 <- read.csv(file1)
  expect_true(file.exists(file1))
  expect_named(f1, c("Species", "pa", names(x@coords), names(x@data)))
})

test_that("Function saves object in two files correctly", {
  expect_silent(swd2csv(x, c(file1, file2)))
  # Presence file
  f1 <- read.csv(file1)
  expect_true(file.exists(file1))
  expect_named(f1, c("Species", names(x@coords), names(x@data)))
  # Absence/background file
  f2 <- read.csv(file2)
  expect_true(file.exists(file2))
  expect_named(f2, c("Species", names(x@coords), names(x@data)))
})

teardown(unlink(c(file1, file2)))
