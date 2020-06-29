skip_on_cran()

files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)
set.seed(25)
x <- dismo::randomPoints(predictors, 9000)

test_that("The function remove coords where there are NA (ex with matrix)", {
  c <- thinData(x, predictors)
  expect_true(nrow(c) < 9000)
  expect_true(inherits(c, "matrix"))
  expect_equal(colnames(c), colnames(x))
})

test_that("The function remove duplicated data (ex with dataframe)", {
  c <- thinData(as.data.frame(rbind(x, x)), predictors)
  expect_true(nrow(c) < 9000)
  expect_true(inherits(c, "data.frame"))
  expect_equal(colnames(c), colnames(x))
})

test_that("The function works with custom dataframe", {
  df <- data.frame(A = x[, "x"], B = x[, "y"], t = rep("a", nrow(x)))
  c <- thinData(df, predictors, x = "A", y = "B")
  expect_true(nrow(c) < 9000)
  expect_true(inherits(c, "data.frame"))
  expect_equal(colnames(c), colnames(df))
})

test_that("The function raises errors", {
  colnames(x) <- c("A", "B")
  expect_error(thinData(x, predictors), "The column 'x' is not present")
  colnames(x) <- c("x", "B")
  expect_error(thinData(x, predictors), "The column 'y' is not present")
})
