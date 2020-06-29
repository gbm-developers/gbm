skip_on_cran()
skip_on_appveyor()

test_that("Show method for ANN class produces the correct output", {
  data <- SDMtune:::t
  data@data <- data@data[, 1:4]
  m <- train("ANN", data = data, size = 10)@model
  expect_output(print(m), "Class: ANN", fixed = TRUE)
  expect_output(print(m), "size : 10", fixed = TRUE)
  expect_output(print(m), "decay: 0", fixed = TRUE)
  expect_output(print(m), "rang : 0.7", fixed = TRUE)
  expect_output(print(m), "maxit: 100", fixed = TRUE)
})
