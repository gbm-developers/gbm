skip_on_cran()
skip_on_appveyor()

test_that("Show method for BRT class produces the correct output", {
  data <- SDMtune:::t
  data@data <- data@data[, 1:4]
  m <- trainBRT(data = data, n.trees = 200, shrinkage = 0.2,
                bag.fraction = 0.6)@model
  expect_output(print(m), "Class            : BRT", fixed = TRUE)
  expect_output(print(m), "distribution     : bernoulli", fixed = TRUE)
  expect_output(print(m), "n.trees          : 200", fixed = TRUE)
  expect_output(print(m), "shrinkage        : 0.2", fixed = TRUE)
  expect_output(print(m), "bag.fraction     : 0.6", fixed = TRUE)
})

