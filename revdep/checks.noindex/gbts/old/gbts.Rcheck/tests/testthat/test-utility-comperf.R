context("utility-comperf")

test_that("utility-comperf: binary response with weights = 1", {
  expect_equal(round(comperf(y = c(1, 0, 1, 0, 1, 1, 1, 1, 0, 1),
                             yhat = c(0.98, 0.09, 0.81, 0.35, 0.39, 0.49, 0.95,
                                      0.78, 0.22, 0.65),
                             pfmc = "dev"), 3), 0.678)
})

test_that("utility-comperf: binary response with weights > 0", {
  expect_equal(round(comperf(y = c(1, 0, 1, 0, 1, 1, 1, 1, 0, 1),
                             yhat = c(0.98, 0.09, 0.81, 0.35, 0.39, 0.49, 0.95,
                                      0.78, 0.22, 0.65),
                             w = c(1, 3, 2, 6, 4, 5, 3, 4, 2, 1),
                             pfmc = "dev"), 3), 0.821)
})

test_that("utility-comperf: binary response with tied scores", {
  expect_equal(round(comperf(y = c(1, 0, 1, 0, 1, 1, 1, 1, 0, 1),
                             yhat = c(0.98, 0.09, 0.81, 0.35, 0.39, 0.49, 0.78,
                                      0.78, 0.22, 0.65),
                             pfmc = "dev"), 3), 0.717)
})

test_that("utility-comperf: continuous response", {
  expect_equal(round(comperf(y = 1:5, yhat = 1:5 + c(0.1, 0.2, 0.2, 0.1, 0.2),
                       w = c(1, 2, 3, 2, 3), pfmc = "mse"), 3), 0.032)
})
