test_that("The function clamps correctly", {
  m <- matrix(c(0, 5, 10, 10, 15, 20), ncol = 2)
  data <- scaleClamp(m, c(2, 12), c(9, 18), TRUE, FALSE)
  expect_equal(data, matrix(c(2, 5, 9, 12, 15, 18), ncol = 2))
  expect_equal(length(data), 6)
  expect_equal(ncol(data), 2)
})

test_that("The function scales correctly", {
  m <- matrix(c(0, 5, 10, 10, 15, 20), ncol = 2)
  data <- scaleClamp(m, c(0, 10), c(10, 20), FALSE, TRUE)
  expect_equal(data, matrix(c(0.0, 0.5, 1.0, 0.0, 0.5, 1.0), ncol = 2))
  expect_equal(length(data), 6)
  expect_equal(ncol(data), 2)
})

test_that("The function scales and clamps correctly", {
  m <- matrix(c(-1, 5, 12, 9, 15, 22), ncol = 2)
  data <- scaleClamp(m, c(0, 10), c(10, 20), TRUE, TRUE)
  expect_equal(data, matrix(c(0.0, 0.5, 1.0, 0.0, 0.5, 1.0), ncol = 2))
  expect_equal(length(data), 6)
  expect_equal(ncol(data), 2)
})
