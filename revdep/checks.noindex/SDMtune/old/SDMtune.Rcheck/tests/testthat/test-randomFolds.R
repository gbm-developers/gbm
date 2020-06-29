train <- SDMtune:::t

test_that("Error are raised", {
  expect_error(folds <- randomFolds(train@data, k = 2,
                       "\"data\" argument is not of class SWD."))
})

test_that("The output is correct for only_presence = FALSE", {
  folds <- randomFolds(train, k = 3, seed = 25)
  expect_length(folds, 2)
  expect_named(folds, c("train", "test"))
  expect_equal(ncol(folds$train), ncol(folds$test), 3)
  expect_equal(nrow(folds$train), nrow(folds$test))
  for (i in 1:3) {
    expect_equal(folds$train[, i], !folds$test[, i])
  }
})

test_that("The output is correct for only presence = TRUE", {
  folds <- randomFolds(train, k = 2, only_presence = TRUE)
  na <- sum(train@pa == 0)
  np <- nrow(folds$train) - na
  n <- np + na
  expect_length(folds, 2)
  expect_named(folds, c("train", "test"))
  expect_equal(ncol(folds$train), ncol(folds$test), 2)
  expect_equal(nrow(folds$train), nrow(folds$test))
  for (i in 1:2) {
    expect_equal(folds$train[, i][1:np], !folds$test[, i][1:np])
    expect_equal(folds$train[, i][(np + 1):n], folds$test[, i][(np + 1):n])
  }
})
