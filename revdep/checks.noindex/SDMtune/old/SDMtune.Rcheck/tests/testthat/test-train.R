skip_on_cran()
skip_on_appveyor()

train <- SDMtune:::t
folds <- randomFolds(train, k = 2)

test_that("Cross validation is executed", {
  cv <- train("Maxnet", data = train, folds = folds, fc = "l")
  expect_s4_class(cv, "SDMmodelCV")
  expect_length(cv@models, 2)
  expect_equal(cv@data, train)
  expect_length(cv@folds, 2)
  expect_equal(ncol(cv@folds[[1]]), 2)
  expect_equal(ncol(cv@folds[[2]]), 2)
})

test_that("Train without cross validation creates the correct output", {
  m <- train("Maxnet", data = train, fc = "l")
  expect_s4_class(m, "SDMmodel")
})

test_that("Train multiple methods creates the correct output", {
  # No errors if argument is not used
  expect_error(m <- train(c("Maxnet", "ANN"), data = train, fc = "l", size = 2,
                          ntree = 100),
               NA)
  expect_type(m, "list")
  expect_named(m, c("Maxnet", "ANN"))
  # Maxent model
  expect_s4_class(m$Maxnet, "SDMmodel")
  expect_s4_class(m$Maxnet@model, "Maxnet")
  expect_equal(m$Maxnet@model@fc, "l")
  expect_equal(m$Maxnet@model@reg, 1)
  # ANN model
  expect_s4_class(m$ANN, "SDMmodel")
  expect_s4_class(m$ANN@model, "ANN")
  expect_equal(m$ANN@model@size, 2)
  expect_equal(m$ANN@model@maxit, 100)
})
