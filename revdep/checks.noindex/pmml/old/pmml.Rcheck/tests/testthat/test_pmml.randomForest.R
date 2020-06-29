library(randomForest)

teardown({
  detach("package:randomForest", unload = TRUE)
})

test_that("error when object is not randomForest", {
  expect_error(pmml.randomForest("foo"), "Not a legitimate randomForest object")
})

test_that("no error occurs in doc example", {
  iris_rf <- randomForest(Species ~ ., data = iris, ntree = 20)
  expect_error(pmml(iris_rf), NA) # expect no error
})
