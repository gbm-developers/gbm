library(neighbr)
data(iris)

teardown({
  detach("package:neighbr", unload = TRUE)
})

test_that("error when transform argument is not null", {
  set.seed(1131231)
  train_set <- iris[1:147, ]
  test_set <- iris[148:150, !names(iris) %in% c("Species")]

  expect_error(
    pmml(knn(
      train_set = train_set, test_set = test_set,
      k = 3,
      categorical_target = "Species",
      comparison_measure = "squared_euclidean"
    ), transforms = "NOT NULL"),
    "transforms currently not supported for knn models"
  )
})
