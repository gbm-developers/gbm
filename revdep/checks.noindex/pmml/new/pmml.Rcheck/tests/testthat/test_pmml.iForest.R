test_that("error when object is not iForest", {
  expect_error(pmml.iForest("foo"), "Not a legitimate iForest object")
})

test_that("confirm non-existent category is still automatically created by iForest function", {
  skip_on_cran()
  skip_on_ci()

  library(isofor)
  data(iris)

  mod <- iForest(iris, nt = 2, phi = 30)
  model_pmml <- pmml(model = mod)
  expect_equal(length(model_pmml[[2]]), 5)
  expect_equal(as.character(model_pmml[[2]][[5]][[4]])[2], ".")
})
