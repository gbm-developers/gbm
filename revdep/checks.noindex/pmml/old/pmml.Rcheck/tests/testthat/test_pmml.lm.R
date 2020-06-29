test_that("pmml.lm error when attempt is made to model with interaction terms", {
  test <- data.frame(x1 = rnorm(100), x2 = sample(c("a", "b"), 100, TRUE), y = rnorm(100), stringsAsFactors = FALSE)
  model <- lm(y ~ x1 * x2, data = test)
  expect_error(pmml(model), "Possible interaction terms detected. Please note that interaction terms for regression models are not presently supported.")
})

test_that("error when object is not lm", {
  a <- "foo"
  expect_error(pmml.lm(a), "Not a legitimate lm object")
})
