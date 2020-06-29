test_that("datadictionary error when features with unsupported class are used in tree model", {
  library(rpart)
  data(kyphosis)
  kyphosis[1, 2] <- ""
  kyphosis[2, 2] <- NA
  expect_error(pmml(rpart(formula = Kyphosis ~ Age + Number + Start, data = kyphosis)), "character class is not supported for features. Supported classes: numeric, logical, factor.")
})

test_that("datadictionary error when features with unsupported class are used in linear model (lm)", {
  test <- data.frame(x1 = rnorm(100), x2 = sample(c("a", "b"), 100, TRUE), y = rnorm(100), stringsAsFactors = FALSE)
  expect_error(pmml(lm(y ~ x1 + x2, data = test)), "character class is not supported for features. Supported classes: numeric, logical, factor.")
})
