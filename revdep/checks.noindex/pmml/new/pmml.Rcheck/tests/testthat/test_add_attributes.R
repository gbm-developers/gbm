test_that("test add_attributes function follows the xml path properly", {
  model0 <- lm(Sepal.Length ~ ., data = iris[, -5])
  model <- pmml(model0)

  # add arbitary attributes on the 3rd 'NumericPredictor' element
  # with an attribute 'exponent'=1
  model <- add_attributes(model,
    "/p:PMML/descendant::p:NumericPredictor[@exponent='1'][3]",
    attributes = c(a = 1, b = "b")
  )

  test <- xmlToList(model)$RegressionModel$RegressionTable[[3]]
  expect_equal(test[[1]], "Petal.Width")
  expect_equal(test[[4]], "1")
  expect_equal(names(test)[4], "a")
  expect_equal(test[[5]], "b")
})
