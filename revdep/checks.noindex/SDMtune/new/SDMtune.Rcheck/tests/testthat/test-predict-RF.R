skip_on_cran()
skip_on_appveyor()

test_that("The ouput is correct", {
  data <- SDMtune:::t
  data@data <- data@data[, 1:4]
  m <- trainRF(data = data, mtry = 2, ntree = 200)
  pred <- predict(m@model, data@data)
  expect_equal(sum(pred >= 0), nrow(data@data))
  expect_equal(sum(pred <= 1), nrow(data@data))
})
