skip_on_cran()
skip_on_appveyor()

test_that("Show method for RF class produces the correct output", {
  data <- SDMtune:::t
  data@data <- data@data[, 1:4]
  m <- train("RF", data = data, mtry = 2, ntree = 200)@model
  expect_output(print(m), "Class   : RF", fixed = TRUE)
  expect_output(print(m), "mtry    : 2", fixed = TRUE)
  expect_output(print(m), "ntree   : 200", fixed = TRUE)
  expect_output(print(m), "nodesize: 1", fixed = TRUE)
})
