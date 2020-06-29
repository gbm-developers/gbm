skip_on_cran()
skip_on_appveyor()

test_that("The function trainANN produces the correct ouput", {
  data <- SDMtune:::t
  data@data <- data@data[, 1:4]
  m <- trainANN(data = data, size = 10)
  expect_s4_class(m, "SDMmodel")
  expect_s4_class(m@model, "ANN")
  expect_s4_class(m@data, "SWD")
  expect_equal(m@model@size, 10)
  expect_equal(m@model@decay, 0)
  expect_equal(m@model@rang, 0.7)
  expect_equal(m@model@maxit, 100)
  expect_equal(m@data, data)
})
