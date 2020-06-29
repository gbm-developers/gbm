skip_on_cran()
skip_on_appveyor()

test_that("The function trainRF produces the correct ouput", {
  data <- SDMtune:::t
  data@data <- data@data[, 1:4]
  m <- trainRF(data = data, mtry = 2, ntree = 200)
  expect_s4_class(m, "SDMmodel")
  expect_s4_class(m@model, "RF")
  expect_s4_class(m@data, "SWD")
  expect_equal(m@model@mtry, 2)
  expect_equal(m@model@ntree, 200)
  expect_equal(m@model@nodesize, 1)
  expect_equal(m@data, data)
})
