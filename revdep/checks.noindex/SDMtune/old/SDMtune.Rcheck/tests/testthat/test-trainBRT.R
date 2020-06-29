skip_on_cran()
skip_on_appveyor()

test_that("The function trainBRT produces the correct ouput", {
  data <- SDMtune:::t
  data@data <- data@data[, 1:4]
  m <- trainBRT(data = data, n.trees = 200, shrinkage = 0.2)
  expect_s4_class(m, "SDMmodel")
  expect_s4_class(m@model, "BRT")
  expect_s4_class(m@data, "SWD")
  expect_equal(m@model@distribution, "bernoulli")
  expect_equal(m@model@n.trees, 200)
  expect_equal(m@model@interaction.depth, 1)
  expect_equal(m@model@shrinkage, 0.2)
  expect_equal(m@model@bag.fraction, 0.5)
  expect_equal(m@data, data)
})
