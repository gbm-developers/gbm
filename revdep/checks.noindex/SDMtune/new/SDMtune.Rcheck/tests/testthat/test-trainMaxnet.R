skip_on_cran()
skip_on_appveyor()

test_that("The function trainMaxent produces the correct ouput", {
  m <- trainMaxnet(data = SDMtune:::t, reg = 1.2, fc = "l")
  expect_s4_class(m, "SDMmodel")
  expect_s4_class(m@model, "Maxnet")
  expect_s4_class(m@data, "SWD")
  expect_equal(m@model@reg, 1.2)
  expect_equal(m@model@fc, "l")
  expect_equal(m@data, SDMtune:::t)
})
