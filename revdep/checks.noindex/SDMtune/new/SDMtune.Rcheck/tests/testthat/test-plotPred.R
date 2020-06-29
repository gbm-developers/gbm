map <- raster::raster(matrix(runif(400, 0, 1), 20, 20))

test_that("The function raises an error if argument is not a raster", {
  expect_error(plotPred(data.frame(a = 1, b = "l")),
               "Prediction must be a RasterLayer object!")
})

test_that("The values are correct", {
  p <- plotPred(map)
  expect_equal(p$plot_env$maxpixels, 50000)
  expect_true(min(p$data$value) >= 0)
  expect_true(max(p$data$value) <= 1)
  # If hr is TRUE it should use the number of pixel in the raster
  p <- plotPred(map, hr = TRUE)
  expect_equal(p$plot_env$maxpixels, 400)
})
