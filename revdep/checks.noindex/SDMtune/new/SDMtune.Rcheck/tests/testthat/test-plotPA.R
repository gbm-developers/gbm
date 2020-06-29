map <- raster::raster(matrix(runif(400, 0, 1), 20, 20))
file <- tempfile(fileext = ".asc")

test_that("The function raises an error if argument is not a raster", {
  expect_error(plotPA(data.frame(a = 1, b = "l")),
               "Prediction must be a RasterLayer object!")
})

test_that("The values are correct and the file is saved with correct format", {
  p <- plotPA(map, th = .8, filename = file, format = "ascii")
  expect_equal(p$plot_env$maxpixels, 50000)
  expect_setequal(p$data$value, c(FALSE, TRUE))
  expect_true(file.exists(file))
  # If hr is TRUE it should use the number of pixel in the raster
  p <- plotPA(map, th = 100, hr = TRUE)
  expect_equal(p$plot_env$maxpixels, 400)
})

teardown(unlink(file))
