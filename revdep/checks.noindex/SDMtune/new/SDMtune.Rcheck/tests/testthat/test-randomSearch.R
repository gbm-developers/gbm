skip_on_cran()
skip_on_appveyor()

files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)
h <- list(fc = c("l", "q", "p"), reg = 1:3)
o <- randomSearch(SDMtune:::bm_maxnet, h, "aicc", pop = 3, env = predictors,
                  seed = 25)

test_that("randomSearch produces the expected output", {
  expect_s4_class(o, "SDMtune")
  expect_s4_class(o@models[[1]], "SDMmodel")
  expect_s3_class(o@results, "data.frame")
  expect_named(o@results, c("fc", "reg", "AICc", "delta_AICc"))
  expect_length(o@models, 3)
  expect_equal(o@results$fc, c("q", "l", "l"))
  expect_equal(o@results$reg, c(2, 2, 3))
})
