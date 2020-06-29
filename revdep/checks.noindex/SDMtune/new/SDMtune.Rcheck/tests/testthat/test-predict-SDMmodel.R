skip_on_cran()
skip_on_appveyor()

m <- SDMtune:::bm_maxent
m1 <- SDMtune:::bm_maxnet
train <- SDMtune:::t
train@data <- train@data[train@pa == 1, ]
train@coords <- train@coords[train@pa == 1, ]
train@pa <- train@pa[train@pa == 1]
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)

test_that("The method works with data frames", {
  p <- predict(m, train@data, "raw", clamp = FALSE)
  expect_length(p, nrow(train@data))
  expect_vector(p)
})

test_that("The method works with SWD objects", {
  p <- predict(m1, train, "logistic")
  expect_length(p, nrow(train@data))
  expect_vector(p)
})

test_that("The method works with raster stack objects", {
  p <- predict(m, predictors, "raw")
  expect_length(p, predictors$bio1@ncols * predictors$bio1@nrows)
  expect_s4_class(p, "RasterLayer")
})

test_that("The method works with raster stack objects and parallel", {
  p <- predict(m, predictors, "raw", parallel = TRUE)
  expect_length(p, predictors$bio1@ncols * predictors$bio1@nrows)
  expect_s4_class(p, "RasterLayer")
})
