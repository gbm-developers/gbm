skip_on_cran()
skip_on_appveyor()

m <- SDMtune:::bm_maxent_cv
train <- SDMtune:::t
train@data <- train@data[train@pa == 1, ]
train@coords <- train@coords[train@pa == 1, ]
train@pa <- train@pa[train@pa == 1]
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)
folder <- tempfile("SDMtune")
dir.create(folder)

test_that("The method works with data frames", {
  p <- predict(m, train@data, type = "raw")
  expect_length(p, nrow(train@data))
  expect_vector(p)
})

test_that("The method works with SWD objects", {
  p <- predict(m, train, type = "raw")
  expect_length(p, nrow(train@data))
  expect_vector(p)
})

test_that("The method works with raster stack objects", {
  p <- predict(m, predictors, type = "raw")
  expect_length(p, predictors$bio1@ncols * predictors$bio1@nrows)
  expect_s4_class(p, "RasterLayer")
})

test_that("The method works with raster stack objects and parallel", {
  p <- predict(m, predictors, type = "raw", parallel = TRUE)
  expect_length(p, predictors$bio1@ncols * predictors$bio1@nrows)
  expect_s4_class(p, "RasterLayer")
  expect_false(getOption("SDMtuneParallel"))
})

test_that("The output is the function applied to the k predictions", {
  train@data <- train@data[1:3, ]

  p <- predict(m, train@data, fun = c("mean", "sd", "min"), type = "raw")

  expect_equal(class(p), "list")
  expect_vector(p, size = 3)
  expect_named(p, c("mean", "sd", "min"))

  # mean
  preds <- matrix(nrow = 3, ncol = 4)
  for (i in 1:4) {
    preds[, i] <- predict(m@models[[i]], train@data, type = "raw")
  }
  for (i in 1:3) {
    expect_equal(p$mean[i], mean(preds[i, ]))
  }

  # sd
  preds <- matrix(nrow = 3, ncol = 4)
  for (i in 1:4) {
    preds[, i] <- predict(m@models[[i]], train@data, type = "raw")
  }
  for (i in 1:3) {
    expect_equal(p$sd[i], sd(preds[i, ]))
  }

  # min
  preds <- matrix(nrow = 3, ncol = 4)
  for (i in 1:4) {
    preds[, i] <- predict(m@models[[i]], train@data, type = "raw")
  }
  for (i in 1:3) {
    expect_equal(p$min[i], min(preds[i, ]))
  }
})

test_that("The function works with raster data and multiple functions", {

  funs <- c("mean", "sd", "min")
  p <- predict(m, predictors, fun = funs, type = "raw",
               filename = file.path(folder, "pred"))

  expect_equal(class(p), "list")
  expect_vector(p, size = 3)
  expect_named(p, c("mean", "sd", "min"))

  # check that files are created
  for (i in 1:3) {
    expect_true(file.exists(paste0(file.path(folder, "pred_"), funs[i],
                                   ".tif")))
    expect_s4_class(p[[i]], "RasterLayer")
  }

})

teardown(unlink(folder))
