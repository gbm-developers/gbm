skip_on_cran()
skip_on_appveyor()

files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
env <- raster::stack(files)
train <- SDMtune:::t

# Reduce observation to simulate k > n
np <- sum(train@pa == 1)
fold <- rep(TRUE, nrow(train@data))
fold[16:np] <- FALSE
train <- .subset_swd(train, fold)
m <- train("Maxnet", data = train, fc = "h")

test_that("NA is returned if k > obs", {
  expect_equal(aicc(m, env), NA)
})

test_that("The correct output is produced", {
  expect_type(aicc(SDMtune:::bm_maxent, env), "double")
})

test_that("Raises an error if called with the wrong model method", {
  data <- SDMtune:::t
  data@data <- data@data[, 1:4]
  m <- trainRF(data = data, mtry = 2, ntree = 200)
  expect_error(aicc(m, env),
               "AICc available only for \"Maxent\" and \"Maxnet\" methods.")
})
