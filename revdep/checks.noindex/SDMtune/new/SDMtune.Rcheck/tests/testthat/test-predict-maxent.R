skip_on_cran()
skip_on_travis()
skip_on_appveyor()
skip_on_covr()

model <- SDMtune:::bm_maxent
maxent_model <- SDMmodel2MaxEnt(model)
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
env <- raster::stack(files)
data <- data_cont <- data_cat <- SDMtune:::t
# Remove categoriacal variable
data_cont@data <- data@data[, 1:7]
# Remove continuous variables
data_cat@data <- data@data[, 8, drop = FALSE]

test_that("The function predicts cloglog correctly", {
  expect_equal(predict(model, data, type = "cloglog"),
               predict(maxent_model, data@data),
               tolerance = 1e-7)
})

test_that("The function predicts logistic correctly", {
  expect_equal(predict(model, data, type = "logistic"),
               predict(maxent_model, data@data, args = "outputformat=logistic"),
               tolerance = 1e-7)
})

test_that("The function predicts raw correctly", {
  expect_equal(predict(model, data, type = "raw"),
               predict(maxent_model, data@data, args = "outputformat=raw"),
               tolerance = 1e-7)
})

test_that("The function predicts raster correctly", {
  expect_equal(predict(model, env, type = "cloglog"),
               predict(maxent_model, env, args = "outputformat=cloglog"),
               tolerance = 1e-7)
})

model <- train("Maxent", data = data_cont, fc = "l")
maxent_model <- SDMmodel2MaxEnt(model)

test_that("The function predicts fc l correctly", {
  expect_equal(predict(model, data_cont, type = "raw"),
               predict(maxent_model, data_cont@data, args = "outputformat=raw"),
               tolerance = 1e-7)
})

model <- train("Maxent", data = data_cont, fc = "q")
maxent_model <- SDMmodel2MaxEnt(model)

test_that("The function predicts fc q correctly", {
  expect_equal(predict(model, data_cont, type = "raw"),
               predict(maxent_model, data_cont@data, args = "outputformat=raw"),
               tolerance = 1e-7)
})

model <- train("Maxent", data = data_cont, fc = "p")
maxent_model <- SDMmodel2MaxEnt(model)

test_that("The function predicts fc p correctly", {
  expect_equal(predict(model, data_cont, type = "raw"),
               predict(maxent_model, data_cont@data, args = "outputformat=raw"),
               tolerance = 1e-7)
})

model <- train("Maxent", data = data_cont, fc = "h")
maxent_model <- SDMmodel2MaxEnt(model)

test_that("The function predicts fc h correctly", {
  expect_equal(predict(model, data_cont, type = "raw"),
               predict(maxent_model, data_cont@data, args = "outputformat=raw"),
               tolerance = 1e-7)
})

model <- train("Maxent", data = data_cont, fc = "t")
maxent_model <- SDMmodel2MaxEnt(model)

test_that("The function predicts fc t correctly", {
  expect_equal(predict(model, data_cont, type = "raw"),
               predict(maxent_model, data_cont@data, args = "outputformat=raw"),
               tolerance = 1e-7)
})

model <- train("Maxent", data = data_cat, fc = "t")
maxent_model <- SDMmodel2MaxEnt(model)

test_that("The function predicts fc categorical correctly", {
  expect_equal(predict(model, data_cat, type = "cloglog"),
               predict(maxent_model, data_cat@data),
               tolerance = 1e-7)
})

data@data <- data@data[, "bio1", drop = FALSE]
m <- train("Maxent", data = data, fc = "l")

test_that("The function works when using a single variable and a single FC", {
  expect_length(predict(m, data, type = "raw"), nrow(data@data))
})
