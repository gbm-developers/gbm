skip_on_cran()
skip_on_appveyor()

folder <- "trash"
m <- SDMtune:::bm_maxnet
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
env <- raster::stack(files)

modelReport(m, type = "cloglog", folder = folder, test = SDMtune:::t,
            permut = 1, env = env)
test_that("The files are created", {
  expect_true(file.exists(file.path(folder, "train.csv")))
  expect_true(file.exists(file.path(folder, "test.csv")))
  expect_true(file.exists(file.path(folder, "model.Rds")))
  expect_true(file.exists(file.path(folder, "map.tif")))
  expect_true(file.exists(file.path(folder, "virtual_species.html")))
  expect_true(file.exists(file.path(folder, "plots", "ROC_curve.png")))
})

teardown(unlink(file.path(getwd(), folder), recursive = TRUE))
