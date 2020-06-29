skip_on_cran()
skip_on_appveyor()

files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)
set.seed(25)
suppressWarnings(bg <- dismo::randomPoints(predictors, 10000))
bg <- prepareSWD(species = "Bgs", a = bg, env = predictors,
                 categorical = "biome")
m <- SDMtune:::bm_maxnet

test_that("Exceptions are thrown", {
  expect_error(varSel(m, metric = "auc", bg4cor = bg, test = t, use_pc = TRUE),
               "Percent contribution cannot be used with model of")
})

test_that("Correlated Variable are removed", {
  set.seed(25, kind = "Mersenne-Twister", sample.kind = "Rejection")
  expect_message(o <- varSel(m, "auc", bg, test = t, cor_th = .9, permut = 1),
                 "Removed variables: bio16, bio6")
  expect_s4_class(o, "SDMmodel")
  expect_s4_class(o@model, "Maxnet")
  expect_false("bio16" %in% colnames(o@data@data))
  expect_false("bio6" %in% colnames(o@data@data))
})
