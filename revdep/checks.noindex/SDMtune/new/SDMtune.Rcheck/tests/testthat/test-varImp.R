skip_on_cran()
skip_on_appveyor()

test_that("The function produces the correct output with SDMmodel objects", {
  x <- varImp(SDMtune:::bm_maxent, permut = 2)
  expect_named(x, c("Variable", "Permutation_importance", "sd"))
  expect_equal(class(x), "data.frame")
  expect_equal(nrow(x), ncol(SDMtune:::bm_maxent@data@data))
  expect_setequal(x$Variable, colnames(SDMtune:::bm_maxent@data@data))
  # Column sd is not present for only one permutation
  expect_named(varImp(SDMtune:::bm_maxent, permut = 1),
               c("Variable", "Permutation_importance"))
})

test_that("The function produces the correct output with SDMmodelCV objects", {

  model <- SDMtune:::bm_maxent_cv
  pis <- vector("numeric", length = 4)
  df <- varImp(model, permut = 1)
  vars <- colnames(model@data@data)
  for (v in vars) {
    for (i in 1:4) {
      x <- varImp(model@models[[i]], permut = 1)
      pis[i] <- x[v, 2]
    }
    expect_equal(df[v, 2], mean(pis))
    expect_equal(df[v, 3], sd(pis))
  }
  expect_s3_class(df, "data.frame")
  expect_named(df,
               c("Variable", "Permutation_importance", "sd"))
  expect_equal(nrow(df), ncol(model@data@data))
  expect_setequal(x$Variable, colnames(model@data@data))
  expect_equal(sum(df[, 2]), 100, tolerance = 0.1)
})
