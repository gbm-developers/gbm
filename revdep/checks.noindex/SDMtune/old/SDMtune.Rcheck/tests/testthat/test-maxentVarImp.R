test_that("Error are raised", {
  expect_error(maxentVarImp(SDMtune:::bm_maxnet))
})

test_that("The output is correct for SDMmodel objects", {
  model <- SDMtune:::bm_maxent
  df <- maxentVarImp(model)
  expect_s3_class(df, "data.frame")
  expect_named(df,
               c("Variable", "Percent_contribution", "Permutation_importance"))
  expect_equal(nrow(df), ncol(model@data@data))
})

test_that("The output is correct for SDMmodelCV objects", {
  model <- SDMtune:::bm_maxent_cv
  pcs <- pis <- vector("numeric", length = 4)
  df <- maxentVarImp(model)
  vars <- colnames(model@data@data)
  for (v in vars) {
    for (i in 1:4) {
      x <- maxentVarImp(model@models[[i]])
      pcs[i] <- x[v, 2]
      pis[i] <- x[v, 3]
    }
    expect_equal(df[v, 2], mean(pcs))
    expect_equal(df[v, 3], mean(pis))
  }
  expect_s3_class(df, "data.frame")
  expect_named(df,
               c("Variable", "Percent_contribution", "Permutation_importance"))
  expect_equal(nrow(df), ncol(model@data@data))
})
