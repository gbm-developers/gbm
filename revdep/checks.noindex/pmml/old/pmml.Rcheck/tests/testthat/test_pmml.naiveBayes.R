test_that("error when object is not naiveBayes", {
  a <- "foo"
  expect_error(pmml.naiveBayes(a), "Not a legitimate naiveBayes object")
})
