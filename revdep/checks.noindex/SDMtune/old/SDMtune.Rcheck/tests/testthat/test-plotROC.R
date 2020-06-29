p <- plotROC(SDMtune:::bm_maxent, test = SDMtune:::t)

test_that("The plot contains the correct data", {
  expect_setequal(unique(p$data$set), c("Train", "Test"))
  expect_equal(unique(p$data$pa), c(1, 0))
  expect_named(p$data, c("set", "pa", "pred"))

})

test_that("The plot has the correct labels", {
  expect_equal(p$labels$x, "False Positive Rate")
  expect_equal(p$labels$y, "True Positive Rate")
})
