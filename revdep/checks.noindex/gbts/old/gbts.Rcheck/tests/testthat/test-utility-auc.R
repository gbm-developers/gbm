context("utility-auc")

test_that("utility-auc", {
  mtdt <- list(y = c(0, 1, 0, 1, 1, 1),
               yhat = c(0.5, 0.9, 0.2, 0.7, 0.6,  0.4),
               w = rep(1, 6),
               pn = c(0, 0.5, 0.5, 1, 1, 1, 1),
               pp = c(0, 0, 0.25, 0.25, 0.5, 0.75, 1),
               pa = c(0, 1/6, 1/3, 1/2, 2/3, 5/6, 1))
  expect_equal(round(auc(mtdt, "fpr", "tpr"), 3), 0.875)
})
