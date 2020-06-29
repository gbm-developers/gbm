context("utility-ks")

test_that("utility-ks", {
  mtdt <- list(y = c(0, 1, 0, 1, 1, 1),
               yhat = c(0.5, 0.9, 0.2, 0.7, 0.6,  0.4),
               w = rep(1, 6),
               pn = c(0, 0.5, 0.5, 1, 1, 1, 1),
               pp = c(0, 0, 0.25, 0.25, 0.5, 0.75, 1))
  expect_equal(round(ks(mtdt), 3), 0.75)
})
