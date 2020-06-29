context("utility-dev")

test_that("utility-dev", {
  mtdt <- list(w = c(1, 1), wcn = 2, wcp = 2, uprd = 0.2)
  expect_equal(dev(mtdt), -2.0 * (sum(2 * log(0.2) + 2 * log(0.8))) / 2)
})
