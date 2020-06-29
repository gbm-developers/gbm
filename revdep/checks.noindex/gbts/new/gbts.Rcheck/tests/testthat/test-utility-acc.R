context("utility-acc")

mtdt <- list(y = c(0, 1, 0, 1, 1, 1),
             yhat = c(0.5, 0.9, 0.2, 0.7, 0.6,0.4),
             w = rep(1, 6),
             nn = c(1, 1, 2 ,2, 2, 2),
             np = c(0, 1, 1, 2, 3, 4),
             na = 1:6,
             pn = c(0, 0.5, 0.5, 1, 1, 1, 1),
             pp = c(0, 0, 0.25, 0.25, 0.5, 0.75, 1),
             pa = c(0, 1/6, 1/3, 1/2, 2/3, 5/6, 1),
             uprd = c(0.2, 0.4, 0.5, 0.6, 0.7, 0.9))

test_that("utility-acc", {

  # All scores within [0, 1]
  expect_equal(round(acc(mtdt), 3), 0.833)

  # Largest score <= 0.5
  mtdt$uprd <- c(0.1, 0.4, 0.3, 0.3, 0.4, 0.5)
  expect_equal(round(acc(mtdt), 3), 0.333)

  # Smallest score > 0.5
  mtdt$uprd <- c(0.6, 0.7, 0.8, 0.9, 0.7, 0.8)
  expect_equal(round(acc(mtdt), 3), 0.667)
})
