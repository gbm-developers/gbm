context('utility-rsq')

test_that("utility-rsq: y = yhat", {
  mtdt <- list(y = 1:10,
               yhat = 1:10,
               w = rep(1, 10),
               rss = 0,
               sow = 10)
  expect_equal(rsq(mtdt), 1)
})

test_that("utility-rsq: y != yhat", {
  mtdt <- list(y = 1:10,
               yhat = c(1:5 - 0.1, 6:10 + 0.1),
               w = rep(1, 10),
               rss = 0.1,
               sow = 10)
  expect_equal(round(rsq(mtdt), 3), 0.999)
})

