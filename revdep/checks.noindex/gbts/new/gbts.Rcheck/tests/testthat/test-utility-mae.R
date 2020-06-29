context('utility-mae')

test_that("utility-mae: y = yhat", {
  mtdt <- list(y = 1:10,
               yhat = 1:10,
               w = rep(1, 10),
               rsa = 0,
               sow = 10)
  expect_equal(mae(mtdt), 0)
})

test_that("utility-mae: y != yhat", {
  mtdt <- list(y = 1:10,
               yhat = c(1:5 - 0.1, 6:10 + 0.1),
               w = rep(1, 10),
               rsa = 1,
               sow = 10)
  expect_equal(round(mae(mtdt), 3), 0.1)
})
