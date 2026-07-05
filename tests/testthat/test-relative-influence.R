test_that("relative.influence identifies the true predictors", {
  set.seed(1234)
  X1 <- matrix(nrow = 1000, ncol = 50)
  X1 <- apply(X1, 2, function(x) rnorm(1000))
  X2 <- matrix(nrow = 1000, ncol = 5)
  X2 <- apply(X2, 2, function(x) c(rnorm(500), rnorm(500, 3)))
  cls <- rep(c(0, 1), ea = 500)
  X <- data.frame(cbind(X1, X2, cls))
  mod <- gbm(
    cls ~ .,
    data = X,
    n.trees = 1000,
    cv.folds = 5,
    n.cores = 1,
    shrinkage = .01,
    interaction.depth = 2
  )
  ri <- rev(sort(relative.influence(mod)))
  wh <- names(ri)[1:5]
  res <- sum(wh %in% paste("V", 51:55, sep = ""))

  expect_identical(res, 5L)
})
