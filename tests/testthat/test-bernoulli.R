test_that("bernoulli predictions recover the signal", {
  set.seed(1)

  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4], N, replace = TRUE))
  mu <- c(-1, 0, 1, 2)[as.numeric(X3)]
  p <- 1 / (1 + exp(-(sin(3 * X1) - 4 * X2 + mu)))
  Y <- rbinom(N, 1, p)
  w <- rexp(N)
  w <- N * w / sum(w)
  data <- data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3)

  gbm1 <- gbm(
    Y ~ X1 + X2 + X3,
    data = data,
    weights = w,
    var.monotone = c(0, 0, 0),
    distribution = "bernoulli",
    n.trees = 3000,
    shrinkage = 0.001,
    interaction.depth = 3,
    bag.fraction = 0.5,
    train.fraction = 0.5,
    cv.folds = 5,
    n.cores = 1,
    n.minobsinnode = 10
  )

  best.iter <- gbm.perf(gbm1, method = "test", plot.it = FALSE)

  set.seed(2)
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4], N, replace = TRUE))
  mu <- c(-1, 0, 1, 2)[as.numeric(X3)]
  p <- 1 / (1 + exp(-(sin(3 * X1) - 4 * X2 + mu)))
  Y <- rbinom(N, 1, p)
  data2 <- data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3)

  f.predict <- predict(gbm1, data2, n.trees = best.iter)
  f.new <- sin(3 * X1) - 4 * X2 + mu

  expect_true(sd(f.new - f.predict) < 1.0)
})
