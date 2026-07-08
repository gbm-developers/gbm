test_that("gaussian predictions recover the signal", {
  set.seed(848)

  N <- 1000
  X1 <- runif(N)
  X2 <- 2 * runif(N)
  X3 <- factor(sample(letters[1:4], N, replace = TRUE))
  X4 <- ordered(sample(letters[1:6], N, replace = TRUE))
  X5 <- factor(sample(letters[1:3], N, replace = TRUE))
  X6 <- 3 * runif(N)
  mu <- c(-1, 0, 1, 2)[as.numeric(X3)]
  SNR <- 10
  Y <- X1 ** 1.5 + 2 * (X2 ** 0.5) + mu
  sigma <- sqrt(var(Y) / SNR)
  Y <- Y + rnorm(N, mean = 0, sd = sigma)
  X1[sample(1:N, size = 100)] <- NA
  X3[sample(1:N, size = 300)] <- NA
  data <- data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, X6 = X6)

  gbm1 <- gbm(
    Y ~ X1 + X2 + X3 + X4 + X5 + X6,
    data = data,
    var.monotone = c(0, 0, 0, 0, 0, 0),
    distribution = "gaussian",
    n.trees = 2000,
    shrinkage = 0.005,
    interaction.depth = 3,
    bag.fraction = 0.5,
    train.fraction = 1,
    n.minobsinnode = 10,
    keep.data = TRUE,
    cv.folds = 10,
    n.cores = 1
  )

  best.iter <- gbm.perf(gbm1, method = "cv", plot.it = FALSE)

  # predict() re-runs the serialized trees, whereas gbm1$fit comes directly
  # from the C++ training loop; they must agree. This guards the tree
  # serialization/prediction path (e.g. the categorical split codes that broke
  # on Linux ARM64, where plain char is unsigned and -1 became 255).
  p.train <- predict(gbm1, data, gbm1$n.trees)
  expect_lt(max(abs(gbm1$fit - p.train)), 1e-8)

  set.seed(223)
  N <- 1000
  X1 <- runif(N)
  X2 <- 2 * runif(N)
  X3 <- factor(sample(letters[1:4], N, replace = TRUE))
  X4 <- ordered(sample(letters[1:6], N, replace = TRUE))
  X5 <- factor(sample(letters[1:3], N, replace = TRUE))
  X6 <- 3 * runif(N)
  mu <- c(-1, 0, 1, 2)[as.numeric(X3)]
  Y <- X1 ** 1.5 + 2 * (X2 ** 0.5) + mu
  data2 <- data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, X6 = X6)

  f.predict <- predict(gbm1, data2, best.iter)

  expect_true(abs(mean(data2$Y - f.predict)) < 0.01)
  expect_true(sd(data2$Y - f.predict) < sigma)
})
