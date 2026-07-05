test_that("coxph predictions recover the signal", {
  set.seed(2)
  N <- 3000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4], N, replace = TRUE))
  mu <- c(-1, 0, 1, 2)[as.numeric(X3)]
  f <- 0.5 * sin(3 * X1 + 5 * X2 ^ 2 + mu / 10)
  tt.surv <- rexp(N, exp(f))
  tt.cens <- rexp(N, 0.5)
  delta <- as.numeric(tt.surv <= tt.cens)
  tt <- apply(cbind(tt.surv, tt.cens), 1, min)

  X1[sample(1:N, size = 100)] <- NA
  X3[sample(1:N, size = 300)] <- NA
  w <- rep(1, N)

  data <- data.frame(
    tt = tt,
    delta = delta,
    X1 = X1,
    X2 = X2,
    X3 = X3
  )

  gbm1 <- gbm(
    survival::Surv(tt, delta) ~ X1 + X2 + X3,
    data = data,
    weights = w,
    var.monotone = c(0, 0, 0),
    distribution = "coxph",
    n.trees = 3000,
    shrinkage = 0.001,
    interaction.depth = 3,
    bag.fraction = 0.5,
    train.fraction = 0.5,
    cv.folds = 5,
    n.cores = 1,
    n.minobsinnode = 10,
    keep.data = TRUE
  )

  best.iter <- gbm.perf(gbm1, method = "test", plot.it = FALSE)

  set.seed(2)
  N <- 1000
  X1 <- runif(N)
  X2 <- runif(N)
  X3 <- factor(sample(letters[1:4], N, replace = TRUE))
  mu <- c(-1, 0, 1, 2)[as.numeric(X3)]
  f <- 0.5 * sin(3 * X1 + 5 * X2 ^ 2 + mu / 10)
  tt.surv <- rexp(N, exp(f))
  tt.cens <- rexp(N, 0.5)

  data2 <- data.frame(
    tt = apply(cbind(tt.surv, tt.cens), 1, min),
    delta = as.numeric(tt.surv <= tt.cens),
    f = f,
    X1 = X1,
    X2 = X2,
    X3 = X3
  )

  f.predict <- predict(gbm1, newdata = data2, n.trees = best.iter)

  expect_true(sd(data2$f - f.predict) < 0.4)
})
