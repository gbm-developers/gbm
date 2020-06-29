# these tests use examples from cv.glmnet documentation

library(glmnet)

teardown({
  detach("package:glmnet", unload = TRUE)
})

test_that("pmml.cv.glmnet throws no error when family is gaussian", {
  x <- matrix(rnorm(100 * 20), 100, 20)
  y <- rnorm(100)
  fit1 <- cv.glmnet(x, y)
  expect_silent(pmml(fit1))
})

test_that("pmml.cv.glmnet throws error when family is mgaussian", {
  x <- matrix(rnorm(100 * 20), 100, 20)
  y <- matrix(rnorm(100 * 3), 100, 3)
  fit1m <- cv.glmnet(x, y, family = "mgaussian")
  expect_error(pmml(fit1m),
    "Only poisson and gaussian family types supported.",
    fixed = TRUE
  )
})

test_that("pmml.cv.glmnet throws error when family is binomial", {
  x <- matrix(rnorm(100 * 20), 100, 20)
  g2 <- sample(1:2, 100, replace = TRUE)
  fit2 <- cv.glmnet(x, g2, family = "binomial")
  expect_error(pmml(fit2),
    "Only poisson and gaussian family types supported.",
    fixed = TRUE
  )
})

test_that("pmml.cv.glmnet throws error when family is multinomial", {
  x <- matrix(rnorm(100 * 20), 100, 20)
  g4 <- sample(1:4, 100, replace = TRUE)
  fit3 <- cv.glmnet(x, g4, family = "multinomial")
  fit3a <- cv.glmnet(x, g4, family = "multinomial", type.multinomial = "grouped")
  expect_error(pmml(fit3),
    "Only poisson and gaussian family types supported.",
    fixed = TRUE
  )
  expect_error(pmml(fit3a),
    "Only poisson and gaussian family types supported.",
    fixed = TRUE
  )
})

test_that("pmml.cv.glmnet throws no error when family is poisson", {
  N <- 500
  p <- 20
  nzc <- 5
  x <- matrix(rnorm(N * p), N, p)
  beta <- rnorm(nzc)
  f <- x[, seq(nzc)] %*% beta
  mu <- exp(f)
  y <- rpois(N, mu)
  fit <- cv.glmnet(x, y, family = "poisson")
  expect_silent(pmml(fit))
})

test_that("pmml.cv.glmnet throws error when family is cox", {
  set.seed(10101)
  N <- 1000
  p <- 30
  nzc <- p / 3
  x <- matrix(rnorm(N * p), N, p)
  beta <- rnorm(nzc)
  fx <- x[, seq(nzc)] %*% beta / 3
  hx <- exp(fx)
  ty <- rexp(N, hx)
  tcens <- rbinom(n = N, prob = .3, size = 1) # censoring indicator
  y <- cbind(time = ty, status = 1 - tcens) # y=Surv(ty,1-tcens) with library(survival)
  fit <- cv.glmnet(x, y, family = "cox")
  expect_error(pmml(fit),
    "Only poisson and gaussian family types supported.",
    fixed = TRUE
  )
})

test_that("pmml.cv.glmnet throws no error when input is sparse matrix", {
  # n=10000;p=200
  n <- 1000
  p <- 20
  nzc <- trunc(p / 10)
  x <- matrix(rnorm(n * p), n, p)
  iz <- sample(1:(n * p), size = n * p * .85, replace = FALSE)
  x[iz] <- 0
  sx <- Matrix(x, sparse = TRUE)
  beta <- rnorm(nzc)
  fx <- x[, seq(nzc)] %*% beta
  eps <- rnorm(n)
  y <- fx + eps
  px <- exp(fx)
  px <- px / (1 + px)
  ly <- rbinom(n = length(px), prob = px, size = 1)
  fit1 <- cv.glmnet(sx, y)
  fit2n <- cv.glmnet(x, y)
  expect_silent(pmml(fit1))
  expect_silent(pmml(fit2n))
})
