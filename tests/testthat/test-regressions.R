test_that("gbm.more works for coxph models fit with keep.data = FALSE", {
  set.seed(10)

  n <- 80
  cox_data <- data.frame(
    time = rexp(n),
    status = rbinom(n, 1, 0.6),
    x = rnorm(n)
  )
  cox_fit <- gbm(
    survival::Surv(time, status) ~ x,
    data = cox_data,
    distribution = "coxph",
    n.trees = 5,
    interaction.depth = 1,
    n.minobsinnode = 5,
    bag.fraction = 0.8,
    keep.data = FALSE,
    verbose = FALSE
  )

  expect_s3_class(
    gbm.more(cox_fit, n.new.trees = 1, data = cox_data, verbose = FALSE),
    "gbm"
  )
})

test_that("gbm.more works for pairwise models fit with keep.data = FALSE", {
  set.seed(10)

  pairwise_data <- data.frame(
    y = rep(c(0, 1), 20),
    x = rnorm(40),
    query = rep(1:20, each = 2)
  )
  pairwise_fit <- gbm(
    y ~ x,
    data = pairwise_data,
    distribution = list(name = "pairwise", group = "query", metric = "ndcg"),
    n.trees = 5,
    interaction.depth = 1,
    n.minobsinnode = 1,
    bag.fraction = 0.8,
    keep.data = FALSE,
    verbose = FALSE
  )

  expect_s3_class(
    gbm.more(pairwise_fit, n.new.trees = 1, data = pairwise_data, verbose = FALSE),
    "gbm"
  )
})

test_that("gbm.more works for multinomial models with kept data", {
  set.seed(10)

  expect_warning(multinomial_fit <- gbm(
    Species ~ .,
    data = iris,
    distribution = "multinomial",
    n.trees = 5,
    interaction.depth = 1,
    n.minobsinnode = 5,
    bag.fraction = 0.8,
    verbose = FALSE
  ), NA)

  extended_fit <- gbm.more(multinomial_fit, n.new.trees = 2, verbose = FALSE)
  predictions <- predict(extended_fit, iris, n.trees = 7, type = "response")

  expect_s3_class(extended_fit, "gbm")
  expect_equal(extended_fit$n.trees, 7)
  expect_equal(dim(extended_fit$fit), c(nrow(iris), nlevels(iris$Species)))
  expect_equal(dim(predictions), c(nrow(iris), nlevels(iris$Species), 1))
  expect_equal(rowSums(predictions[, , 1]), rep(1, nrow(iris)))
})

test_that("gbm.more works for multinomial models fit with keep.data = FALSE", {
  set.seed(10)

  expect_warning(multinomial_fit <- gbm(
    Species ~ .,
    data = iris,
    distribution = "multinomial",
    n.trees = 5,
    interaction.depth = 1,
    n.minobsinnode = 5,
    bag.fraction = 0.8,
    keep.data = FALSE,
    verbose = FALSE
  ), NA)

  extended_fit <- gbm.more(multinomial_fit, n.new.trees = 2, data = iris,
                           verbose = FALSE)
  predictions <- predict(extended_fit, iris, n.trees = 7, type = "response")

  expect_s3_class(extended_fit, "gbm")
  expect_equal(extended_fit$n.trees, 7)
  expect_equal(dim(extended_fit$fit), c(nrow(iris), nlevels(iris$Species)))
  expect_equal(dim(predictions), c(nrow(iris), nlevels(iris$Species), 1))
  expect_equal(rowSums(predictions[, , 1]), rep(1, nrow(iris)))
})

test_that("gbmDoFold reorders weights to match fold rows", {
  fold_x <- data.frame(x = seq_len(40))
  fold_y <- seq_len(40)
  fold_w <- seq_len(40)
  fold_group <- rep(1:2, each = 20)
  cv_group <- rep(1:2, length.out = 40)
  i <- order(cv_group == 1)
  fold_model <- gbmDoFold(
    X = 1,
    i.train = 1:40,
    x = fold_x,
    y = fold_y,
    offset = NULL,
    distribution = list(name = "gaussian"),
    w = fold_w,
    var.monotone = 0,
    n.trees = 3,
    interaction.depth = 1,
    n.minobsinnode = 2,
    shrinkage = 0.1,
    bag.fraction = 0.8,
    cv.group = cv_group,
    var.names = "x",
    response.name = "y",
    group = fold_group,
    s = 1L
  )
  set.seed(1L)
  expected_model <- gbm.fit(
    x = fold_x[i, , drop = FALSE],
    y = fold_y[i],
    offset = NULL,
    distribution = list(name = "gaussian"),
    w = fold_w[i],
    var.monotone = 0,
    n.trees = 3,
    interaction.depth = 1,
    n.minobsinnode = 2,
    shrinkage = 0.1,
    bag.fraction = 0.8,
    nTrain = sum(cv_group != 1),
    keep.data = FALSE,
    verbose = FALSE,
    response.name = "y",
    group = fold_group[i]
  )

  expect_equal(fold_model$train.error, expected_model$train.error)
  expect_equal(fold_model$valid.error, expected_model$valid.error)
})

test_that("plot.gbm default palette works without requiring viridis", {
  set.seed(10)

  plot_data <- data.frame(
    y = rnorm(100),
    x1 = rnorm(100),
    x2 = rnorm(100)
  )
  plot_fit <- gbm(
    y ~ x1 + x2,
    data = plot_data,
    distribution = "gaussian",
    n.trees = 5,
    interaction.depth = 1,
    n.minobsinnode = 5,
    bag.fraction = 0.8,
    verbose = FALSE
  )

  expect_s3_class(plot(plot_fit, i.var = 1:2, n.trees = 1), "trellis")
})

test_that("quantile regression supports non-constant weights", {
  quantile_data <- data.frame(
    y = 1:20,
    x = 1:20,
    w = c(rep(1, 19), 100)
  )

  expect_warning(
    quantile_fit <- gbm(
      y ~ x,
      data = quantile_data,
      weights = w,
      distribution = list(name = "quantile", alpha = 0.5),
      n.trees = 1,
      interaction.depth = 1,
      n.minobsinnode = 2,
      bag.fraction = 1,
      train.fraction = 1,
      verbose = FALSE
    ),
    NA
  )

  expect_equal(quantile_fit$initF, 20)

  terminal_data <- data.frame(
    y = c(rep(0, 10), 100, 150, 200),
    x = c(rep(0, 10), 1, 1, 1),
    w = c(rep(100, 10), 1, 1, 100)
  )
  terminal_fit <- gbm(
    y ~ x,
    data = terminal_data,
    weights = w,
    distribution = list(name = "quantile", alpha = 0.5),
    n.trees = 1,
    interaction.depth = 1,
    n.minobsinnode = 1,
    bag.fraction = 1,
    train.fraction = 1,
    shrinkage = 1,
    verbose = FALSE
  )

  expect_equal(unique(predict(terminal_fit, terminal_data[terminal_data$x == 1, ],
                              n.trees = 1)), 200)
})
