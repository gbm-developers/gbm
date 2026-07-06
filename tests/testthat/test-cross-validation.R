test_that("cross-validation fitted values are stored for the selected iteration", {
  set.seed(101)

  cv_data <- data.frame(
    y = rnorm(60),
    x1 = runif(60),
    x2 = factor(sample(letters[1:3], 60, replace = TRUE))
  )

  cv_fit <- gbm(
    y ~ x1 + x2,
    data = cv_data,
    distribution = "gaussian",
    n.trees = 10,
    interaction.depth = 1,
    n.minobsinnode = 3,
    shrinkage = 0.05,
    bag.fraction = 0.8,
    train.fraction = 1,
    cv.folds = 3,
    n.cores = 1,
    verbose = FALSE
  )

  expect_length(cv_fit$cv.error, cv_fit$n.trees)
  expect_length(cv_fit$cv.fitted, nrow(cv_data))
  expect_false(anyNA(cv_fit$cv.fitted))
})

test_that("cross-validation fitted values preserve multinomial columns", {
  set.seed(102)

  cv_data <- iris[iris$Species != "virginica", ]
  cv_data$Species <- droplevels(cv_data$Species)

  cv_fit <- gbm(
    Species ~ Sepal.Length + Sepal.Width,
    data = cv_data,
    distribution = "multinomial",
    n.trees = 5,
    interaction.depth = 1,
    n.minobsinnode = 3,
    shrinkage = 0.05,
    bag.fraction = 0.8,
    train.fraction = 1,
    cv.folds = 2,
    n.cores = 1,
    verbose = FALSE
  )

  expect_equal(dim(cv_fit$cv.fitted), c(nrow(cv_data), cv_fit$num.classes))
  expect_false(anyNA(cv_fit$cv.fitted))
})

test_that("cross-validation fitted values are stored for pairwise models", {
  set.seed(103)

  n <- 60
  cv_data <- data.frame(
    y = runif(n),
    query = sample(rep(seq_len(12), each = 5)),
    x1 = runif(n),
    x2 = rnorm(n)
  )

  cv_fit <- gbm(
    y ~ x1 + x2,
    data = cv_data,
    distribution = list(name = "pairwise", group = "query", metric = "ndcg"),
    n.trees = 5,
    interaction.depth = 1,
    n.minobsinnode = 2,
    shrinkage = 0.05,
    bag.fraction = 0.8,
    train.fraction = 1,
    cv.folds = 3,
    n.cores = 1,
    verbose = FALSE
  )

  expect_length(cv_fit$cv.fitted, nrow(cv_data))
  expect_false(anyNA(cv_fit$cv.fitted))
})

test_that("one cross-validation fold is treated as no cross-validation", {
  set.seed(104)

  cv_data <- data.frame(
    y = rnorm(40),
    x = runif(40)
  )

  expect_warning(
    cv_fit <- gbm(
      y ~ x,
      data = cv_data,
      distribution = "gaussian",
      n.trees = 5,
      interaction.depth = 1,
      n.minobsinnode = 3,
      shrinkage = 0.05,
      bag.fraction = 0.8,
      train.fraction = 1,
      cv.folds = 1,
      verbose = FALSE
    ),
    "cv.folds = 1 is not meaningful"
  )

  expect_equal(cv_fit$cv.folds, 0)
  expect_null(cv_fit$cv.error)
  expect_null(cv_fit$cv.fitted)
})
