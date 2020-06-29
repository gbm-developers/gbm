# invisible(capture.output()) is used in catch xgboost's message about train-error

library(xgboost)
data(iris)

tmp_01_save <- tempfile()
tmp_01_dump <- tempfile()

teardown(unlink(c(tmp_01_save, tmp_01_dump), recursive = TRUE))


test_that("discrete variables are one-hot-encoded", {
  data(agaricus.train, package = "xgboost")
  train <- agaricus.train

  invisible(capture.output(model_fit <- xgboost(
    data = train$data, label = train$label,
    max_depth = 2, eta = 1, nthread = 2, nrounds = 2, objective = "binary:logistic",
    save_name = tmp_01_save
  )))

  xgb.dump(model_fit, tmp_01_dump)

  model_pmml <- pmml(
    model = model_fit, input_feature_names = colnames(train$data),
    output_label_name = "f", output_categories = c("0", "1"),
    xgb_dump_file = tmp_01_dump
  )

  expect_equal(length(model_pmml[[3]]), 4)
  expect_equal(xmlToList(model_pmml[[3]][[1]])[[1]], "odor")
  expect_equal(xmlToList(model_pmml[[3]][[1]])[[2]], "none")
  expect_equal(names(model_pmml)[[3]], "TransformationDictionary")
})


test_that("error is thrown when objective = reg:linear", {
  modX <- xgboost(
    data = as.matrix(iris[, 1:3]), label = iris[, 4],
    max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
    objective = "reg:linear", verbose = 0,
    save_name = tmp_01_save
  )

  xgb.dump(modX, tmp_01_dump)

  expect_error(
    pmml(
      model = modX, input_feature_names = colnames(iris[, 1:3]),
      output_label_name = "Petal.Width",
      xgb_dump_file = tmp_01_dump
    ),
    "Only the following objectives are supported: multi:softprob, multi:softmax, binary:logistic."
  )
})


test_that("error is thrown when objective = reg:logistic", {
  data(iris)

  dat_07 <- as.matrix(iris[1:100, 1:4])
  label_07 <- as.numeric(iris[1:100, 5]) - 1

  modX <- xgboost(
    data = dat_07, label = label_07, max_depth = 2,
    nrounds = 2, objective = "reg:logistic", verbose = 0,
    save_name = tmp_01_save
  )
  xgb.dump(modX, tmp_01_dump)

  expect_error(
    pmml(
      model = modX, input_feature_names = colnames(iris[, 1:4]),
      output_label_name = "Species",
      output_categories = c(0, 1),
      xgb_dump_file = tmp_01_dump
    ),
    "Only the following objectives are supported: multi:softprob, multi:softmax, binary:logistic."
  )
})

test_that("error is thrown when objective = binary:logitraw", {
  data(iris)
  ir <- iris[1:100, ]
  ir[, 5] <- as.character(ir[, 5])
  ir[, 5] <- as.factor(ir[, 5])

  model9 <- xgboost(
    data = as.matrix(ir[, 1:4]), label = as.numeric(ir[, 5]) - 1,
    max_depth = 3, nrounds = 3, objective = "binary:logitraw", verbose = 0,
    save_name = tmp_01_save
  )

  xgb.dump(model9, tmp_01_dump)

  expect_error(
    pmml(model9,
      input_feature_names = colnames(as.matrix(ir[, 1:4])),
      output_label_name = "Species",
      output_categories = c(1, 2),
      xgb_dump_file = tmp_01_dump
    ),
    "Only the following objectives are supported: multi:softprob, multi:softmax, binary:logistic."
  )
})
