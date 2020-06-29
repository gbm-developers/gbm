library(xgboost)
library(randomForest)
data(agaricus.train, package = "xgboost")

tmp_02_save <- tempfile()
tmp_02_dump <- tempfile()
tmp_03_save <- tempfile()
tmp_03_dump <- tempfile()
teardown(unlink(c(tmp_02_save, tmp_02_dump, tmp_03_save, tmp_03_dump), recursive = TRUE))


test_that("default invalidValueTreatment attribute is exported correctly for linear models", {
  lm_model_0 <- lm(Sepal.Length ~ ., data = iris[, -5])
  lm_model_1 <- pmml(lm_model_0)

  ms <- xmlToList(lm_model_1)$RegressionModel$MiningSchema

  expect_equal(unlist(ms), c(
    "Sepal.Length", "predicted", "returnInvalid", "Sepal.Width",
    "active", "returnInvalid", "Petal.Length", "active",
    "returnInvalid", "Petal.Width", "active", "returnInvalid"
  ))
})

test_that("invalidValueTreatment attribute is exported correctly for xgboost models", {
  train <- agaricus.train
  invisible(capture.output(model_fit <- xgboost(
    data = train$data, label = train$label,
    max_depth = 2, eta = 1, nthread = 2, nrounds = 2, objective = "binary:logistic",
    save_name = tmp_02_save
  )))
  xgb.dump(model_fit, tmp_02_dump)

  # default invalidValueTreatment arguments
  model_pmml <- pmml(
    model = model_fit, input_feature_names = colnames(train$data),
    output_label_name = "f", output_categories = c("0", "1"),
    xgb_dump_file = tmp_02_dump
  )
  # parent segment
  ms2 <- unlist(xmlToList(model_pmml)$MiningModel$MiningSchema)
  expect_equal(ms2, c(
    "odor", "active", "returnInvalid", "stalk-root", "active", "returnInvalid",
    "spore-print-color", "active", "returnInvalid", "f", "predicted", "returnInvalid"
  ))
  # child segment 0
  ms3 <- unlist(xmlToList(model_pmml)$MiningModel$Segmentation[[2]]$MiningSchema)
  expect_equal(ms3, c(
    "odor", "active", "asIs", "stalk-root", "active", "asIs",
    "spore-print-color", "active", "asIs", "f", "predicted", "asIs"
  ))
  # child segment 1
  ms4 <- unlist(xmlToList(model_pmml)$MiningModel$Segmentation[[5]]$MiningSchema)
  expect_equal(ms4, c(
    "odor", "active", "asIs", "stalk-root", "active", "asIs",
    "spore-print-color", "active", "asIs", "f", "predicted", "asIs"
  ))

  # child segment 2 - the regression model segment
  ms5 <- unlist(xmlToList(model_pmml)$MiningModel$Segmentation[[8]]$MiningSchema)
  expect_equal(ms5, c(
    "predictedValueTree0", "active", "continuous", "asIs",
    "predictedValueTree1", "active", "continuous", "asIs"
  ))

  # non-default invalidValueTreatment arguments - 1
  model_pmml_2 <- pmml(
    model = model_fit, input_feature_names = colnames(train$data),
    output_label_name = "f", output_categories = c("0", "1"),
    xgb_dump_file = tmp_02_dump,
    parent_invalid_value_treatment = "returnInvalid",
    child_invalid_value_treatment = "returnInvalid"
  )
  # parent segment
  ms22 <- xmlToList(model_pmml_2)$MiningModel$MiningSchema
  expect_equal(unlist(ms22), c(
    "odor", "active", "returnInvalid", "stalk-root", "active", "returnInvalid",
    "spore-print-color", "active", "returnInvalid", "f", "predicted", "returnInvalid"
  ))
  # child segment 0
  ms23 <- unlist(xmlToList(model_pmml_2)$MiningModel$Segmentation[[2]]$MiningSchema)
  expect_equal(ms23, c(
    "odor", "active", "returnInvalid", "stalk-root", "active", "returnInvalid",
    "spore-print-color", "active", "returnInvalid", "f", "predicted", "returnInvalid"
  ))
  # child segment 1
  ms24 <- unlist(xmlToList(model_pmml_2)$MiningModel$Segmentation[[5]]$MiningSchema)
  expect_equal(ms24, c(
    "odor", "active", "returnInvalid", "stalk-root", "active", "returnInvalid",
    "spore-print-color", "active", "returnInvalid", "f", "predicted", "returnInvalid"
  ))

  # child segment 2 - the regression model segment
  ms25 <- unlist(xmlToList(model_pmml_2)$MiningModel$Segmentation[[8]]$MiningSchema)
  expect_equal(ms25, c(
    "predictedValueTree0", "active", "continuous", "returnInvalid",
    "predictedValueTree1", "active", "continuous", "returnInvalid"
  ))

  # non-default invalidValueTreatment arguments - 2
  model_pmml_3 <- pmml(
    model = model_fit, input_feature_names = colnames(train$data),
    output_label_name = "f", output_categories = c("0", "1"),
    xgb_dump_file = tmp_02_dump,
    parent_invalid_value_treatment = "asIs"
  )
  # parent segment
  ms32 <- xmlToList(model_pmml_3)$MiningModel$MiningSchema
  expect_equal(unlist(ms32), c(
    "odor", "active", "asIs", "stalk-root", "active", "asIs",
    "spore-print-color", "active", "asIs", "f", "predicted", "asIs"
  ))

  # child segment 0
  ms33 <- unlist(xmlToList(model_pmml_3)$MiningModel$Segmentation[[2]]$MiningSchema)
  expect_equal(ms33, c(
    "odor", "active", "asIs", "stalk-root", "active", "asIs",
    "spore-print-color", "active", "asIs", "f", "predicted", "asIs"
  ))
  # child segment 1
  ms34 <- unlist(xmlToList(model_pmml_3)$MiningModel$Segmentation[[5]]$MiningSchema)
  expect_equal(ms34, c(
    "odor", "active", "asIs", "stalk-root", "active", "asIs",
    "spore-print-color", "active", "asIs", "f", "predicted", "asIs"
  ))

  # child segment 2 - the regression model segment
  ms35 <- unlist(xmlToList(model_pmml_3)$MiningModel$Segmentation[[8]]$MiningSchema)
  expect_equal(ms35, c(
    "predictedValueTree0", "active", "continuous", "asIs",
    "predictedValueTree1", "active", "continuous", "asIs"
  ))
})


test_that("invalidValueTreatment attribute is exported correctly for randomForest models", {
  rf_fit <- randomForest(Species ~ ., data = iris, ntree = 3)

  # default invalidValueTreatment arguments
  rf_fit_pmml_1 <- pmml(rf_fit)

  # parent segment
  expect_equal(
    unlist(xmlToList(rf_fit_pmml_1)$MiningModel$MiningSchema),
    c(
      "Species", "predicted", "returnInvalid", "Sepal.Length", "active", "returnInvalid",
      "Sepal.Width", "active", "returnInvalid", "Petal.Length", "active", "returnInvalid",
      "Petal.Width", "active", "returnInvalid"
    )
  )

  # child segment 1
  expect_equal(
    unlist(xmlToList(rf_fit_pmml_1)$MiningModel$Segmentation[[2]]$MiningSchema),
    c(
      "Species", "predicted", "asIs", "Sepal.Length", "active", "asIs",
      "Sepal.Width", "active", "asIs", "Petal.Length", "active", "asIs",
      "Petal.Width", "active", "asIs"
    )
  )

  # child segment 2
  expect_equal(
    unlist(xmlToList(rf_fit_pmml_1)$MiningModel$Segmentation[[5]]$MiningSchema),
    c(
      "Species", "predicted", "asIs", "Sepal.Length", "active", "asIs",
      "Sepal.Width", "active", "asIs", "Petal.Length", "active", "asIs",
      "Petal.Width", "active", "asIs"
    )
  )

  # child segment 3
  expect_equal(
    unlist(xmlToList(rf_fit_pmml_1)$MiningModel$Segmentation[[8]]$MiningSchema),
    c(
      "Species", "predicted", "asIs", "Sepal.Length", "active", "asIs",
      "Sepal.Width", "active", "asIs", "Petal.Length", "active", "asIs",
      "Petal.Width", "active", "asIs"
    )
  )

  # non-default invalidValueTreatment arguments - 1
  rf_fit_pmml_2 <- pmml(rf_fit,
    parent_invalid_value_treatment = "returnInvalid",
    child_invalid_value_treatment = "returnInvalid"
  )

  # parent segment
  expect_equal(
    unlist(xmlToList(rf_fit_pmml_2)$MiningModel$MiningSchema),
    c(
      "Species", "predicted", "returnInvalid", "Sepal.Length", "active", "returnInvalid",
      "Sepal.Width", "active", "returnInvalid", "Petal.Length", "active", "returnInvalid",
      "Petal.Width", "active", "returnInvalid"
    )
  )

  # child segment 1
  expect_equal(
    unlist(xmlToList(rf_fit_pmml_2)$MiningModel$Segmentation[[2]]$MiningSchema),
    c(
      "Species", "predicted", "returnInvalid", "Sepal.Length", "active", "returnInvalid",
      "Sepal.Width", "active", "returnInvalid", "Petal.Length", "active", "returnInvalid",
      "Petal.Width", "active", "returnInvalid"
    )
  )

  # child segment 2
  expect_equal(
    unlist(xmlToList(rf_fit_pmml_2)$MiningModel$Segmentation[[5]]$MiningSchema),
    c(
      "Species", "predicted", "returnInvalid", "Sepal.Length", "active", "returnInvalid",
      "Sepal.Width", "active", "returnInvalid", "Petal.Length", "active", "returnInvalid",
      "Petal.Width", "active", "returnInvalid"
    )
  )

  # child segment 3
  expect_equal(
    unlist(xmlToList(rf_fit_pmml_2)$MiningModel$Segmentation[[8]]$MiningSchema),
    c(
      "Species", "predicted", "returnInvalid", "Sepal.Length", "active", "returnInvalid",
      "Sepal.Width", "active", "returnInvalid", "Petal.Length", "active", "returnInvalid",
      "Petal.Width", "active", "returnInvalid"
    )
  )


  # non-default invalidValueTreatment arguments - 2
  rf_fit_pmml_3 <- pmml(rf_fit,
    parent_invalid_value_treatment = "asIs"
  )

  # parent segment
  expect_equal(
    unlist(xmlToList(rf_fit_pmml_3)$MiningModel$MiningSchema),
    c(
      "Species", "predicted", "asIs", "Sepal.Length", "active", "asIs",
      "Sepal.Width", "active", "asIs", "Petal.Length", "active", "asIs",
      "Petal.Width", "active", "asIs"
    )
  )

  # child segment 1
  expect_equal(
    unlist(xmlToList(rf_fit_pmml_3)$MiningModel$Segmentation[[2]]$MiningSchema),
    c(
      "Species", "predicted", "asIs", "Sepal.Length", "active", "asIs",
      "Sepal.Width", "active", "asIs", "Petal.Length", "active", "asIs",
      "Petal.Width", "active", "asIs"
    )
  )

  # child segment 2
  expect_equal(
    unlist(xmlToList(rf_fit_pmml_3)$MiningModel$Segmentation[[5]]$MiningSchema),
    c(
      "Species", "predicted", "asIs", "Sepal.Length", "active", "asIs",
      "Sepal.Width", "active", "asIs", "Petal.Length", "active", "asIs",
      "Petal.Width", "active", "asIs"
    )
  )

  # child segment 3
  expect_equal(
    unlist(xmlToList(rf_fit_pmml_3)$MiningModel$Segmentation[[8]]$MiningSchema),
    c(
      "Species", "predicted", "asIs", "Sepal.Length", "active", "asIs",
      "Sepal.Width", "active", "asIs", "Petal.Length", "active", "asIs",
      "Petal.Width", "active", "asIs"
    )
  )
})


test_that("error is thrown if invalidValueTreatment argument is incorrect", {
  data(agaricus.train, package = "xgboost")
  train <- agaricus.train
  invisible(capture.output(model_fit_2 <- xgboost(
    data = train$data, label = train$label,
    max_depth = 2, eta = 1, nthread = 2, nrounds = 2, objective = "binary:logistic",
    save_name = tmp_03_save
  )))
  xgb.dump(model_fit_2, tmp_03_dump)

  # default invalidValueTreatment arguments
  model_pmml_5 <- pmml(
    model = model_fit_2, input_feature_names = colnames(train$data),
    output_label_name = "f", output_categories = c("0", "1"),
    xgb_dump_file = tmp_03_dump
  )

  expect_error(
    pmml(
      model = model_fit_2, input_feature_names = colnames(train$data),
      output_label_name = "f", output_categories = c("0", "1"),
      xgb_dump_file = tmp_03_dump,
      parent_invalid_value_treatment = "foobar"
    ),
    "\"foobar\" is not a valid enumeration value for parent_invalid_value_treatment. Use one of the following: returnInvalid, asIs, asMissing."
  )
  expect_error(
    pmml(
      model = model_fit_2, input_feature_names = colnames(train$data),
      output_label_name = "f", output_categories = c("0", "1"),
      xgb_dump_file = tmp_03_dump,
      child_invalid_value_treatment = "asis"
    ),
    "\"asis\" is not a valid enumeration value for child_invalid_value_treatment. Use one of the following: returnInvalid, asIs, asMissing."
  )
})
