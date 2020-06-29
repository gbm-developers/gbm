data(iris)
data(audit)

test_that("ClusteringModel/stats kmeans pmml() output contains MiningSchema, Output, ClusteringField, and Cluster nodes", {
  library(clue)
  fit <- kmeans(iris[, 1:4], 3)
  pmml_fit <- pmml(fit)
  expect_named(names(pmml_fit[[3]]), c("MiningSchema", "Output", "ComparisonMeasure", "ClusteringField", "ClusteringField", "ClusteringField", "ClusteringField", "Cluster", "Cluster", "Cluster"))
})

test_that("GeneralRegressionModel/glmnet pmml() output contains Extension, MiningSchema, Output, ParameterList, CovariateList, PPMatrix, and ParamMatrix nodes", {
  library(glmnet)
  x <- data.matrix(iris[1:4])
  y <- data.matrix(iris[5])
  fit <- cv.glmnet(x, y)
  pmml_fit <- pmml(fit)
  expect_named(names(pmml_fit[[3]]), c("Extension", "MiningSchema", "Output", "ParameterList", "CovariateList", "PPMatrix", "ParamMatrix"))
})

test_that("GeneralRegressionModel/stats pmml() output contains MiningSchema, Output, ParameterList, FactorList, CovariateList, PPMatrix, and ParamMatrix nodes", {
  suppressWarnings(fit <- glm(as.factor(Adjusted) ~ Age + Employment + Education + Marital + Occupation + Income + Sex + Deductions + Hours,
    family = binomial(link = logit), audit
  ))
  pmml_fit <- pmml(fit)
  expect_named(names(pmml_fit[[3]]), c("MiningSchema", "Output", "ParameterList", "FactorList", "CovariateList", "PPMatrix", "ParamMatrix"))
})

test_that("MiningModel/randomForest pmml() output contains MiningSchema, Output, and Segmentation nodes", {
  library(randomForest)
  fit <- randomForest(Species ~ ., data = iris, ntree = 3)
  a <- capture.output(pmml_fit <- pmml(fit))
  expect_named(names(pmml_fit[[3]]), c("MiningSchema", "Output", "Segmentation"))
})

test_that("NaiveBayesModel/e1071 pmml() output contains MiningSchema, Output, BayesInputs, and BayesOutput nodes", {
  library(e1071)
  fit <- naiveBayes(as.factor(Adjusted) ~ Employment + Education + Marital + Occupation + Sex, data = audit)
  pmml_fit <- pmml(fit, predicted_field = "Adjusted")
  expect_named(names(pmml_fit[[3]]), c("MiningSchema", "Output", "BayesInputs", "BayesOutput"))
})

test_that("NeuralNetwork/nnet pmml() output contains MiningSchema, Output, NeuralInputs, NeuralLayer, and NeuralOutputs nodes", {
  library(nnet)
  fit <- nnet(Species ~ ., data = iris, size = 4, trace = F)
  pmml_fit <- pmml(fit)
  expect_named(names(pmml_fit[[3]]), c("MiningSchema", "Output", "NeuralInputs", "NeuralLayer", "NeuralLayer", "NeuralOutputs"))
})

test_that("RegressionModel/stats pmml() output contains MiningSchema, Output, and RegressionTable nodes", {
  fit <- lm(Sepal.Length ~ ., data = iris)
  pmml_fit <- pmml(fit)
  expect_named(names(pmml_fit[[3]]), c("MiningSchema", "Output", "RegressionTable"))
})

test_that("RegressionModel/nnet pmml() output contains MiningSchema, Output, and RegressionTable nodes", {
  library(nnet)
  fit <- multinom(Species ~ ., data = iris, trace = F)
  pmml_fit <- pmml(fit)
  expect_named(names(pmml_fit[[3]]), c("MiningSchema", "Output", "RegressionTable", "RegressionTable", "RegressionTable"))
})

test_that("SupportVectorMachineModel/e1071 pmml() output contains MiningSchema, Output, LocalTransformations, RadialBasisKernelType, VectorDictionary, and SupportVectorMachine nodes", {
  library(e1071)
  fit <- svm(Species ~ ., data = iris)
  pmml_fit <- pmml(fit)
  expect_named(names(pmml_fit[[3]]), c("MiningSchema", "Output", "LocalTransformations", "RadialBasisKernelType", "VectorDictionary", "SupportVectorMachine", "SupportVectorMachine", "SupportVectorMachine"))
})

test_that("SupportVectorMachineModel/kernlab pmml() output contains MiningSchema, Output, LocalTransformations, RadialBasisKernelType, VectorDictionary, and SupportVectorMachine nodes", {
  library(kernlab)
  a <- capture.output(fit <- ksvm(Species ~ ., data = iris, kernel = "rbfdot"))
  pmml_fit <- pmml(fit, data = iris)
  expect_named(names(pmml_fit[[3]]), c("MiningSchema", "Output", "LocalTransformations", "RadialBasisKernelType", "VectorDictionary", "SupportVectorMachine", "SupportVectorMachine", "SupportVectorMachine"))
})

test_that("Transformations pmml() output contains MiningSchema, Output, LocalTransformations, and RegressionTable nodes", {
  dataBox <- xform_wrap(iris)
  dataBox <- xform_min_max(dataBox, "1")
  dataBox <- xform_z_score(dataBox, "1", map_missing_to = 999)
  dataBox <- xform_norm_discrete(dataBox, input_var = "Species")
  dataBox <- xform_function(dataBox, orig_field_name = "Sepal.Width", new_field_name = "a_derived_field", expression = "sqrt(Sepal.Width^2 - 3)")

  fit <- lm(Petal.Width ~ ., data = dataBox$data)
  pmml_fit <- pmml(fit, transform = dataBox)

  expect_named(names(pmml_fit[[3]]), c("MiningSchema", "Output", "LocalTransformations", "RegressionTable"))
})

test_that("Rpart pmml() output contains MiningSchema, Output, and Node nodes", {
  library(rpart)
  fit <- rpart(Species ~ ., data = iris)
  pmml_fit <- pmml(fit)
  expect_named(names(pmml_fit[[3]]), c("MiningSchema", "Output", "Node"))
})
