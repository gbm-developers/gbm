# library(amap)
library(clue)
library(data.table)
library(glmnet)
library(ada)
library(gbm)
library(caret)
library(randomForest)
library(xgboost)
library(Matrix)
library(e1071)
library(neighbr)
library(nnet)
library(kernlab)
library(forecast)

data(iris)
data(audit)
data("WWWusage")
iris_p <- read.csv("iris.csv", stringsAsFactors = TRUE)
audit <- na.omit(audit)
audit_factor <- audit
audit_factor[, 13] <- as.factor(audit_factor[, 13])
elnino <- read.csv("elnino.csv", stringsAsFactors = TRUE)
heart <- read.csv("heart.csv", stringsAsFactors = TRUE)
glm_issue3543_data <- read.csv("glm_issue3543_data.csv", stringsAsFactors = TRUE)
credit_class <- read.csv("credit_class.csv", stringsAsFactors = TRUE)
covtype2 <- read.csv("covtype2.csv", header = TRUE, stringsAsFactors = TRUE)
credit <- read.csv("credit.csv", stringsAsFactors = TRUE)
credit_class_01 <- read.csv("credit_class_01.csv", stringsAsFactors = TRUE)
audit_nor_logical <- na.omit(read.csv("audit_nor_logical.csv", stringsAsFactors = TRUE))
audit_nor <- na.omit(read.csv("audit_nor.csv", stringsAsFactors = TRUE))
audit_nor_fake_logical <- na.omit(read.csv("audit_nor_fake_logical.csv", stringsAsFactors = TRUE))
random_data_small <- read.csv("random_data_small.csv", stringsAsFactors = TRUE)
iris_nor <- read.csv("iris_nor.csv", stringsAsFactors = TRUE)
bank <- na.omit(read.csv("bank.csv", stringsAsFactors = TRUE))
audit_r_build_in <- na.omit(read.csv("audit_r_build_in.csv", stringsAsFactors = TRUE))
insurance <- na.omit(read.csv("insurance.csv", stringsAsFactors = TRUE))
iris_bin <- read.csv("iris_bin.csv", stringsAsFactors = TRUE)
house_votes <- na.omit(read.csv("house_votes_84.csv", stringsAsFactors = TRUE))
iris_mini_dot <- read.csv("iris_mini_dot.csv", stringsAsFactors = TRUE)
petfood <- read.csv("petfood.csv", stringsAsFactors = TRUE)
job_cat <- read.csv("job_cat.csv", stringsAsFactors = TRUE)
job_cat_index <- read.csv("job_cat_index.csv", stringsAsFactors = TRUE)
iris_nor_logical <- read.csv("iris_nor_logical.csv", stringsAsFactors = TRUE)
factor_40k <- read.csv("factor_40k.csv", stringsAsFactors = TRUE)
numeric_10k <- na.omit(read.csv("numeric_10k.csv", stringsAsFactors = TRUE))
factor_10k <- read.csv("factor_10k.csv", stringsAsFactors = TRUE)
numeric_no_na_10k <- read.csv("numeric_no_na_10k.csv", stringsAsFactors = TRUE)


xgb_tmp_01_save <- tempfile()
xgb_tmp_01_dump <- tempfile()

teardown(unlink(c(xgb_tmp_01_save, xgb_tmp_01_dump), recursive = TRUE))

validate_pmml <- function(pmml_doc, schema) {

  # Convert pmml_doc from XMLNode to XMLInternalDocument.
  # Necessary to be able to use xmlSchemaValidate.
  pmml_string <- toString(pmml_doc)
  pmml_parsed <- xmlTreeParse(pmml_string, useInternalNodes = TRUE)

  val_result <- XML::xmlSchemaValidate(schema, pmml_parsed)
  if (length(val_result$errors) == 0) {
    return(0)
  } else {
    paste(unlist(lapply(
      val_result$errors,
      function(x) {
        paste(x$line, x$msg, sep = ": ")
      }
    )), collapse = " -- ")
  }
}

zmz_transform_iris <- function(box_obj) {
  # Apply tranforms to box_obj for iris dataset
  box_obj <- xform_z_score(box_obj, "column1->d1")
  box_obj <- xform_z_score(box_obj, "column2->d2")
  box_obj <- xform_z_score(box_obj, "column3->d3")
  box_obj <- xform_z_score(box_obj, "column4->d4")
  box_obj <- xform_min_max(box_obj, "d1->dd1")
  box_obj <- xform_min_max(box_obj, "d2->dd2")
  box_obj <- xform_min_max(box_obj, "d3->dd3")
  box_obj <- xform_min_max(box_obj, "d4->dd4")
  box_obj <- xform_z_score(box_obj, "dd1->ddd1")
  box_obj <- xform_z_score(box_obj, "dd2->ddd2")
  box_obj <- xform_z_score(box_obj, "dd3->ddd3")
  box_obj <- xform_z_score(box_obj, "dd4->ddd4")
  return(box_obj)
}

zmz_transform_elnino <- function(box_obj) {
  # Apply tranforms to box_obj for elnino dataset
  box_obj <- xform_z_score(box_obj, xform_info = "column1->d1")
  box_obj <- xform_z_score(box_obj, xform_info = "column2->d2")
  box_obj <- xform_z_score(box_obj, xform_info = "column3->d3")
  box_obj <- xform_z_score(box_obj, xform_info = "column4->d4")
  box_obj <- xform_z_score(box_obj, xform_info = "column5->d5")
  box_obj <- xform_z_score(box_obj, xform_info = "column6->d6")
  box_obj <- xform_z_score(box_obj, xform_info = "column7->d7")
  box_obj <- xform_min_max(box_obj, xform_info = "d1->dd1")
  box_obj <- xform_min_max(box_obj, xform_info = "d2->dd2")
  box_obj <- xform_min_max(box_obj, xform_info = "d3->dd3")
  box_obj <- xform_min_max(box_obj, xform_info = "d4->dd4")
  box_obj <- xform_min_max(box_obj, xform_info = "d5->dd5")
  box_obj <- xform_min_max(box_obj, xform_info = "d6->dd6")
  box_obj <- xform_min_max(box_obj, xform_info = "d7->dd7")
  box_obj <- xform_z_score(box_obj, xform_info = "dd1->ddd1")
  box_obj <- xform_z_score(box_obj, xform_info = "dd2->ddd2")
  box_obj <- xform_z_score(box_obj, xform_info = "dd3->ddd3")
  box_obj <- xform_z_score(box_obj, xform_info = "dd4->ddd4")
  box_obj <- xform_z_score(box_obj, xform_info = "dd5->ddd5")
  box_obj <- xform_z_score(box_obj, xform_info = "dd6->ddd6")
  box_obj <- xform_z_score(box_obj, xform_info = "dd7->ddd7")
  return(box_obj)
}

zmz_transform_audit <- function(box_obj) {
  # Apply tranforms to box_obj for audit dataset
  box_obj <- xform_z_score(box_obj, "column2->d_Age")
  box_obj <- xform_z_score(box_obj, "column7->d_Income")
  box_obj <- xform_z_score(box_obj, "column9->d_Deductions")
  box_obj <- xform_z_score(box_obj, "column10->d_Hours")
  box_obj <- xform_min_max(box_obj, "d_Age->dd_Age")
  box_obj <- xform_min_max(box_obj, "d_Income->dd_Income")
  box_obj <- xform_min_max(box_obj, "d_Deductions->dd_Deductions")
  box_obj <- xform_min_max(box_obj, "d_Hours->dd_Hours")
  box_obj <- xform_z_score(box_obj, "dd_Age->ddd_Age")
  box_obj <- xform_z_score(box_obj, "dd_Income->ddd_Income")
  box_obj <- xform_z_score(box_obj, "dd_Deductions->ddd_Deductions")
  box_obj <- xform_z_score(box_obj, "dd_Hours->ddd_Hours")
  return(box_obj)
}


# schema <- XML::xmlSchemaParse("pmml-4-4_xslt_20190830_10.5.0.0.xsd")
# schema <- XML::xmlSchemaParse("pmml-4-4.xsd")
schema <- XML::xmlSchemaParse("pmml-4-4_statespace.xsd")


test_that("TimeSeries/Arima PMML validates against schema", {
  fit <- Arima(WWWusage, order = c(1, 0, 1))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- Arima(WWWusage, order = c(0, 0, 0))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- Arima(WWWusage, order = c(3, 1, 1))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- Arima(JohnsonJohnson, order = c(0, 1, 0), seasonal = c(0, 1, 2))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- Arima(AirPassengers, order = c(0, 1, 1), seasonal = c(0, 1, 1))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- Arima(JohnsonJohnson, order = c(0, 1, 0), seasonal = c(0, 1, 2))
  expect_equal(validate_pmml(pmml(fit, exact_least_squares = TRUE), schema), 0)

  fit <- Arima(AirPassengers, order = c(0, 1, 1), seasonal = c(0, 1, 1))
  expect_equal(validate_pmml(pmml(fit, exact_least_squares = TRUE), schema), 0)

  fit <- Arima(JohnsonJohnson, order = c(2, 1, 3), seasonal = c(0, 1, 2))
  expect_equal(validate_pmml(pmml(fit, exact_least_squares = TRUE), schema), 0)

  fit <- Arima(AirPassengers, order = c(4, 2, 1), seasonal = c(1, 1, 1))
  expect_equal(validate_pmml(pmml(fit, exact_least_squares = TRUE), schema), 0)
})


test_that("AnomalyDetectioneModel/iForest PMML validates against schema", {
  skip_on_cran()
  skip_on_ci()

  library(isofor)

  fit <- iForest(iris, nt = 10, phi = 30)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- iForest(as.matrix(iris[, 1:4]), nt = 10, phi = 30)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  box_obj <- xform_wrap(audit[, -1])
  box_obj <- xform_norm_discrete(box_obj, xform_info = "Sex")
  box_obj <- xform_function(box_obj, orig_field_name = "Age,Hours", new_field_name = "Agrs", expression = "Age/Hours")
  fit <- iForest(box_obj$data[, -c(1, 7, 9, 10)], nt = 5, phi = 420)
  expect_equal(validate_pmml(pmml(fit, transforms = box_obj), schema), 0)
})


test_that("ClusteringModel/stats kmeans PMML validates against schema", {
  fit <- kmeans(audit[, c(2, 7, 9, 10, 12)], 2)
  expect_equal(validate_pmml(pmml(fit), schema), 0)


  fit <- kmeans(iris[, 1:4], 3)
  expect_equal(validate_pmml(pmml(fit), schema), 0)


  box_obj <- xform_wrap(iris)
  box_obj <- zmz_transform_iris(box_obj)
  fit <- kmeans(box_obj$data[, 14:17], 3)
  p_fit <- pmml(fit, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris)
  box_obj <- zmz_transform_iris(box_obj)
  box_obj <- xform_map(box_obj, xform_info = "[Species->d_Species][string->double]", table = "iris_class_table.csv", default_value = "-1", map_missing_to = "1")
  box_obj <- xform_norm_discrete(box_obj, input_var = "Species")
  fit <- kmeans(box_obj$data[, 14:21], 3)
  p_fit <- pmml(fit, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris)
  box_obj <- zmz_transform_iris(box_obj)
  box_obj <- xform_map(box_obj,
    xform_info = "[Species->d_Species][string->double]",
    table = "iris_class_full_name_table.csv", default_value = "-1", map_missing_to = "1"
  )
  box_obj <- xform_norm_discrete(box_obj, input_var = "Species")
  fit <- kmeans(box_obj$data[, 14:21], 3)
  p_fit <- pmml(fit, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)
})


test_that("GeneralRegressionModel/glmnet PMML validates against schema", {
  x <- data.matrix(audit[, c(2, 7, 9:10)])
  y <- data.matrix(audit[, 13])
  fit <- cv.glmnet(x, y)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  x <- data.matrix(iris[1:4])
  y <- data.matrix(iris[5]) # changes string categories to numeric
  fit <- cv.glmnet(x, y)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  x <- data.matrix(elnino[1:6])
  y <- data.matrix(elnino[7])
  fit <- cv.glmnet(x, y, family = "poisson")
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- cv.glmnet(x, y)
  expect_equal(validate_pmml(pmml(fit), schema), 0)


  box_obj <- xform_wrap(elnino)
  rownames(box_obj$field_data)[7] <- "predictedScore"
  box_obj$field_data["predictedScore", "dataType"] <- "numeric"
  names(box_obj$data)[7] <- "predictedScore"
  box_obj <- zmz_transform_elnino(box_obj)
  x <- data.matrix(box_obj$data[1:6])
  y <- data.matrix(box_obj$data[7])
  fit <- cv.glmnet(x, y)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris)
  box_obj <- xform_map(box_obj,
    xform_info = "[Species->d_Species][string->double]",
    table = "iris_class_full_name_table.csv", default_value = "-1", map_missing_to = "1"
  )
  x <- data.matrix(box_obj$data[, c(1:3, 6)])
  y <- data.matrix(box_obj$data[4])
  fit <- cv.glmnet(x, y)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(elnino)
  box_obj <- rename_wrap_var(wrap_object = box_obj, xform_info = "temp->predictedScore")
  box_obj <- zmz_transform_elnino(box_obj)
  x <- data.matrix(box_obj$data[1:6])
  y <- data.matrix(box_obj$data[7])
  fit <- cv.glmnet(x, y)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  x <- data.matrix(iris_p[, c(1:3, 5)])
  y <- data.matrix(iris_p[4])
  box_obj <- xform_wrap(x)
  box_obj <- xform_map(box_obj,
    xform_info = "[class->d_class][string->double]",
    table = "iris_class_index_table.csv", default_value = "-1", map_missing_to = "1"
  )
  box_obj <- xform_norm_discrete(box_obj, input_var = "class")
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_length->dis_pl][double->integer]",
    table = "iris_discretize_pl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_width->dis_pw][double->integer]",
    table = "iris_discretize_pw.csv", map_missing_to = "0", default_value = "1"
  )
  fit <- cv.glmnet(as.matrix(box_obj$data), y)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  x <- data.frame(replicate(20, rnorm(1000)), stringsAsFactors = TRUE)
  y <- rnorm(1000)
  box_obj <- xform_wrap(x)
  box_obj <- xform_min_max(box_obj, xform_info = "column1->d_X1")
  box_obj <- xform_min_max(box_obj, xform_info = "X2->d_X2")
  box_obj <- xform_min_max(box_obj, xform_info = "X3->myDerived_X3[10,20]")
  box_obj <- xform_z_score(box_obj, xform_info = "column4->d_X4")
  box_obj <- xform_z_score(box_obj, xform_info = "X5->d_X5")
  fit <- cv.glmnet(as.matrix(box_obj$data), y)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)
})


test_that("GeneralRegressionModel/stats PMML validates against schema", {
  # Suppress warning: "glm.fit: fitted probabilities numerically 0 or 1 occurred"
  suppressWarnings(fit <- glm(
    formula = as.factor(Adjusted) ~ Age + Employment + Education + Marital + Occupation + Income + Sex + Deductions + Hours,
    family = binomial(link = logit), audit
  ))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- glm(Out ~ ., data = glm_issue3543_data)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  suppressWarnings(fit <- glm(
    formula = as.factor(Adjusted) ~ Age + Employment + Education + Marital + Occupation + Income + Sex + Deductions + Hours,
    family = binomial(link = logit), audit
  ))
  pmml_fit <- pmml(fit)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- glm(formula = target ~ A1 + A2 + A3, family = binomial(link = logit), data = credit_class)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- glm(
    formula = Income ~ Age + Employment + Education + Marital + Occupation + Sex + Hours,
    family = Gamma(link = inverse), audit_nor
  )
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- glm(
    formula = Adjusted ~ Age + Employment + Education + Marital + Occupation + Income + Sex + Deductions + Hours,
    family = gaussian(link = identity), audit
  )
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- glm(formula = as.factor(fbs) ~ ., family = binomial(link = logit), heart)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  suppressWarnings(fit <- glm(
    formula = Adjusted ~ Age + Employment + Education + Marital + Occupation + Income + Sex + Deductions + Hours,
    family = poisson(link = log), audit
  ))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  iris_binom <- iris
  iris_binom$y <- I(iris$Species == "setosa")
  # Suppress warning: "glm.fit: algorithm did not converge"
  suppressWarnings(fit <- glm(y ~ ., data = iris_binom, family = binomial))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  box_obj <- xform_wrap(audit)
  box_obj <- zmz_transform_audit(box_obj)
  suppressWarnings(fit <- glm(
    formula = as.factor(Adjusted) ~ ddd_Age + Employment +
      Education + Marital + Occupation + ddd_Income + Sex + ddd_Deductions + ddd_Hours,
    family = binomial(link = logit), box_obj$data
  ))
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  audit$Adjusted <- as.factor(audit$Adjusted)
  box_obj <- xform_wrap(audit)
  box_obj <- zmz_transform_audit(box_obj)
  box_obj <- xform_norm_discrete(box_obj, input_var = "Employment")
  box_obj <- xform_map(box_obj,
    xform_info = "[Marital-> d_Marital][string->double]",
    table = "audit_marital_table.csv", default_value = "-1", map_missing_to = "1"
  )

  suppressWarnings(fit <- glm(
    formula = Adjusted ~ ddd_Age + ddd_Income + ddd_Deductions +
      ddd_Hours + d_Marital + Employment_Private + Employment_Consultant +
      Employment_SelfEmp + Employment_PSLocal + Employment_PSState +
      Employment_PSFederal + Employment_Volunteer + Sex + Occupation + Education,
    family = binomial(link = logit), box_obj$data
  ))
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)
})


test_that("MiningModel/ada PMML validates against schema", {
  fit <- ada(Adjusted ~ Employment + Education + Hours + Income, iter = 3, audit)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- ada(as.factor(fbs) ~ ., iter = 5, data = heart)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- ada(target ~ ., iter = 11, data = credit_class)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  iris_binom_2 <- iris[iris[, 5] != "setosa", ]
  iris_binom_2[, 5] <- as.factor(levels(iris[, 5])[2:3])[as.numeric(iris[, 5]) - 1]
  fit <- ada(Species ~ ., data = iris_binom_2, iter = 20, nu = 0.9, type = "discrete")
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- ada(as.factor(Adjusted) ~ Employment + Education + Hours + Income, iter = 3, audit)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  box_obj <- xform_wrap(audit)
  box_obj <- zmz_transform_audit(box_obj)
  fit <- ada(as.factor(Adjusted) ~ ddd_Age + ddd_Income + Sex + ddd_Deductions + ddd_Hours, iter = 3, box_obj$data)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)
})


test_that("MiningModel/gbm PMML validates against schema", {
  audit_dat <- audit[, -c(1, 4, 6, 9, 10, 11, 12)]

  fit <- gbm(Adjusted ~ ., data = audit_dat, n.trees = 3, interaction.depth = 4, distribution = "bernoulli")
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- gbm(Adjusted ~ ., data = audit_dat, n.trees = 3, interaction.depth = 4, distribution = "gaussian")
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  covtype2_matrix <- as.matrix(covtype2)
  y0 <- as.vector(covtype2_matrix[, "X3"])
  invisible(capture.output(fit <- gbm.fit(covtype2_matrix[, 1:11], y0,
    distribution = "multinomial", n.trees = 3, interaction.depth = 4
  )))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- gbm(target ~ ., data = credit, n.trees = 4, interaction.depth = 4, distribution = "gaussian")
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- gbm(target ~ ., data = credit_class, n.trees = 5, distribution = "multinomial", interaction.depth = 4)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- gbm(target ~ ., data = credit_class_01, n.trees = 3, interaction.depth = 4, distribution = "bernoulli")
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- gbm(Species ~ ., data = iris, n.trees = 2, interaction.depth = 3, distribution = "multinomial")
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  box_obj <- xform_wrap(iris_p)
  box_obj <- zmz_transform_iris(box_obj)
  fit <- gbm(class ~ ddd1 + ddd2 + ddd3 + ddd4, data = box_obj$data, n.trees = 2, interaction.depth = 3, distribution = "multinomial")
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)
})


test_that("MiningModel/randomForest PMML validates against schema", {
  audit_nor_logical[, "Sex"] <- as.factor(audit_nor_logical[, "Sex"])
  suppressWarnings(fit <- randomForest(Adjusted ~ ., audit_nor_logical[, -1], ntree = 8))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  suppressWarnings(fit <- randomForest(Adjusted ~ ., audit_nor, ntree = 4))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  suppressWarnings(fit <- randomForest(Adjusted ~ ., audit_nor_fake_logical, ntree = 5))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  suppressWarnings(fit <- randomForest(predictedClass ~ ., random_data_small, ntree = 7))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  suppressWarnings(fit <- randomForest(temp ~ ., elnino, ntree = 6))
  expect_equal(validate_pmml(pmml(fit), schema), 0)


  fit <- randomForest(SEPAL_LE ~ ., data = iris_nor, ntree = 9)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  iris_nor_logical <- read.csv("iris_nor_logical.csv", stringsAsFactors = TRUE)
  iris_nor_logical[, 5] <- as.factor(iris_nor_logical[, 5])
  fit <- randomForest(SEPAL_LE ~ ., iris_nor_logical, ntree = 7)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  box_obj <- xform_wrap(iris)
  box_obj <- zmz_transform_iris(box_obj)
  set.seed(123)
  fit <- randomForest(Species ~ Petal.Length + ddd2 + ddd3 + ddd4, box_obj$data, ntree = 7)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- zmz_transform_iris(box_obj)
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_width->dis_sw][double->boolean]",
    table = "iris_discretize_bool_sw.csv", map_missing_to = "0", default_value = "1"
  )
  fit <- randomForest(class ~ petal_length + ddd2 + ddd3 + dis_sw, box_obj$data, ntree = 7)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(audit_factor)
  box_obj <- zmz_transform_audit(box_obj)
  set.seed(14)
  fit <- randomForest(as.factor(Adjusted) ~ ddd_Age + Employment + Education + Marital +
    Occupation + ddd_Income + Sex + ddd_Deductions + ddd_Hours,
  box_obj$data,
  ntree = 7
  )
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(audit_factor[, c("Age", "Marital", "Sex", "Adjusted")])
  box_obj <- xform_norm_discrete(box_obj, xform_info = "Marital")
  box_obj <- xform_norm_discrete(box_obj, xform_info = "Sex")
  audit_box_features <- box_obj$data[!names(box_obj$data) %in% c("Marital", "Sex")]
  fit <- randomForest(Adjusted ~ ., audit_box_features, ntree = 3)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(audit_factor[, c(
    "Age", "Employment", "Income",
    "Deductions", "Marital", "Sex", "Adjusted"
  )])
  box_obj <- xform_norm_discrete(box_obj, xform_info = "Marital")
  box_obj <- xform_norm_discrete(box_obj, xform_info = "Sex")
  box_obj <- xform_function(box_obj,
    orig_field_name = "Age",
    new_field_name = "Age.log",
    expression = "log(Age)"
  )
  box_obj <- xform_function(box_obj,
    orig_field_name = "Income,Age,Deductions",
    new_field_name = "Inc.Ded.Age",
    expression = "(Income-Deductions)/Age"
  )
  audit_box_features <- box_obj$data[!names(box_obj$data) %in% c("Marital", "Sex")]
  fit <- randomForest(Adjusted ~ ., audit_box_features, ntree = 3)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris)
  box_obj <- zmz_transform_iris(box_obj)
  set.seed(335)
  fit <- randomForest(Species ~ ddd1 + ddd2 + ddd3 + ddd4, box_obj$data[1:120, ], ntree = 5)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)
})


test_that("MiningModel/xgboost PMML validates against schema", {
  invisible(capture.output(fit <- xgboost(
    data = as.matrix(iris[, 1:4]), label = as.numeric(iris[, 5]) - 1,
    max_depth = 2, eta = 1, nthread = 2, nrounds = 2, objective = "multi:softprob", num_class = 3,
    save_name = xgb_tmp_01_save
  )))
  xgb.dump(fit, xgb_tmp_01_dump)

  pmml_fit <- pmml(
    model = fit, input_feature_names = colnames(iris[, 1:4]), output_label_name = "Species",
    output_categories = c(1, 2, 3), xgb_dump_file = xgb_tmp_01_dump
  )
  expect_equal(validate_pmml(pmml_fit, schema), 0)

  audit_factor <- audit
  audit_factor[, 13] <- as.factor(audit_factor[, 13])
  invisible(capture.output(fit <- xgboost(
    data = as.matrix(audit_factor[, c(2, 7, 9, 10, 12)]),
    label = as.numeric(audit_factor[, 13]) - 1, max_depth = 2, nrounds = 2,
    objective = "binary:logistic", save_name = xgb_tmp_01_save
  )))
  xgb.dump(fit, xgb_tmp_01_dump)
  pmml_fit <- pmml(fit,
    input_feature_names = colnames(audit_factor[, c(2, 7, 9, 10, 12)]), output_label_name = "Adjusted",
    output_categories = c(0, 1), xgb_dump_file = xgb_tmp_01_dump
  )
  expect_equal(validate_pmml(pmml_fit, schema), 0)


  invisible(capture.output(fit <- xgboost(
    data = as.matrix(audit_factor[, c(2, 7, 9, 10, 12)]),
    label = as.numeric(audit_factor[, 13]) - 1,
    max_depth = 2, nrounds = 2,
    objective = "binary:logistic",
    save_name = xgb_tmp_01_save
  )))
  xgb.dump(fit, xgb_tmp_01_dump)
  pmml_fit <- pmml(fit,
    input_feature_names = colnames(audit_factor[, c(2, 7, 9, 10, 12)]), output_label_name = "Adjusted",
    output_categories = c(0, 1),
    xgb_dump_file = xgb_tmp_01_dump
  )
  expect_equal(validate_pmml(pmml_fit, schema), 0)


  sparse_mat <- as.matrix(sparse.model.matrix(Adjusted ~ . - 1, data = audit[, c("Marital", "Sex", "Adjusted")]))
  invisible(capture.output(fit <- xgboost(
    data = sparse_mat, label = audit[, c("Adjusted")], max_depth = 2,
    eta = 1, nthread = 2, nrounds = 2, objective = "binary:logistic",
    save_name = xgb_tmp_01_save
  )))
  xgb.dump(fit, xgb_tmp_01_dump)
  pmml_fit <- pmml(fit,
    input_feature_names = colnames(sparse_mat), output_label_name = "Adjusted",
    output_categories = c(0, 1), xgb_dump_file = xgb_tmp_01_dump
  )
  expect_equal(validate_pmml(pmml_fit, schema), 0)


  # The next 5 tests check that the naming convention where field name strings are
  # subsets of each other does not cause issues. E.g., V1 is a subset of V11 and V112.
  iris_string_subsets <- iris[1:100, ]
  iris_string_subsets[, 5] <- as.factor(as.character(iris_string_subsets[, 5]))
  colnames(iris_string_subsets) <- c("V11", "V112", "V128", "V22", "V1")
  invisible(capture.output(fit <- xgboost(
    data = as.matrix(iris_string_subsets[, 1:4]),
    label = as.numeric(iris_string_subsets[, 5]) - 1,
    max_depth = 3, eta = 1, nthread = 1, nrounds = 3, objective = "binary:logistic",
    save_name = xgb_tmp_01_save
  )))
  xgb.dump(fit, xgb_tmp_01_dump)
  pmml_fit <- pmml(fit,
    input_feature_names = colnames(as.matrix(iris_string_subsets[, 1:4])),
    output_label_name = "V1",
    output_categories = c(1, 2),
    xgb_dump_file = xgb_tmp_01_dump
  )
  expect_equal(validate_pmml(pmml_fit, schema), 0)


  iris_string_subsets <- iris
  colnames(iris_string_subsets) <- c("V11", "V112", "V128", "V1281", "V1")
  invisible(capture.output(fit <- xgboost(
    data = as.matrix(iris_string_subsets[, 1:4]), label = as.numeric(iris_string_subsets[, 5]) - 1,
    max_depth = 4, eta = 1, nthread = 1, nrounds = 3, num_class = 3,
    objective = "multi:softprob",
    save_name = xgb_tmp_01_save
  )))
  xgb.dump(fit, xgb_tmp_01_dump)
  pmml_fit <- pmml(fit,
    input_feature_names = colnames(as.matrix(iris_string_subsets[, 1:4])), output_label_name = "V1",
    output_categories = c(1, 2, 3), xgb_dump_file = xgb_tmp_01_dump
  )
  expect_equal(validate_pmml(pmml_fit, schema), 0)

  # Use larger number of trees (nrounds) so that some are created with no branches.
  invisible(capture.output(fit <- xgboost(
    data = as.matrix(iris_string_subsets[, 1:4]), label = as.numeric(iris_string_subsets[, 5]) - 1,
    max_depth = 4, eta = 1, nthread = 1, nrounds = 18, num_class = 3,
    objective = "multi:softprob",
    save_name = xgb_tmp_01_save
  )))
  xgb.dump(fit, xgb_tmp_01_dump)
  pmml_fit <- pmml(fit,
    input_feature_names = colnames(as.matrix(iris_string_subsets[, 1:4])), output_label_name = "V1",
    output_categories = c(1, 2, 3), xgb_dump_file = xgb_tmp_01_dump
  )
  expect_equal(validate_pmml(pmml_fit, schema), 0)

  # Multinomial model with one tree each
  invisible(capture.output(fit <- xgboost(
    data = as.matrix(iris_string_subsets[, 1:4]), label = as.numeric(iris_string_subsets[, 5]) - 1,
    max_depth = 4, eta = 1, nthread = 1, nrounds = 1, num_class = 3,
    objective = "multi:softprob",
    save_name = xgb_tmp_01_save
  )))
  xgb.dump(fit, xgb_tmp_01_dump)
  pmml_fit <- pmml(fit,
    input_feature_names = colnames(as.matrix(iris_string_subsets[, 1:4])), output_label_name = "V1",
    output_categories = c(1, 2, 3), xgb_dump_file = xgb_tmp_01_dump
  )
  expect_equal(validate_pmml(pmml_fit, schema), 0)


  iris_matrix <- as.matrix(iris[, 1:4])
  invisible(capture.output(model8 <- xgboost(
    data = iris_matrix, label = as.numeric(iris[, 5]) - 1,
    max_depth = 4, eta = 1, nthread = 1, nrounds = 2, num_class = 3,
    objective = "multi:softmax", save_name = xgb_tmp_01_save
  )))
  xgb.dump(model8, xgb_tmp_01_dump)
  pmml_fit <- pmml(model8,
    input_feature_names = colnames(iris_matrix), output_label_name = "Species",
    output_categories = c(0, 1, 2),
    xgb_dump_file = xgb_tmp_01_dump
  )
  expect_equal(validate_pmml(pmml_fit, schema), 0)


  box_obj <- xform_wrap(audit_factor[, c("Marital", "Sex", "Adjusted")])
  box_obj <- xform_norm_discrete(box_obj, xform_info = "Marital")
  box_obj <- xform_norm_discrete(box_obj, xform_info = "Sex")
  output_vector <- as.numeric(audit_factor$Adjusted) - 1
  audit_box_filt <- box_obj$data[!names(box_obj$data) %in% c("Marital", "Sex", "Adjusted")]
  set.seed(234)
  invisible(capture.output(fit <- xgboost(
    data = as.matrix(audit_box_filt),
    label = output_vector, max_depth = 2,
    eta = 1, nthread = 2, nrounds = 2, objective = "binary:logistic",
    save_name = xgb_tmp_01_save
  )))
  xgb.dump(fit, xgb_tmp_01_dump)
  p_fit <- pmml(fit,
    input_feature_names = colnames(audit_box_filt),
    output_label_name = "Adjusted",
    output_categories = c(0, 1),
    xgb_dump_file = xgb_tmp_01_dump,
    transform = box_obj
  )
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(audit_factor[, c("Marital", "Sex", "Adjusted")])
  box_obj <- xform_norm_discrete(box_obj, xform_info = "Marital", levelSeparator = "_")
  box_obj <- xform_norm_discrete(box_obj, xform_info = "Sex", levelSeparator = "_")
  output_vector <- as.numeric(audit_factor$Adjusted) - 1
  audit_box_filt <- box_obj$data[!names(box_obj$data) %in% c("Marital", "Sex", "Adjusted")]
  set.seed(234)
  invisible(capture.output(fit <- xgboost(
    data = as.matrix(audit_box_filt),
    label = output_vector, max_depth = 2,
    eta = 1, nthread = 2, nrounds = 2, objective = "binary:logistic",
    save_name = xgb_tmp_01_save
  )))
  xgb.dump(fit, xgb_tmp_01_dump)
  p_fit <- pmml(fit,
    input_feature_names = colnames(audit_box_filt),
    output_label_name = "Adjusted",
    output_categories = c(0, 1),
    xgb_dump_file = xgb_tmp_01_dump,
    transform = box_obj
  )
  expect_equal(validate_pmml(p_fit, schema), 0)
})


test_that("NaiveBayesModel/e1071 PMML validates against schema", {
  fit <- naiveBayes(as.factor(Adjusted) ~ Employment + Education + Marital + Occupation + Sex, data = audit_nor)
  expect_equal(validate_pmml(pmml(fit, predicted_field = "Adjusted"), schema), 0)

  fit <- naiveBayes(BANKCARD ~ GENDER + MARITAL_STATUS + PROFESSION + SAVINGS_ACCOUNT + ONLINE_ACCESS + JOINED_ACCOUNTS,
    data = bank
  )
  expect_equal(validate_pmml(pmml(fit, predicted_field = "BANKCARD"), schema), 0)

  fit <- naiveBayes(CLASS ~ ., data = iris_nor)
  expect_equal(validate_pmml(pmml(fit, predicted_field = "CLASS"), schema), 0)

  fit <- naiveBayes(Marital ~ ., data = audit[, c(2:8, 10)])
  expect_equal(validate_pmml(pmml(fit, predicted_field = "Marital"), schema), 0)

  fit <- naiveBayes(Marital ~ ., data = audit_r_build_in[, c(2:8, 10)])
  expect_equal(validate_pmml(pmml(fit, predicted_field = "Marital"), schema), 0)

  fit <- naiveBayes(as.factor(amount_of_claims) ~ gender + domicile, data = insurance)
  expect_equal(validate_pmml(pmml(fit, predicted_field = "amount_of_claims"), schema), 0)

  fit <- naiveBayes(as.factor(amount_of_claims) ~ gender + domicile + no_of_claims, data = insurance)
  expect_equal(validate_pmml(pmml(fit, predicted_field = "amount_of_claims"), schema), 0)

  fit <- naiveBayes(class ~ ., data = iris_bin)
  expect_equal(validate_pmml(pmml(fit, predicted_field = "class"), schema), 0)

  fit <- naiveBayes(target ~ ., data = credit_class)
  expect_equal(validate_pmml(pmml(fit, predicted_field = "target"), schema), 0)




  box_obj <- xform_wrap(iris_p)
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_length->dis_pl][double->integer]",
    table = "iris_discretize_pl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_width->dis_pw][double->integer]",
    table = "iris_discretize_pw.csv", map_missing_to = "0", default_value = "1"
  )
  fit <- naiveBayes(class ~ dis_pl + dis_pw + sepal_length + sepal_width, data = box_obj$data)
  p_fit <- pmml(fit, predicted_field = "class", transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_length->dis_pl][double->string]",
    table = "iris_discretize_pl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_width->dis_pw][double->string]",
    table = "iris_discretize_pw.csv", map_missing_to = "0", default_value = "1"
  )

  fit <- naiveBayes(class ~ dis_pl + dis_pw + sepal_length + sepal_width, data = box_obj$data)
  p_fit <- pmml(fit, predicted_field = "class", transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_length->dis_pl][double->string]",
    table = "iris_discretize_pl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_width->dis_pw][double->string]",
    table = "iris_discretize_pw.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_length->dis_sl][double->string]",
    table = "iris_discretize_sl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_width->dis_sw][double->string]",
    table = "iris_discretize_sw.csv", map_missing_to = "0", default_value = "1"
  )

  fit <- naiveBayes(class ~ dis_pl + dis_pw + dis_sl + dis_sw, data = box_obj$data)
  p_fit <- pmml(fit, predicted_field = "class", transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_length->dis_pl][double->string]",
    table = "iris_discretize_pl.csv"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_width->dis_pw][double->string]",
    table = "iris_discretize_pw.csv"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_length->dis_sl][double->string]",
    table = "iris_discretize_sl.csv"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_width->dis_sw][double->string]",
    table = "iris_discretize_sw.csv"
  )

  fit <- naiveBayes(class ~ dis_pl + dis_pw + dis_sl + dis_sw, data = box_obj$data)
  p_fit <- pmml(fit, predicted_field = "class", transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_length->dis_pl][double->string]",
    table = "iris_discretize_pl.csv"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_width->dis_pw][double->string]",
    table = "iris_discretize_pw.csv"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_length->dis_sl][double->string]",
    table = "iris_discretize_sl.csv"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_width->dis_sw][double->string]",
    table = "iris_discretize_sw.csv"
  )
  fit <- naiveBayes(class ~ dis_pl + dis_pw + dis_sl + dis_sw, data = box_obj$data)

  p_fit <- pmml(fit, predicted_field = "class", transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)
})


test_that("NearestNeighborModel/neighbr PMML validates against schema", {
  iris_train <- iris[1:140, ]
  iris_test <- iris[141:150, -c(4, 5)]
  fit <- knn(
    train_set = iris_train, test_set = iris_test, k = 3, categorical_target = "Species",
    continuous_target = "Petal.Width", comparison_measure = "euclidean"
  )
  expect_equal(validate_pmml(pmml(fit), schema), 0)


  iris_with_id <- iris
  iris_with_id$ID <- c(1:150)
  iris_train <- iris_with_id[1:130, -c(4, 5)]
  iris_test <- iris_with_id[132:150, -c(4, 5, 6)]
  fit <- knn(
    train_set = iris_train, test_set = iris_test, k = 5, comparison_measure = "euclidean",
    return_ranked_neighbors = 3, id = "ID"
  )
  expect_equal(validate_pmml(pmml(fit), schema), 0)


  iris_train <- iris_with_id[1:130, ]
  iris_test <- iris_with_id[132:150, -c(5, 6)]
  fit <- knn(
    train_set = iris_train, test_set = iris_test, k = 5, categorical_target = "Species",
    comparison_measure = "squared_euclidean", return_ranked_neighbors = 4, id = "ID"
  )
  expect_equal(validate_pmml(pmml(fit), schema), 0)


  house_votes_nbr <- house_votes
  feature_names <- names(house_votes_nbr)[!names(house_votes_nbr) %in% c("Class", "ID")]
  for (n in feature_names) {
    levels(house_votes_nbr[, n])[levels(house_votes_nbr[, n]) == "n"] <- 0
    levels(house_votes_nbr[, n])[levels(house_votes_nbr[, n]) == "y"] <- 1
  }
  for (n in feature_names) {
    house_votes_nbr[, n] <- as.numeric(levels(house_votes_nbr[, n]))[house_votes_nbr[, n]]
  }
  house_votes_nbr$ID <- c(1:nrow(house_votes_nbr))

  house_votes_train <- house_votes_nbr[1:100, ]
  house_votes_test <- house_votes_nbr[212:232, !names(house_votes_nbr) %in% c("Class", "ID")]
  fit <- knn(
    train_set = house_votes_train, test_set = house_votes_test, k = 7, categorical_target = "Class",
    comparison_measure = "jaccard", return_ranked_neighbors = 3, id = "ID"
  )
  expect_equal(validate_pmml(pmml(fit), schema), 0)


  house_votes_train <- house_votes_nbr[1:30, ]
  house_votes_test <- house_votes_nbr[105:232, !names(house_votes_nbr) %in% c("Class", "ID")]
  fit <- knn(
    train_set = house_votes_train, test_set = house_votes_test, k = 4, categorical_target = "Class",
    comparison_measure = "simple_matching", return_ranked_neighbors = 4, id = "ID"
  )
  expect_equal(validate_pmml(pmml(fit), schema), 0)


  house_votes_train <- house_votes_nbr[2:90, !names(house_votes_nbr) %in% c("Class")]
  house_votes_test <- house_votes_nbr[195:232, !names(house_votes_nbr) %in% c("Class", "ID")]
  fit <- knn(
    train_set = house_votes_train, test_set = house_votes_test, k = 4, comparison_measure = "tanimoto",
    return_ranked_neighbors = 4, id = "ID"
  )
  expect_equal(validate_pmml(pmml(fit), schema), 0)
})


test_that("NeuralNetwork/nnet PMML validates against schema", {
  audit_nor_factor <- audit_nor
  audit_nor_factor[, 13] <- as.factor(audit_nor[, 13])
  invisible(capture.output(fit <- nnet(Marital ~ ., data = audit_nor_factor[, c(2, 5, 7, 8, 10, 13)], size = 4)))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  invisible(capture.output(fit <- nnet(CLASS ~ ., data = iris_nor, size = 4)))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  invisible(capture.output(fit <- nnet(Adjusted ~ ., data = audit_nor, size = 4)))
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  box_obj <- xform_wrap(iris)
  box_obj <- zmz_transform_iris(box_obj)
  box_obj <- xform_map(box_obj,
    xform_info = "[Species->d_Species][string->double]",
    table = "iris_class_table.csv", default_value = "-1", map_missing_to = "1"
  )
  invisible(capture.output(fit <- nnet(Species ~ ddd1 + ddd2 + ddd3 + ddd4, box_obj$data, size = 5)))
  p_fit <- pmml(fit, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris)
  box_obj <- zmz_transform_iris(box_obj)
  invisible(capture.output(fit <- nnet(Species ~ ddd1 + ddd2 + ddd3 + ddd4, box_obj$data, size = 3)))
  p_fit <- pmml(fit, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)
})


test_that("RegressionModel/nnet PMML validates against schema", {
  fit <- multinom(as.factor(Adjusted) ~ ., data = audit_nor, trace = F)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- multinom(CLASS ~ ., data = iris_nor, trace = F)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  box_obj <- xform_wrap(iris_p)
  box_obj <- zmz_transform_iris(box_obj)
  fit <- multinom(class ~ ddd1 + ddd2 + ddd3 + ddd4, data = box_obj$data, trace = F)
  p_fit <- pmml(fit, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris)
  box_obj <- zmz_transform_iris(box_obj)
  box_obj <- xform_map(box_obj,
    xform_info = "[Species->d_Species][string->double]",
    table = "iris_class_table.csv", default_value = "-1", map_missing_to = "1"
  )
  fit <- multinom(Species ~ ddd1 + ddd2 + ddd3 + ddd4, data = box_obj$data, trace = F)
  p_fit <- pmml(fit, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)
})


test_that("RegressionModel/stats PMML validates against schema", {
  fit <- lm(Sepal.Length ~ ., data = iris)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- lm(temp ~ ., data = elnino)
  expect_equal(validate_pmml(pmml(fit), schema), 0)


  box_obj <- xform_wrap(audit)
  box_obj <- xform_map(box_obj,
    xform_info = "[Employment,Education,Sex-> d_E]",
    table = "audit_3to1_table.csv", default_value = "X", map_missing_to = "Y"
  )
  fit <- lm(Adjusted ~ d_E + Income + Hours, data = box_obj$data)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- zmz_transform_iris(box_obj)
  box_obj <- xform_norm_discrete(box_obj, input_var = "class")
  box_obj <- xform_map(box_obj,
    xform_info = "[class->d_class][string->double]",
    table = "iris_p_class_table.csv", default_value = "-1", map_missing_to = "1"
  )
  fit <- lm(sepal_width ~ ddd1 + ddd2 + ddd3 + d_class + class_Iris_setosa +
    class_Iris_versicolor + class_Iris_virginica, box_obj$data)
  p_fit <- pmml(fit, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_length->dis_pl][double->string]",
    table = "iris_discretize_pl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_width->dis_pw][double->integer]",
    table = "iris_discretize_pw.csv", map_missing_to = "0", default_value = "1"
  )

  fit <- lm(sepal_width ~ dis_pl + dis_pw + class, data = box_obj$data)
  p_fit <- pmml(fit, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_length->dis_pl][double->string]",
    table = "iris_discretize_pl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_width->dis_pw][double->integer]",
    table = "iris_discretize_pw.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_norm_discrete(box_obj, input_var = "class")
  box_obj <- xform_map(box_obj,
    xform_info = "[class->d_class][string->double]",
    table = "iris_p_class_table.csv", default_value = "-1", map_missing_to = "1"
  )
  fit <- lm(sepal_width ~ dis_pl + dis_pw + d_class, data = box_obj$data)
  p_fit <- pmml(fit, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)
})



test_that("AnomalyDetectionModel/e1071 one-classification PMML validates against schema", {
  fit <- svm(iris[, 1:3], y = NULL, type = "one-classification", scale = TRUE)
  expect_equal(validate_pmml(pmml(fit, dataset = iris[, 1:3], model_name = "radial_iris_ocsvm"), schema), 0)

  fit <- svm(iris[, 1:4], y = NULL, type = "one-classification", nu = 0.10, scale = FALSE, kernel = "linear")
  expect_equal(validate_pmml(pmml(fit, dataset = iris[, 1:4]), schema), 0)

  fit <- svm(iris[, 1:4], y = NULL, type = "one-classification", nu = 0.11, scale = TRUE, kernel = "polynomial")
  expect_equal(validate_pmml(pmml(fit, dataset = iris[, 1:4]), schema), 0)

  fit <- svm(iris[, 1:4], y = NULL, type = "one-classification", nu = 0.21, kernel = "sigmoid")
  expect_equal(validate_pmml(pmml(fit, dataset = iris[, 1:4]), schema), 0)

  iris_y <- as.numeric(iris$Species == "setosa")
  fit <- svm(iris[, 1:4], y = iris_y, type = "one-classification", nu = 0.15, kernel = "sigmoid")
  expect_equal(validate_pmml(pmml(fit, dataset = iris[, 1:4]), schema), 0)

  fit <- svm(audit[100:400, c("Income", "Deductions")],
    y = NULL, type = "one-classification",
    nu = 0.10, scale = TRUE, kernel = "linear"
  )
  expect_equal(validate_pmml(pmml(fit, dataset = audit[, c("Income", "Deductions")]), schema), 0)

  audit_numeric <- audit[1:500, c("Age", "Income", "Deductions", "Hours", "Adjustment", "Adjusted")]
  audit_numeric$Age <- as.numeric(audit_numeric$Age)
  audit_numeric$Hours <- as.numeric(audit_numeric$Hours)
  audit_numeric$Adjustment <- as.numeric(audit_numeric$Adjustment)
  audit_numeric$Adjusted <- as.numeric(audit_numeric$Adjusted)
  fit <- svm(audit_numeric, y = NULL, type = "one-classification", nu = 0.10, scale = FALSE, kernel = "radial")
  expect_equal(validate_pmml(pmml(fit, dataset = audit_numeric), schema), 0)

  audit_numeric <- audit[600:900, c("Age", "Income", "Deductions", "Hours", "Adjustment", "Adjusted")]
  audit_numeric$Age <- as.numeric(audit_numeric$Age)
  audit_numeric$Hours <- as.numeric(audit_numeric$Hours)
  audit_numeric$Adjustment <- as.numeric(audit_numeric$Adjustment)
  audit_numeric$Adjusted <- as.numeric(audit_numeric$Adjusted)
  fit <- svm(audit_numeric, y = NULL, type = "one-classification", nu = 0.10, scale = FALSE, kernel = "radial")
  expect_equal(validate_pmml(pmml(fit, dataset = audit_numeric), schema), 0)
})


test_that("SupportVectorMachineModel/e1071 PMML validates against schema", {
  fit <- svm(Petal.Width ~ ., data = iris[, 1:4], kernel = "linear")
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- svm(Adjusted ~ Age + Income + Hours, data = audit[1:900, ])
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- svm(Sex ~ ., data = audit[200:700, 2:9], scale = FALSE)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  audit_logical <- audit[1:800, c(2, 8, 13)]
  audit_logical$Adjusted <- as.logical(audit_logical$Adjusted)
  fit <- svm(Sex ~ ., data = audit_logical, scale = FALSE)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- svm(as.factor(Adjusted) ~ Age + Income + Deductions + Hours, data = audit[1:800, ])
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- svm(as.factor(Adjusted) ~ ., data = audit[1:700, ])
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- svm(Marital ~ Income + Deductions, data = audit[1:700, ], kernel = "polynomial")
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- svm(Species ~ ., data = iris)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- svm(sepal_length ~ ., data = iris_mini_dot)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- svm(Species ~ ., data = iris, scale = FALSE, probability = TRUE)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- svm(Species ~ ., data = iris, kernel = "linear")
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- svm(Species ~ ., data = iris, kernel = "polynomial")
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- svm(Species ~ ., data = iris, kernel = "sigmoid")
  expect_equal(validate_pmml(pmml(fit), schema), 0)


  box_obj <- xform_wrap(iris[, 1:4])
  fit <- svm(box_obj$data, y = NULL, type = "one-classification")
  p_fit <- pmml(fit, dataset = iris[, 1:4], transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris[, 1:4])
  box_obj <- xform_z_score(box_obj)
  fit <- svm(box_obj$data[, 5:8], y = NULL, type = "one-classification")
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris[, 1:4])
  box_obj <- zmz_transform_iris(box_obj)
  fit <- svm(box_obj$data[, 13:16], y = NULL, type = "one-classification")
  p_fit <- pmml(fit, dataset = box_obj$data[, 13:16], transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p[, 1:4])
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_length->dis_pl][double->integer]",
    table = "iris_discretize_pl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_width->dis_pw][double->integer]",
    table = "iris_discretize_pw.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_length->dis_sl][double->integer]",
    table = "iris_discretize_sl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_width->dis_sw][double->integer]",
    table = "iris_discretize_sw.csv", map_missing_to = "0", default_value = "1"
  )
  suppressWarnings(fit <- svm(box_obj$data[, 5:8], y = NULL, type = "one-classification"))
  p_fit <- pmml(fit, dataset = box_obj$data[, 5:8], transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(audit)
  box_obj <- zmz_transform_audit(box_obj)
  box_obj <- xform_norm_discrete(box_obj, input_var = "Employment")
  box_obj <- xform_map(box_obj,
    xform_info = "[Marital-> d_Marital][string->double]",
    table = "audit_marital_table.csv",
    default_value = "-1", map_missing_to = "1"
  )
  fit <- svm(box_obj$data[, c(22, 23, 25)],
    y = NULL, type = "one-classification", nu = 0.10,
    scale = TRUE, kernel = "linear"
  )
  p_fit <- pmml(fit, dataset = box_obj$data[, c(22, 23, 25)], transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(audit[, c("Income", "Deductions")])
  fit <- svm(box_obj$data, y = NULL, type = "one-classification")
  p_fit <- pmml(fit, dataset = box_obj$data, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_length->dis_pl][double->integer]",
    table = "iris_discretize_pl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_width->dis_pw][double->integer]",
    table = "iris_discretize_pw.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_length->dis_sl][double->integer]",
    table = "iris_discretize_sl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_width->dis_sw][double->integer]",
    table = "iris_discretize_sw.csv", map_missing_to = "0", default_value = "1"
  )
  suppressWarnings(fit <- svm(class ~ dis_pl + dis_pw + dis_sl + dis_sw, data = box_obj$data))
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(audit)
  box_obj <- zmz_transform_audit(box_obj)
  box_obj <- xform_norm_discrete(box_obj, input_var = "Employment")
  box_obj <- xform_map(box_obj,
    xform_info = "[Marital-> d_Marital][string->double]",
    table = "audit_marital_table.csv", default_value = "-1", map_missing_to = "1"
  )
  fit <- svm(Adjusted ~ ddd_Age + ddd_Income + ddd_Hours, data = box_obj$data)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(audit_factor)
  box_obj <- zmz_transform_audit(box_obj)
  box_obj <- xform_norm_discrete(box_obj, input_var = "Employment")
  box_obj <- xform_map(box_obj,
    xform_info = "[Marital-> d_Marital][string->double]",
    table = "audit_marital_table.csv", default_value = "-1", map_missing_to = "1"
  )
  fit <- svm(Adjusted ~ ., data = box_obj$data[, -c(1, 2, 7, 9, 10, 3, 5)])
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  fit <- svm(as.factor(Adjusted) ~ ddd_Age + ddd_Income + ddd_Deductions + ddd_Hours, data = box_obj$data)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_length->dis_pl][double->integer]",
    table = "iris_discretize_pl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_width->dis_pw][double->integer]",
    table = "iris_discretize_pw.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_length->dis_sl][double->integer]",
    table = "iris_discretize_sl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_width->dis_sw][double->integer]",
    table = "iris_discretize_sw.csv", map_missing_to = "0", default_value = "1"
  )
  suppressWarnings(fit <- svm(class ~ dis_pl + dis_pw + dis_sl + dis_sw, data = box_obj$data))
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- zmz_transform_iris(box_obj)
  fit <- svm(class ~ ddd1 + ddd2 + ddd3 + ddd4, data = box_obj$data)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- xform_z_score(box_obj)
  fit <- svm(class ~ derived_petal_length + derived_petal_width + derived_sepal_length + derived_sepal_width,
    data = box_obj$data
  )
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- zmz_transform_iris(box_obj)
  box_obj <- xform_map(box_obj,
    xform_info = "[class->d_class][string->double]",
    table = "iris_p_class_table.csv", default_value = "-1", map_missing_to = "1"
  )
  box_obj <- xform_norm_discrete(box_obj, input_var = "class")
  fit <- svm(class ~ ddd1 + ddd2 + ddd3 + ddd4, data = box_obj$data)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_length->dis_pl][double->integer]",
    table = "iris_discretize_pl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_width->dis_pw][double->integer]",
    table = "iris_discretize_pw.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_length->dis_sl][double->integer]",
    table = "iris_discretize_sl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_width->dis_sw][double->integer]",
    table = "iris_discretize_sw.csv", map_missing_to = "0", default_value = "1"
  )
  suppressWarnings(fit <- svm(class ~ dis_pl + dis_pw + dis_sl + dis_sw, data = box_obj$data))
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)
})


test_that("SupportVectorMachineModel/kernlab PMML validates against schema", {
  fit <- ksvm(target ~ ., data = credit, kernel = "rbfdot")
  expect_equal(validate_pmml(pmml(fit, data = credit), schema), 0)

  fit <- ksvm(CLASS ~ ., data = iris_nor, kernel = "rbfdot")
  expect_equal(validate_pmml(pmml(fit, data = iris_nor), schema), 0)

  invisible(capture.output(fit <- ksvm(CLASS ~ ., data = iris_nor, kernel = "vanilladot")))
  expect_equal(validate_pmml(pmml(fit, data = iris_nor), schema), 0)

  fit <- ksvm(Adjusted ~ ., data = audit[1:900, ], kernel = "rbfdot")
  expect_equal(validate_pmml(pmml(fit, data = audit[1:900, ]), schema), 0)

  fit <- ksvm(as.factor(purchase) ~ ., data = petfood)
  expect_equal(validate_pmml(pmml(fit, data = petfood), schema), 0)

  fit <- ksvm(as.factor(Adjusted) ~ ., data = audit[1:900, ], kernel = "rbfdot")
  expect_equal(validate_pmml(pmml(fit, data = audit[1:900, ]), schema), 0)

  fit <- ksvm(PRE_1 ~ ., data = job_cat, kernel = "rbfdot")
  expect_equal(validate_pmml(pmml(fit, data = job_cat), schema), 0)

  fit <- ksvm(as.factor(PRE_1) ~ ., data = job_cat_index, kernel = "rbfdot")
  expect_equal(validate_pmml(pmml(fit, data = job_cat_index), schema), 0)


  box_obj <- xform_wrap(iris)
  box_obj <- zmz_transform_iris(box_obj)
  fit <- ksvm(Species ~ ddd1 + ddd2 + ddd3 + ddd4, data = box_obj$data)
  p_fit <- pmml(fit, dataset = box_obj$data, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- zmz_transform_iris(box_obj)
  box_obj <- xform_map(box_obj,
    xform_info = "[class->d_class][string->double]",
    table = "iris_p_class_table.csv", default_value = "-1", map_missing_to = "1"
  )
  box_obj <- xform_norm_discrete(box_obj, xform_info = "class")
  fit <- ksvm(class ~ ddd1 + ddd2 + ddd3 + ddd4, data = box_obj$data)
  p_fit <- pmml(fit, dataset = box_obj$data, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- zmz_transform_iris(box_obj)
  box_obj <- xform_map(box_obj,
    xform_info = "[class->d_class][string->double]",
    table = "iris_p_class_table.csv", default_value = "-1", map_missing_to = "1"
  )
  box_obj <- xform_norm_discrete(box_obj, input_var = "class")
  fit <- ksvm(class ~ ddd1 + ddd2 + ddd3 + ddd4, data = box_obj$data)
  p_fit <- pmml(fit, dataset = box_obj$data, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(audit_factor)
  box_obj <- zmz_transform_audit(box_obj)
  box_obj <- xform_norm_discrete(box_obj, input_var = "Employment")
  box_obj <- xform_map(box_obj,
    xform_info = "[Marital-> d_Marital][string->double]",
    table = "audit_marital_table.csv", default_value = "-1", map_missing_to = "1"
  )
  fit <- ksvm(Adjusted ~ ddd_Age + ddd_Income + ddd_Deductions +
    ddd_Hours + d_Marital + Employment_Private + Employment_Consultant +
    Employment_SelfEmp + Employment_PSLocal + Employment_PSState +
    Employment_PSFederal + Employment_Volunteer + Sex + Occupation +
    Education, data = box_obj$data)

  p_fit <- pmml(fit, dataset = box_obj$data, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)
})


test_that("TreeModel/rpart PMML validates against schema", {
  fit <- rpart(as.factor(Adjusted) ~ Employment + Education + Marital + Occupation + Sex, data = audit_nor)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- rpart(temp ~ ., data = elnino)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- rpart(CLASS ~ ., data = iris_nor)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- rpart(fbs ~ ., data = heart)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  fit <- rpart(as.factor(fbs) ~ ., data = heart)
  expect_equal(validate_pmml(pmml(fit), schema), 0)

  box_obj <- xform_wrap(iris_p)
  box_obj <- zmz_transform_iris(box_obj)
  fit <- rpart(class ~ ddd1 + ddd2 + ddd3 + ddd4, data = box_obj$data)
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)
})


test_that("Transformations PMML validates against schema", {
  skip_on_cran()
  skip_on_ci()

  box_obj <- xform_wrap(iris_p)
  box_obj <- xform_function(box_obj,
    orig_field_name = "sepal_length",
    new_field_name = "a_derived_field",
    expression = "sqrt(sepal_length^2 + 3)"
  )
  box_obj <- xform_function(box_obj,
    orig_field_name = list("sepal_length, sepal_width"),
    new_field_name = "two_field_formula",
    expression = "sepal_length * sepal_width"
  )
  fit <- lm(petal_width ~ ., data = box_obj$data)
  p_fit <- pmml(fit, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- xform_min_max(box_obj, "1")
  box_obj <- xform_z_score(box_obj, "1", map_missing_to = 999)
  box_obj <- xform_norm_discrete(box_obj, input_var = "class")
  box_obj <- xform_function(box_obj,
    orig_field_name = "sepal_width",
    new_field_name = "a_derived_field",
    expression = "sqrt(sepal_width^2 - 3)"
  )
  fit <- lm(petal_width ~ ., data = box_obj$data)
  p_fit <- pmml(fit, transform = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(iris_p)
  box_obj <- xform_discretize(box_obj,
    xform_info = "[petal_width->dis_pw][double->string]",
    table = "iris_discretize_pw.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_length->dis_sl][double->string]",
    table = "iris_discretize_sl.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_discretize(box_obj,
    xform_info = "[sepal_width->dis_sw][double->string]",
    table = "iris_discretize_sw.csv", map_missing_to = "0", default_value = "1"
  )
  box_obj <- xform_map(box_obj,
    xform_info = "[class->d_class][string->double]",
    table = "iris_p_class_table.csv", default_value = "-1", map_missing_to = "1"
  )
  fit <- lm(petal_length ~ ., data = box_obj$data[, -c(2, 3, 4, 5, 7)])
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  box_obj <- xform_wrap(audit_factor)
  box_obj <- zmz_transform_audit(box_obj)
  box_obj <- xform_norm_discrete(box_obj, input_var = "Employment")
  box_obj <- xform_map(box_obj,
    xform_info = "[Marital-> d_Marital][string->double]",
    table = "audit_marital_table.csv", default_value = "-1", map_missing_to = "1"
  )
  fit <- rpart(Adjusted ~ ., data = box_obj$data[, -1])
  p_fit <- pmml(fit, transforms = box_obj)
  expect_equal(validate_pmml(p_fit, schema), 0)


  factor_40k_box <- xform_wrap(factor_40k)
  factor_40k_box <- xform_norm_discrete(factor_40k_box, xform_info = "CateA")
  factor_40k_box <- xform_norm_discrete(factor_40k_box, xform_info = "CateB")
  fit <- rpart(letter ~ ., data = factor_40k_box$data[, -c(2, 3)])
  p_fit <- pmml(fit, transforms = factor_40k_box)
  expect_equal(validate_pmml(p_fit, schema), 0)


  numeric_10k_box <- xform_wrap(numeric_10k)
  numeric_10k_box <- xform_min_max(numeric_10k_box,
    xform_info = "var_10->d_var_10", map_missing_to = "0"
  )
  numeric_10k_box <- xform_min_max(numeric_10k_box,
    xform_info = "var_11->d_var_11", map_missing_to = "0"
  )
  numeric_10k_box <- xform_min_max(numeric_10k_box,
    xform_info = "var_12->d_var_12", map_missing_to = "0"
  )
  numeric_10k_box <- xform_min_max(numeric_10k_box,
    xform_info = "var_13->d_var_13", map_missing_to = "0"
  )
  fit <- lm(var_14 ~ ., data = numeric_10k_box$data)
  p_fit <- pmml(fit, transforms = numeric_10k_box)
  expect_equal(validate_pmml(p_fit, schema), 0)


  numeric_10k_box <- xform_wrap(numeric_10k)
  numeric_10k_box <- xform_z_score(numeric_10k_box, xform_info = "var_0->d_var_0", map_missing_to = "0")
  numeric_10k_box <- xform_z_score(numeric_10k_box, xform_info = "var_1->d_var_1", map_missing_to = "0")
  numeric_10k_box <- xform_z_score(numeric_10k_box, xform_info = "var_2->d_var_2", map_missing_to = "0")
  numeric_10k_box <- xform_z_score(numeric_10k_box, xform_info = "var_3->d_var_3", map_missing_to = "0")
  fit <- lm(var_14 ~ ., data = numeric_10k_box$data)
  p_fit <- pmml(fit, transforms = numeric_10k_box)
  expect_equal(validate_pmml(p_fit, schema), 0)


  factor_10k_box <- xform_wrap(factor_10k)
  factor_10k_box <- xform_norm_discrete(factor_10k_box, input_var = "CateA")
  factor_10k_box <- xform_norm_discrete(factor_10k_box, input_var = "CateB")
  fit <- rpart(letter ~ ., data = factor_10k_box$data[, -c(2, 3)])
  p_fit <- pmml(fit, transforms = factor_10k_box)
  expect_equal(validate_pmml(p_fit, schema), 0)


  a <- which(factor_10k[, 1] == "A")
  b <- which(factor_10k[, 1] == "B")
  y <- which(factor_10k[, 1] == "Y")
  z <- which(factor_10k[, 1] == "Z")
  factor_10k_smp <- factor_10k[sample(c(a, b, y, z), length(c(a, b, y, z))), ]
  factor_10k_smp[, 1] <- as.character(factor_10k_smp[, 1])
  levels(factor_10k_smp[, 1]) <- c("A", "B", "Y", "Z")
  factor_10k_smp[, 1] <- as.factor(factor_10k_smp[, 1])
  factor_10k_box <- xform_wrap(factor_10k_smp)
  factor_10k_box <- xform_map(factor_10k_box,
    xform_info = "[letter,CateA->d_CateB][string,string->string]",
    table = "map_factor_400.csv", default_value = "-1", map_missing_to = "1"
  )
  fit <- rpart(letter ~ ., data = factor_10k_box$data[, -2])
  p_fit <- pmml(fit, transforms = factor_10k_box)
  expect_equal(validate_pmml(p_fit, schema), 0)


  numeric_no_na_10k_box <- xform_wrap(numeric_no_na_10k)
  numeric_no_na_10k_box <- xform_discretize(numeric_no_na_10k_box,
    xform_info = "[var_0->d_var_0][double->integer]",
    table = "numeric_discretize_var.csv",
    map_missing_to = "0", default_value = "1"
  )
  numeric_no_na_10k_box <- xform_discretize(numeric_no_na_10k_box,
    xform_info = "[var_1->d_var_1][double->integer]",
    table = "numeric_discretize_var.csv",
    map_missing_to = "0", default_value = "1"
  )
  numeric_no_na_10k_box <- xform_discretize(numeric_no_na_10k_box,
    xform_info = "[var_2->d_var_2][double->integer]",
    table = "numeric_discretize_var.csv",
    map_missing_to = "0", default_value = "1"
  )
  numeric_no_na_10k_box <- xform_discretize(numeric_no_na_10k_box,
    xform_info = "[var_3->d_var_3][double->integer]",
    table = "numeric_discretize_var.csv",
    map_missing_to = "0", default_value = "1"
  )
  fit <- lm(var_14 ~ ., data = numeric_no_na_10k_box$data[1:600, ])
  p_fit <- pmml(fit, transforms = numeric_no_na_10k_box)
  expect_equal(validate_pmml(p_fit, schema), 0)


  numeric_no_na_10k_box <- xform_wrap(numeric_no_na_10k)
  numeric_no_na_10k_box <- xform_min_max(numeric_no_na_10k_box,
    xform_info = "var_0->d_var_0", map_missing_to = "0"
  )
  numeric_no_na_10k_box <- xform_min_max(numeric_no_na_10k_box,
    xform_info = "var_1->d_var_1", map_missing_to = "0"
  )
  numeric_no_na_10k_box <- xform_min_max(numeric_no_na_10k_box,
    xform_info = "var_2->d_var_2", map_missing_to = "0"
  )
  numeric_no_na_10k_box <- xform_min_max(numeric_no_na_10k_box,
    xform_info = "var_3->d_var_3", map_missing_to = "0"
  )
  fit <- lm(var_14 ~ ., data = numeric_no_na_10k_box$data)
  p_fit <- pmml(fit, transforms = numeric_no_na_10k_box)
  expect_equal(validate_pmml(p_fit, schema), 0)


  numeric_no_na_10k_box <- xform_wrap(numeric_no_na_10k)
  numeric_no_na_10k_box <- xform_z_score(numeric_no_na_10k_box,
    xform_info = "var_0->d_var_0"
  )
  numeric_no_na_10k_box <- xform_z_score(numeric_no_na_10k_box,
    xform_info = "var_1->d_var_1", map_missing_to = "0"
  )
  numeric_no_na_10k_box <- xform_z_score(numeric_no_na_10k_box,
    xform_info = "var_2->d_var_2", map_missing_to = "0"
  )
  numeric_no_na_10k_box <- xform_z_score(numeric_no_na_10k_box,
    xform_info = "var_3->d_var_3", map_missing_to = "0"
  )
  fit <- lm(var_14 ~ ., data = numeric_no_na_10k_box$data)
  p_fit <- pmml(fit, transforms = numeric_no_na_10k_box)
  expect_equal(validate_pmml(p_fit, schema), 0)


  numeric_10k_box <- xform_wrap(numeric_10k)
  numeric_10k_box <- xform_discretize(numeric_10k_box,
    xform_info = "[var_0->d_var_0][double->integer]",
    table = "numeric_discretize_var.csv", map_missing_to = "0", default_value = "1"
  )
  numeric_10k_box <- xform_discretize(numeric_10k_box,
    xform_info = "[var_1->d_var_1][double->integer]",
    table = "numeric_discretize_var.csv", map_missing_to = "0", default_value = "1"
  )
  numeric_10k_box <- xform_discretize(numeric_10k_box,
    xform_info = "[var_2->d_var_2][double->integer]",
    table = "numeric_discretize_var.csv", map_missing_to = "0", default_value = "1"
  )
  numeric_10k_box <- xform_discretize(numeric_10k_box,
    xform_info = "[var_3->d_var_3][double->integer]",
    table = "numeric_discretize_var.csv", map_missing_to = "0", default_value = "1"
  )
  fit <- lm(var_14 ~ ., data = numeric_10k_box$data[1:500, ])
  p_fit <- pmml(fit, transforms = numeric_10k_box)
  expect_equal(validate_pmml(p_fit, schema), 0)
})
