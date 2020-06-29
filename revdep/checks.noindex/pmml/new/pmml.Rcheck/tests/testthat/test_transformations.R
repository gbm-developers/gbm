data(iris)
data(audit)

tmp_file <- tempfile()
teardown(unlink(tmp_file))

test_that("xform_wrap box$field_data contains specific column names", {
  iris_box <- xform_wrap(iris)
  expect_equal(names(iris_box$field_data), c(
    "type", "dataType", "orig_field_name",
    "sampleMin", "sampleMax", "xformedMin",
    "xformedMax", "centers", "scales",
    "fieldsMap", "transform", "default",
    "missingValue", "xform_function"
  ))
})

test_that("xform_z_score centers and scales for derived fields equal specific values", {
  iris_box <- xform_wrap(iris)
  iris_box <- xform_z_score(iris_box, "1")
  iris_box <- xform_z_score(iris_box, "2")
  expect_equal(iris_box$field_data["derived_Sepal.Length", "transform"], "zxform")
  expect_equal(iris_box$field_data["derived_Sepal.Width", "transform"], "zxform")
  expect_equal(iris_box$field_data["derived_Sepal.Length", "centers"], 5.843333, tolerance = 1e-6)
  expect_equal(iris_box$field_data["derived_Sepal.Length", "scales"], 0.8280661, tolerance = 1e-6)
  expect_equal(iris_box$field_data["derived_Sepal.Width", "centers"], 3.057333, tolerance = 1e-6)
  expect_equal(iris_box$field_data["derived_Sepal.Width", "scales"], 0.4358663, tolerance = 1e-6)
})

test_that("xform_min_max normalizes all data to be between 0 and 1", {
  iris_box <- xform_wrap(iris)
  iris_box <- xform_min_max(iris_box)
  expect_equal(iris_box$field_data["derived_Sepal.Length", "transform"], "minmax")
  expect_equal(iris_box$field_data["derived_Sepal.Width", "transform"], "minmax")
  expect_equal(iris_box$field_data["derived_Petal.Length", "transform"], "minmax")
  expect_equal(iris_box$field_data["derived_Petal.Width", "transform"], "minmax")
  expect_equal(min(iris_box$data$derived_Sepal.Length), 0)
  expect_equal(max(iris_box$data$derived_Sepal.Length), 1)
  expect_equal(min(iris_box$data$derived_Sepal.Width), 0)
  expect_equal(min(iris_box$data$derived_Petal.Length), 0)
  expect_equal(max(iris_box$data$derived_Petal.Length), 1)
  expect_equal(min(iris_box$data$derived_Petal.Width), 0)
  expect_equal(max(iris_box$data$derived_Petal.Width), 1)
})

test_that("xform_norm_discrete box$field_data$fieldsMap contains values of setosa, versicolor, virginica", {
  iris_box <- xform_wrap(iris)
  iris_box <- xform_norm_discrete(iris_box, xform_info = "Species")
  expect_equal(iris_box$field_data["Species_setosa", "transform"], "NormDiscrete")
  expect_equal(iris_box$field_data["Species_versicolor", "transform"], "NormDiscrete")
  expect_equal(iris_box$field_data["Species_virginica", "transform"], "NormDiscrete")
  expect_equal(iris_box$field_data["Species_setosa", "fieldsMap"][[1]], "setosa")
  expect_equal(iris_box$field_data["Species_versicolor", "fieldsMap"][[1]], "versicolor")
  expect_equal(iris_box$field_data["Species_virginica", "fieldsMap"][[1]], "virginica")
})

test_that("rename_wrap_var box$field_data and box$data contain renamed variable", {
  iris_box <- xform_wrap(iris)
  iris_box <- rename_wrap_var(wrap_object = iris_box, xform_info = "column1->SL")
  expect_equal(row.names(iris_box$field_data)[[1]], "SL")
  expect_equal(names(iris_box$data)[[1]], "SL")
})

test_that("xform_discretize produces correct field_data and data values", {
  iris_box <- xform_wrap(iris)
  t <- list()
  m <- data.frame(rbind(
    c("Petal.Length", "dis_pl", "leftInterval", "leftValue", "rightInterval", "rightValue"),
    c("double", "integer", "string", "double", "string", "double"),
    c("0)", 0, "open", NA, "Open", 0),
    c(NA, 1, "closed", 0, "Open", 1),
    c(NA, 2, "closed", 1, "Open", 2),
    c(NA, 3, "closed", 2, "Open", 3),
    c(NA, 4, "closed", 3, "Open", 4),
    c("[4", 5, "closed", 4, "Open", NA)
  ), stringsAsFactors = TRUE)
  t[[1]] <- m
  def <- c(11)
  mis <- c(22)
  iris_box <- xform_discretize(iris_box, xform_info = t, default_value = def, map_missing_to = mis)

  expect_equal(iris_box$field_data["dis_pl", "transform"], "discretize")
  expect_equal(iris_box$field_data["dis_pl", "default"], 11)
  expect_equal(iris_box$field_data["dis_pl", "missingValue"], 22)
  expect_true(iris_box$data$dis_pl[[1]] == 2)
})


test_that("xform_discretize produces correct discretization for a closed left interval", {
  iris_box <- xform_wrap(iris)
  t <- list()
  m <- data.frame(rbind(
    c("Petal.Length", "dis_pl", "leftInterval", "leftValue", "rightInterval", "rightValue"),
    c("double", "integer", "string", "double", "string", "double"),
    c("0)", 0, "open", NA, "Open", 0),
    c(NA, 1, "closed", 0, "Open", 1),
    c(NA, 2, "closed", 1, "Open", 2),
    c(NA, 3, "closed", 2, "Open", 3),
    c(NA, 4, "closed", 3, "Open", 4),
    c("[4", 5, "closed", 4, "Open", NA)
  ))
  t[[1]] <- m
  def <- c(11)
  mis <- c(22)
  iris_box <- xform_discretize(iris_box, xform_info = t, default_value = def, map_missing_to = mis)

  f <- iris_box$data$dis_pl[iris_box$data$Petal.Length == 4]
  expect_equal(as.numeric(levels(f))[f], c(5, 5, 5, 5, 5)) # test that value 4 is transformed to 5
})

test_that("xform_discretize produces correct discretization for a closed right interval", {
  iris_box <- xform_wrap(iris)
  t <- list()
  m <- data.frame(rbind(
    c("Petal.Length", "dis_pl", "leftInterval", "leftValue", "rightInterval", "rightValue"),
    c("double", "integer", "string", "double", "string", "double"),
    c("0]", 0, "open", NA, "Closed", 0),
    c(NA, 1, "open", 0, "Closed", 1),
    c(NA, 2, "open", 1, "Closed", 2),
    c(NA, 3, "open", 2, "Closed", 3),
    c(NA, 4, "open", 3, "Closed", 4),
    c("(4", 5, "open", 4, "Open", NA)
  ))
  t[[1]] <- m
  def <- c(11)
  mis <- c(22)
  iris_box <- xform_discretize(iris_box, xform_info = t, default_value = def, map_missing_to = mis)

  f <- iris_box$data$dis_pl[iris_box$data$Petal.Length == 4]
  expect_equal(as.numeric(levels(f))[f], c(4, 4, 4, 4, 4)) # test that value 4 is transformed to 4
})


test_that("xform_map produces correct mapping for a data point", {
  audit_box <- xform_wrap(audit)
  t <- list()
  m <- data.frame(
    c("Sex", "string", "Male", "Female"), c("Employment", "string", "PSLocal", "PSState"),
    c("d_sex", "integer", 1, 0)
  )
  t[[1]] <- m
  audit_box <- xform_map(audit_box, xform_info = t, default_value = c(3), map_missing_to = 2)

  expect_equal(audit_box$field_data["d_sex", "transform"], "MapValues")
  expect_equal(audit_box$field_data["d_sex", "default"], 3)
  expect_equal(audit_box$field_data["d_sex", "missingValue"], 2)
  expect_true(audit_box$data$d_sex[[1]] == 3)
})

test_that("xform_function works correctly", {
  iris_box <- xform_wrap(iris)
  iris_box <- xform_function(iris_box,
    orig_field_name = "Sepal.Length",
    new_field_name = "Sepal.Length.Transformed",
    expression = "(Sepal.Length^2)/100"
  )
  expect_equal(iris_box$field_data["Sepal.Length.Transformed", "orig_field_name"], "Sepal.Length")
  expect_equal(iris_box$field_data["Sepal.Length.Transformed", "xform_function"], "(Sepal.Length^2)/100")
  expect_equal(iris_box$data$Sepal.Length.Transformed[[1]], 0.2601)
})

test_that("xform_discretize does not give error when 1st column of data matrix is a factor", {
  iris2 <- iris
  iris2[, 6] <- iris2[, 1]
  colnames(iris2)[6] <- "Sepal.Length"
  iris2[, 1] <- iris2[, 5]
  iris2[, 5] <- NULL
  colnames(iris2)[1] <- "Species"
  iris_box <- xform_wrap(iris2)

  t <- list()
  m <- data.frame(rbind(
    c("Petal.Length", "dis_pl", "leftInterval", "leftValue", "rightInterval", "rightValue"),
    c("double", "integer", "string", "double", "string", "double"),
    c(NA, 0, "open", NA, "Open", 0),
    c(NA, 1, "closed", 0, "Closed", 1),
    c(NA, 2, "open", 1, "Closed", 2),
    c(NA, 3, "open", 2, "Open", 3),
    c(NA, 4, "closed", 3, "Open", 4),
    c(NA, 5, "closed", 4, "Open", NA)
  ))
  t[[1]] <- m
  def <- c(11)
  mis <- c(22)

  iris_box <- xform_discretize(iris_box, xform_info = t, default_value = def, map_missing_to = mis)
  expect_equal(iris_box$field_data[6, 11], "discretize")
})


test_that(".init_wrap_params adds NA columns", {
  iris_box <- xform_wrap(iris)

  iris_box$field_data$xformedMax <- NULL
  iris_box$field_data$centers <- NULL
  iris_box$field_data$fieldsMap <- NULL
  iris_box$field_data$transform <- NULL
  iris_box$field_data$default <- NULL
  iris_box$field_data$xform_function <- NULL

  iris_box_updated <- .init_wrap_params(iris_box)

  expect_equal(iris_box_updated$field_data$xformedMax, rep(NA, 5))


  # iris_box <- rename_wrap_var(wrap_object = iris_box,
  #                             xform_info = "column1->SL")
})
