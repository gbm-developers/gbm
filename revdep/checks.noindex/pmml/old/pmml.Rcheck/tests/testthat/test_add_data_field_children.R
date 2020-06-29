test_that("add_data_field_children works correctly", {
  model0 <- lm(Sepal.Length ~ ., data = iris[, -5])
  model <- pmml(model0)

  mi <- make_intervals(list("openClosed", "openOpen", "closedOpen"), list(NULL, 1, 2), list(1, 2, NULL))
  mv <- make_values(list(1.1, 2.2, 3.3), list(NULL, NULL, NULL), list("valid", NULL, "invalid"))
  model <- add_data_field_children(model, field = "Sepal.Length", intervals = mi, values = mv)

  expect_equal(xmlToList(mi[[1]])[[1]], "openClosed")
  expect_equal(names(xmlToList(mi[[1]]))[2], "rightMargin")

  test <- xmlToList(model)[[2]][[1]][[6]]
  expect_equal(names(test)[1], "value")
  expect_equal(test[[2]], "invalid")
})
