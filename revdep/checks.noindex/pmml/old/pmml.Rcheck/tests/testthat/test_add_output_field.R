data(iris)

test_that("add_output_field works correctly 1", {
  mod <- lm(Sepal.Length ~ ., iris)
  pmod <- pmml(mod)

  onodes0 <- make_output_nodes(
    name = list("OutputField", "OutputField"),
    attributes = list(list(name = "dbl", optype = "continuous"), NULL),
    expression = list("ln(x)", "ln(x/(1-x))")
  )

  test0 <- xmlToList(onodes0[[1]])$Apply
  expect_equal(test0[[2]][[1]], "ln")
})

test_that("add_output_field works correctly  2", {
  mod <- lm(Sepal.Length ~ ., iris)
  pmod <- pmml(mod)

  onodes0 <- make_output_nodes(
    name = list("OutputField", "OutputField"),
    attributes = list(list(name = "dbl", optype = "continuous"), NULL),
    expression = list("ln(x)", "ln(x/(1-x))")
  )

  test0 <- xmlToList(onodes0[[1]])$Apply
  expect_equal(names(test0[[2]])[1], "function")
})

test_that("add_output_field works correctly 3", {
  mod <- lm(Sepal.Length ~ ., iris)
  pmod <- pmml(mod)

  onodes0 <- make_output_nodes(
    name = list("OutputField", "OutputField"),
    attributes = list(list(name = "dbl", optype = "continuous"), NULL),
    expression = list("ln(x)", "ln(x/(1-x))")
  )

  expect_warning(pmod2 <- add_output_field(
    xml_model = pmod, outputNodes = onodes0, at = "End",
    xformText = NULL, nodeName = NULL, attributes = NULL, whichOutput = 1
  ))

  test1 <- xmlToList(pmod2)[[3]][[2]][[3]]$Apply$Apply$Apply
  expect_equal(test1[[3]][[1]], "-")
})

test_that("add_output_field works correctly 4", {
  onodes1 <- make_output_nodes(name = list("OutputField"), attributes = list(name = "name3", dataType = "double", optype = "continuous"))
  expect_equal(length(onodes1), 1)
  expect_equal(names(xmlToList(onodes1[[1]]))[1], "name")
  expect_equal(as.character(xmlToList(onodes1[[1]])[1]), "name3")
  expect_equal(names(xmlToList(onodes1[[1]]))[2], "dataType")
  expect_equal(as.character(xmlToList(onodes1[[1]]))[2], "double")
  expect_equal(names(xmlToList(onodes1[[1]]))[3], "optype")
  expect_equal(as.character(xmlToList(onodes1[[1]]))[3], "continuous")
})
