data(iris)

fit <- lm(Sepal.Length ~ ., data = iris[, -5])
fit_pmml <- pmml(fit)
attributes <- data.frame(c("FlowerWidth", 1), c("FlowerLength", 0),
  stringsAsFactors = FALSE
)
rownames(attributes) <- c("displayName", "isCyclic")
colnames(attributes) <- c("Sepal.Width", "Petal.Length")
attributes[] <- lapply(attributes, as.character)

test_that("Error when attribute is NULL", {
  expect_error(add_data_field_attributes(
    xml_model = fit_pmml, attributes = NULL,
    field = "Sepal.Width"
  ), "attribute must be a data.frame or vector.")
})

test_that("Error when field additions are not in model", {
  expect_error(
    add_data_field_attributes(
      fit_pmml,
      c(
        displayName = "FlowerWidth",
        isCyclic = 1
      ), "Foo.Bar"
    ),
    "The following field additions are not in the model: Foo.Bar"
  )


  attributes_2 <- data.frame(c("FlowerWidth", 1), c("FlowerLength", 0),
    stringsAsFactors = FALSE
  )
  rownames(attributes_2) <- c("displayName", "isCyclic")
  colnames(attributes_2) <- c("FooBar", "raBooF")
  expect_error(
    add_data_field_attributes(fit_pmml, attributes_2,
      namespace = "4_4"
    ),
    "The following field additions are not in the model: FooBar, raBooF"
  )
})
