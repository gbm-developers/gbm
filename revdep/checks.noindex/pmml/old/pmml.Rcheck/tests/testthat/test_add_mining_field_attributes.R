test_that("add_mining_field_attributes writes MiningField attributes properly", {
  model0 <- lm(Sepal.Length ~ ., data = iris[, -5])
  model <- pmml(model0)

  attributes <- data.frame(
    c("active", 1.1, "asIs"),
    c("active", 2.2, "asIs"),
    c("active", NA, "asMissing"),
    stringsAsFactors = TRUE
  )
  rownames(attributes) <- c(
    "usageType", "missingValueReplacement",
    "invalidValueTreatment"
  )
  colnames(attributes) <- c(
    "Sepal.Width", "Petal.Length",
    "Petal.Width"
  )

  # Although not needed in this first try, necessary to easily
  # add new values later
  for (k in 1:ncol(attributes)) {
    attributes[[k]] <- as.character(attributes[[k]])
  }

  model <- add_mining_field_attributes(model, attributes, namespace = "4_4")

  ms <- xmlToList(model)$RegressionModel$MiningSchema
  ms1 <- ms[[1]]
  ms2 <- ms[[2]]
  ms3 <- ms[[3]]
  ms4 <- ms[[4]]

  expect_equal(length(ms1), 3)
  expect_equal(length(ms2), 4)
  expect_equal(names(ms2)[3], "invalidValueTreatment")
  expect_equal(names(ms2)[4], "missingValueReplacement")
  expect_equal(as.character(ms2)[4], "1.1")
  expect_equal(length(ms3), 4)
  expect_equal(names(ms2)[3], "invalidValueTreatment")
  expect_equal(names(ms3)[4], "missingValueReplacement")
  expect_equal(as.character(ms3)[3], "asIs")
  expect_equal(length(ms4), 3)
  expect_equal(names(ms4)[3], "invalidValueTreatment")
  expect_equal(as.character(ms4)[3], "asMissing")
})
