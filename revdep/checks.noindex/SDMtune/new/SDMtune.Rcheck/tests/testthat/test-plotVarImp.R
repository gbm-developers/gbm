vi <- data.frame(Variable = c("bio1", "bio5", "bio7"),
                 Percent_contribution = c(36.0, 12.4, 27.3),
                 stringsAsFactors = FALSE)

test_that("Output is correct", {
  p <- plotVarImp(vi)
  expect_equal(p$labels$y, "Percent contribution")
  # The function orders the values
  expect_equivalent(p$data$Variable, as.factor(c("bio5", "bio7", "bio1")))
  colnames(vi) <- c("Variable", "Permutation_importance")
})
