x <- SDMtune:::t

test_that("The function raises error", {
  expect_error(corVar(x@data), "\"bg\" must be an SWD object!")
})

test_that("The function creates the correct output", {
  cm <- corVar(x)
  expect_equal(class(cm), "data.frame")
  expect_named(cm, c("Var1", "Var2", "value"))
  expect_true(abs(min(corVar(x, cor_th = 0.8)$value)) >= 0.8)
  # The output is ordered if order = TRUE
  expect_true(abs(cm$value[1]) >= abs(cm$value[2]))
  cm1 <- corVar(x, order = FALSE)
  expect_false(abs(cm1$value[1]) >= abs(cm1$value[2]))
  # The diagonal is removed if remove_diagonal = TRUE
  expect_true(cm[1, 3] != 1)
  cm1 <- corVar(x, remove_diagonal = FALSE)
  expect_equal(cm1[1, 3], 1)
  # The first 2 columns are not factors
  expect_false(is.factor(cm[, 1]))
  expect_false(is.factor(cm[, 2]))
})
