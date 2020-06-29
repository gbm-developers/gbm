context("utility-roc")

mtdt = list(pn = c(0, 0.5, 0.5, 1, 1, 1, 1),
            pp = c(0, 0, 0.25, 0.25, 0.5, 0.75, 1),
            pa = c(0, 1/6, 1/3, 1/2, 2/3, 5/6, 1))

test_that("utility-roc: common ROC curves", {
  expect_equal(roc(mtdt, "rpp", "tpr"),
               list(x = c(0, 1/6, 1/3, 1/2, 2/3, 5/6, 1),
                    y = c(0, 0.25, 0.5, 0.75, 0.75, 1, 1)))
})

test_that("utility-roc: evaluate ROC curve at cutoff", {

  # Cutoff = 0 or 1
  expect_equal(roc(mtdt, "rpp", "tpr", c(0, 1)), c(0, 1))

  # Cutoff is a data point
  expect_equal(roc(mtdt, "rpp", "tpr", 2/3), 0.75)

  # Cutoff is not a data point
  expect_equal(roc(mtdt, "rpp", "tpr", 2/3 + 1/10), 0.9)
})

test_that("utility-roc: cutoff(s) outside of [0, 1]", {
  expect_warning(roc(mtdt, "rpp", "tpr", c(-0.1, 0, 0.1)))
  expect_warning(roc(mtdt, "rpp", "tpr", c(0.9, 1, 1.1)))
})

test_that("utility-roc: invalid 'cdfx' or 'cdfy'", {
  expect_error(roc(mtdt, "fpp", "tpr"))
  expect_error(roc(mtdt, "fpr", "tpp"))
})
