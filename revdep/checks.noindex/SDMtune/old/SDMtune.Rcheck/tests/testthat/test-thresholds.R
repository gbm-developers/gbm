m <- SDMtune:::bm_maxent
test = SDMtune:::t

test_that("The output is correct when test argument is not given", {
  th <- thresholds(m, type = "cloglog")
  expect_named(th, c("Threshold", "Cloglog value", "Fractional predicted area",
                     "Training omission rate"))
  expect_equal(th$Threshold, c("Minimum training presence",
                               "Equal training sensitivity and specificity",
                               "Maximum training sensitivity plus specificity"))
  expect_equal(class(th), "data.frame")
})

test_that("The output is correct when test argument is given", {
  th <- thresholds(m, type = "cloglog", test = test)
  expect_named(th, c("Threshold", "Cloglog value", "Fractional predicted area",
                     "Training omission rate", "Test omission rate",
                     "P-values"))
  expect_equal(th$Threshold, c("Minimum training presence",
                               "Equal training sensitivity and specificity",
                               "Maximum training sensitivity plus specificity",
                               "Equal test sensitivity and specificity",
                               "Maximum test sensitivity plus specificity"))
  expect_equal(class(th), "data.frame")
})
