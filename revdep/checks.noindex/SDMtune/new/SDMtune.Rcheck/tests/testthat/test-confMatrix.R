skip_on_cran()
skip_on_appveyor()

m <- SDMtune:::bm_maxnet
test <- SDMtune:::t
test <- .subset_swd(test, test@pa == 1)
cm <- confMatrix(m, type = "cloglog", th = 0.4)

test_that("There is the correct number of thresholds when passing th", {
  expect_length(cm$th, 1)
})

test_that("The output as the correct column names", {
  expect_named(cm, c("th", "tp", "fp", "fn", "tn"))
})

test_that("The output is correct", {
  # Threshold is correct
  expect_equal(cm$th, 0.4)
  # Sum of tp, fp, fn, tn is equal to sum of presence and background locations
  expect_equal(sum(cm[1, 2:5]), nrow(m@data@data))
  # Sum of tp and fn is equal to number of presence locations
  expect_equal(sum(cm$tp, cm$fn), nrow(m@data@data[m@data@pa == 1, ]))
  # Sum of tn and fp is equal to number of background locations
  expect_equal(sum(cm$tn, cm$fp), nrow(m@data@data[m@data@pa == 0, ]))
  # Correct output with test argument
  # 402 is the number of unique prediction values plus 0 and 1
  expect_equal(nrow(confMatrix(m, test = test, type = "cloglog")), 402)
})

test_that("The thresholds start with 0 and end with 1 when th is not passed", {
  cm <- confMatrix(m, type = "logistic")
  expect_equal(cm$th[1], 0)
  expect_equal(cm$th[nrow(cm)], 1)
})

test_that("Exception is raised", {
  expect_error(confMatrix(SDMtune:::bm_maxent_cv, type = "cloglog"),
               "Function available only for SDMmodel objects.")
})
