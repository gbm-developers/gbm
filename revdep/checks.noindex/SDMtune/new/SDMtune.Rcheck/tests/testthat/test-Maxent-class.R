test_that("Show method for Maxent class produces the correct output", {
  m <- SDMtune:::bm_maxent@model
  expect_output(print(m), "Class     : Maxent", fixed = TRUE)
  expect_output(print(m), "Reg       : 1", fixed = TRUE)
  expect_output(print(m), "FCs       : lqph", fixed = TRUE)
  expect_output(print(m), "Iterations: 500", fixed = TRUE)
})
