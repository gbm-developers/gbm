test_that("Show method for Maxnet class produces the correct output", {
  m <- SDMtune:::bm_maxnet@model
  expect_output(print(m), "Class: Maxnet", fixed = TRUE)
  expect_output(print(m), "Reg  : 1", fixed = TRUE)
  expect_output(print(m), "FCs  : lqph", fixed = TRUE)
})
