m <- SDMtune:::bm_maxent

test_that("The function creates the correct output", {
  expect_s4_class(x <- SDMmodel2MaxEnt(m), "MaxEnt")
  expect_equal(x@presence, .get_presence(m@data))
  expect_equal(x@absence, .get_absence(m@data))
  expect_equal(x@lambdas, m@model@lambdas)
  expect_true(x@hasabsence)
})

test_that("The function raises errors", {
  expect_error(SDMmodel2MaxEnt(SDMtune:::bm_maxnet))
})
