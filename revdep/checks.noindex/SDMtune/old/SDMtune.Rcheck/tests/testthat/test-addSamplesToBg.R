x <- SDMtune:::t

test_that("Add all samples to background", {
  o <- addSamplesToBg(x, all = TRUE)
  expect_s4_class(o, "SWD")
  expect_equal(nrow(o@data), 5800)
  expect_equal(nrow(o@coords), 5800)
  expect_equal(length(o@pa), 5800)
  expect_equivalent(o@data[5401:5800, ], x@data[x@pa == 1, ])
  expect_equivalent(o@coords[5401:5800, ], x@coords[x@pa == 1, ])
})

test_that("Add samples to background only if not duplicates", {
  o <- addSamplesToBg(x)
  expect_s4_class(o, "SWD")
  expect_equal(nrow(o@data), 5607)
  expect_equal(nrow(o@coords), 5607)
  expect_equal(length(o@pa), 5607)
})

test_that("Error are raised", {
  expect_error(addSamplesToBg(x@data), "The function accepts only SWD objects.")
})
