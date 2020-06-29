t <- SDMtune:::t

test_that("Exception are raised", {
  expect_error(trainValTest(t@data, .2), "x must be an SWD object!")
})

test_that("The output is correct for train/test split", {
  x <- trainValTest(t, test = 0.3, only_presence = FALSE, seed = 25)
  np <- length(which(t@pa == 1))
  na <- length(which(t@pa == 0))
  expect_type(x, "list")
  expect_length(x, 2)
  # Training SWD
  expect_s4_class(x[[1]], "SWD")
  expect_equal(nrow(x[[1]]@data[x[[1]]@pa == 1, ]), round(np * 0.7, 0))
  expect_equal(nrow(x[[1]]@coords[x[[1]]@pa == 1, ]), round(np * 0.7, 0))
  expect_equal(nrow(x[[1]]@data[x[[1]]@pa == 0, ]), round(na * 0.7, 0))
  expect_equal(nrow(x[[1]]@coords[x[[1]]@pa == 0, ]), round(na * 0.7, 0))
  expect_equal(nrow(x[[1]]@data), nrow(x[[1]]@coords), length(x[[1]]@pa))
  # Testing SWD
  expect_s4_class(x[[2]], "SWD")
  expect_equal(nrow(x[[2]]@data[x[[2]]@pa == 1, ]), round(np * 0.3, 0))
  expect_equal(nrow(x[[2]]@coords[x[[2]]@pa == 1, ]), round(np * 0.3, 0))
  expect_equal(nrow(x[[2]]@data[x[[2]]@pa == 0, ]), round(na * 0.3, 0))
  expect_equal(nrow(x[[2]]@coords[x[[2]]@pa == 0, ]), round(na * 0.3, 0))
  expect_equal(nrow(x[[2]]@data), nrow(x[[2]]@coords), length(x[[2]]@pa))

  expect_equal(nrow(rbind(x[[1]]@data, x[[2]]@data)), nrow(t@data))
  expect_equal(nrow(rbind(x[[1]]@coords, x[[2]]@coords)), nrow(t@coords))
  expect_equal(length(c(x[[1]]@pa, x[[2]]@pa)), length(t@pa))
})

test_that("The output is correct for train/test split for only presence", {
  # Train/test split, only presence
  x <- trainValTest(t, test = 0.3, only_presence = TRUE)
  np <- length(which(t@pa == 1))
  na <- length(which(t@pa == 0))
  expect_type(x, "list")
  expect_length(x, 2)
  expect_s4_class(x[[1]], "SWD")
  expect_equal(nrow(x[[1]]@data[x[[1]]@pa == 1, ]), round(np * 0.7, 0))
  expect_equal(nrow(x[[1]]@coords[x[[1]]@pa == 1, ]), round(np * 0.7, 0))
  expect_equal(nrow(x[[1]]@data[x[[1]]@pa == 0, ]), na)
  expect_equal(nrow(x[[1]]@coords[x[[1]]@pa == 0, ]), na)
  expect_equal(nrow(x[[1]]@data), nrow(x[[1]]@coords), length(x[[1]]@pa))
  expect_s4_class(x[[2]], "SWD")
  expect_equal(nrow(x[[2]]@data[x[[2]]@pa == 1, ]), round(np * 0.3, 0))
  expect_equal(nrow(x[[2]]@coords[x[[2]]@pa == 1, ]), round(np * 0.3, 0))
  expect_equal(nrow(x[[2]]@data[x[[2]]@pa == 0, ]), na)
  expect_equal(nrow(x[[2]]@coords[x[[2]]@pa == 0, ]), na)
  expect_equal(nrow(x[[2]]@data), nrow(x[[2]]@coords), length(x[[2]]@pa))

  expect_equal(nrow(rbind(x[[1]]@data, x[[2]]@data)),
               nrow(rbind(t@data, t@data[t@pa == 0, ])))
  expect_equal(nrow(rbind(x[[1]]@coords, x[[2]]@coords)),
               nrow(rbind(t@coords, t@coords[t@pa == 0, ])))
  expect_equal(length(c(x[[1]]@pa, x[[2]]@pa)),
               length(c(t@pa, t@pa[t@pa == 0])))
})

test_that("The output is correct for train/val/test split", {
  x <- trainValTest(t, val = 0.2, test = 0.2, only_presence = FALSE)
  np <- length(which(t@pa == 1))
  na <- length(which(t@pa == 0))
  expect_type(x, "list")
  expect_length(x, 3)
  expect_s4_class(x[[1]], "SWD")
  expect_equal(nrow(x[[1]]@data[x[[1]]@pa == 1, ]), round(np * 0.6, 0))
  expect_equal(nrow(x[[1]]@coords[x[[1]]@pa == 1, ]), round(np * 0.6, 0))
  expect_equal(nrow(x[[1]]@data[x[[1]]@pa == 0, ]), round(na * 0.6, 0))
  expect_equal(nrow(x[[1]]@coords[x[[1]]@pa == 0, ]), round(na * 0.6, 0))
  expect_equal(nrow(x[[1]]@data), nrow(x[[1]]@coords), length(x[[1]]@pa))
  expect_s4_class(x[[2]], "SWD")
  expect_equal(nrow(x[[2]]@data[x[[2]]@pa == 1, ]), round(np * 0.2, 0))
  expect_equal(nrow(x[[2]]@coords[x[[2]]@pa == 1, ]), round(np * 0.2, 0))
  expect_equal(nrow(x[[2]]@data[x[[2]]@pa == 0, ]), round(na * 0.2, 0))
  expect_equal(nrow(x[[2]]@coords[x[[2]]@pa == 0, ]), round(na * 0.2, 0))
  expect_equal(nrow(x[[2]]@data), nrow(x[[2]]@coords), length(x[[2]]@pa))
  expect_s4_class(x[[3]], "SWD")
  expect_equal(nrow(x[[3]]@data[x[[3]]@pa == 1, ]), round(np * 0.2, 0))
  expect_equal(nrow(x[[3]]@coords[x[[3]]@pa == 1, ]), round(np * 0.2, 0))
  expect_equal(nrow(x[[3]]@data[x[[3]]@pa == 0, ]), round(na * 0.2, 0))
  expect_equal(nrow(x[[3]]@coords[x[[3]]@pa == 0, ]), round(na * 0.2, 0))
  expect_equal(nrow(x[[3]]@data), nrow(x[[3]]@coords), length(x[[3]]@pa))

  expect_equal(nrow(rbind(x[[1]]@data, x[[2]]@data, x[[3]]@data)), nrow(t@data))
  expect_equal(nrow(rbind(x[[1]]@coords, x[[2]]@coords, x[[3]]@coords)),
               nrow(t@coords))
  expect_equal(length(c(x[[1]]@pa, x[[2]]@pa, x[[3]]@pa)), length(t@pa))
})

test_that("The output is correct for train/val/test split for only presence", {
  x <- trainValTest(t, val = 0.2, test = 0.2, only_presence = TRUE)
  np <- length(which(t@pa == 1))
  na <- length(which(t@pa == 0))
  expect_type(x, "list")
  expect_length(x, 3)
  expect_s4_class(x[[1]], "SWD")
  expect_equal(nrow(x[[1]]@data[x[[1]]@pa == 1, ]), round(np * 0.6, 0))
  expect_equal(nrow(x[[1]]@coords[x[[1]]@pa == 1, ]), round(np * 0.6, 0))
  expect_equal(nrow(x[[1]]@data[x[[1]]@pa == 0, ]), na)
  expect_equal(nrow(x[[1]]@coords[x[[1]]@pa == 0, ]), na)
  expect_equal(nrow(x[[1]]@data), nrow(x[[1]]@coords), length(x[[1]]@pa))
  expect_s4_class(x[[2]], "SWD")
  expect_equal(nrow(x[[2]]@data[x[[2]]@pa == 1, ]), round(np * 0.2, 0))
  expect_equal(nrow(x[[2]]@coords[x[[2]]@pa == 1, ]), round(np * 0.2, 0))
  expect_equal(nrow(x[[2]]@data[x[[2]]@pa == 0, ]), na)
  expect_equal(nrow(x[[2]]@coords[x[[2]]@pa == 0, ]), na)
  expect_equal(nrow(x[[2]]@data), nrow(x[[2]]@coords), length(x[[2]]@pa))
  expect_s4_class(x[[3]], "SWD")
  expect_equal(nrow(x[[3]]@data[x[[3]]@pa == 1, ]), round(np * 0.2, 0))
  expect_equal(nrow(x[[3]]@coords[x[[3]]@pa == 1, ]), round(np * 0.2, 0))
  expect_equal(nrow(x[[3]]@data[x[[3]]@pa == 0, ]), na)
  expect_equal(nrow(x[[3]]@coords[x[[3]]@pa == 0, ]), na)
  expect_equal(nrow(x[[3]]@data), nrow(x[[3]]@coords), length(x[[3]]@pa))

  expect_equal(nrow(rbind(x[[1]]@data, x[[2]]@data, x[[3]]@data)),
               nrow(rbind(t@data, t@data[t@pa == 0, ], t@data[t@pa == 0, ])))
  expect_equal(nrow(rbind(x[[1]]@coords, x[[2]]@coords, x[[3]]@coords)),
               nrow(rbind(t@coords, t@coords[t@pa == 0, ],
                          t@coords[t@pa == 0, ])))
  expect_equal(length(c(x[[1]]@pa, x[[2]]@pa, x[[3]]@pa)),
               length(c(t@pa, t@pa[t@pa == 0], t@pa[t@pa == 0])))
})
