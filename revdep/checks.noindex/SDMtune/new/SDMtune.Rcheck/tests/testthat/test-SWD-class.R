test_that("Show method for SWD class produces the correct output", {
  data <- SDMtune:::t
  expect_output(print(data), "Object of class SWD", fixed = TRUE)
  expect_output(print(data), "Species: Virtual species", fixed = TRUE)
  expect_output(print(data), "Presence locations: 400", fixed = TRUE)
  expect_output(print(data), "Absence locations: 5000", fixed = TRUE)
  expect_output(print(data), "bio1 bio12 bio16 bio17 bio5 bio6 bio7 bio8",
                fixed = TRUE)
  expect_output(print(data), "biome")
  data@data <- data@data[, 1:3]
  expect_output(print(data), "Continuous: bio1 bio12 bio16", fixed = TRUE)
  expect_output(print(data), "Categorical: NA", fixed = TRUE)
  data@data <- SDMtune:::t@data[, "biome", drop = FALSE]
  expect_output(print(data), "Continuous: NA", fixed = TRUE)
  expect_output(print(data), "Categorical: biome", fixed = TRUE)
})
