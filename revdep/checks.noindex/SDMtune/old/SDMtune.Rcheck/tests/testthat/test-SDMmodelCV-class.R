test_that("Show method for SDMmodelCV class produces the correct output", {
  m <- SDMtune:::bm_maxent_cv
  expect_output(print(m), "Object of class SDMmodelCV", fixed = TRUE)
  expect_output(print(m), "Method: Maxent", fixed = TRUE)
  expect_output(print(m), "Replicates: 4", fixed = TRUE)
  expect_output(print(m), "Species: Virtual species", fixed = TRUE)
  expect_output(print(m), "Presence locations: 400", fixed = TRUE)
  expect_output(print(m), "Absence locations: 5000", fixed = TRUE)
  expect_output(print(m), "fc: lqph", fixed = TRUE)
  expect_output(print(m), "reg: 1", fixed = TRUE)
  expect_output(print(m), "iter: 500", fixed = TRUE)
  expect_output(print(m), "bio1 bio12 bio16 bio17 bio5 bio6 bio7 bio8",
                fixed = TRUE)
  expect_output(print(m), "biome")
  m@data@data <- SDMtune:::t@data[, 1:3]
  expect_output(print(m), "Continuous: bio1 bio12 bio16", fixed = TRUE)
  expect_output(print(m), "Categorical: NA", fixed = TRUE)
  m@data@data <- SDMtune:::t@data[, "biome", drop = FALSE]
  expect_output(print(m), "Continuous: NA", fixed = TRUE)
  expect_output(print(m), "Categorical: biome", fixed = TRUE)
  m <- SDMtune:::bm_maxnet_cv
  expect_output(print(m), "Object of class SDMmodelCV", fixed = TRUE)
  expect_output(print(m), "Method: Maxnet", fixed = TRUE)
  expect_output(print(m), "Replicates: 4", fixed = TRUE)
  expect_output(print(m), "Species: Virtual species", fixed = TRUE)
  expect_output(print(m), "Presence locations: 400", fixed = TRUE)
  expect_output(print(m), "Absence locations: 5000", fixed = TRUE)
  expect_output(print(m), "fc: lqph", fixed = TRUE)
  expect_output(print(m), "reg: 1", fixed = TRUE)
  expect_output(print(m), "bio1 bio12 bio16 bio17 bio5 bio6 bio7 bio8",
                fixed = TRUE)
  expect_output(print(m), "biome")
})
