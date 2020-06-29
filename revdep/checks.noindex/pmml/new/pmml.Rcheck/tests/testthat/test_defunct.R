dummy_file <- tempfile()

teardown(unlink(c(dummy_file), recursive = TRUE))

test_that("Error message is displayed when using defunct function name", {
  data(iris)

  iris_box <- xform_wrap(iris)

  expect_error(
    DiscretizeXform(iris_box,
      xform_info = "[Sepal.Length -> dsl][double -> string]",
      table = dummy_file, map_missing_to = "0"
    ),
    "This function is defunct; use xform_discretize instead."
  )

  expect_error(
    FunctionXform(iris_box,
      orig_field_name = "Sepal.Length",
      new_field_name = "Sepal.Length.Transformed",
      expression = "(Sepal.Length^2)/100"
    ),
    "This function is defunct; use xform_function instead."
  )

  expect_error(
    MapXform(audit_box,
      xform_info = "[Sex -> d_sex][string->integer]",
      table = dummy_file, map_missing_to = "0"
    ),
    "This function is defunct; use xform_map instead."
  )

  expect_error(
    MinMaxXform(iris_box),
    "This function is defunct; use xform_min_max instead."
  )

  expect_error(
    NormDiscreteXform(iris_box, xform_info = "Species"),
    "This function is defunct; use xform_norm_discrete instead."
  )

  expect_error(
    RenameVar(wrap_object = iris_box, xform_info = "Sepal.Width->SW"),
    "This function is defunct; use rename_wrap_var instead."
  )

  expect_error(WrapData(iris), "This function is defunct; use xform_wrap instead.")

  expect_error(ZScoreXform(iris_box), "This function is defunct; use xform_z_score instead.")
})
