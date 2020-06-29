# Temp files for xgboost
test_pmml_file <- tempfile()

teardown(unlink(c(test_pmml_file), recursive = TRUE))


test_that("function returns a parsed XML node", {
  iris_box <- xform_wrap(iris)
  iris_box <- xform_z_score(iris_box, xform_info = "column1->d1")
  iris_box <- xform_z_score(iris_box, xform_info = "column2->d2")

  # Make a LocalTransformations element and save it to an external file:
  pmml_trans <- pmml(NULL, transforms = iris_box)
  write(toString(pmml_trans), file = test_pmml_file)
  # write(pmml_trans, file = test_pmml_file)

  # Later, we may need to read in the PMML model into R
  # 'lt' below is now a XML Node, as opposed to a string:
  lt <- file_to_xml_node(test_pmml_file)
  expect_equal(toString(lt), "<LocalTransformations>\n <DerivedField name=\"d1\" dataType=\"double\" optype=\"continuous\">\n  <NormContinuous field=\"Sepal.Length\">\n   <LinearNorm orig=\"5.84333333333333\" norm=\"0\"/>\n   <LinearNorm orig=\"6.6713994613112\" norm=\"1\"/>\n  </NormContinuous>\n </DerivedField>\n <DerivedField name=\"d2\" dataType=\"double\" optype=\"continuous\">\n  <NormContinuous field=\"Sepal.Width\">\n   <LinearNorm orig=\"3.05733333333333\" norm=\"0\"/>\n   <LinearNorm orig=\"3.49319961827003\" norm=\"1\"/>\n  </NormContinuous>\n </DerivedField>\n</LocalTransformations>")
})
