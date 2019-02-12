Sys.setenv("R_TESTS" = "")
library(testthat, quietly = TRUE)
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  # NOT_CRAN
  # run all tests
  test_check("AzureML")
  
} else {
  # CRAN
  # skip some tests on CRAN, to comply with timing directive and other policy
  test_check("AzureML", filter = "1-workspace-no-config")
  test_check("AzureML", filter = "7-discover-schema")
}
