Sys.setenv("R_TESTS" = "")
library(testthat)
library(personalized)

test_check("personalized")
