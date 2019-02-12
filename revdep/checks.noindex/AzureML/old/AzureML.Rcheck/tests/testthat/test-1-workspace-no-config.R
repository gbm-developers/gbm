if(interactive()) library("testthat")


settingsFile <- AzureML.config.default
workspace <- function(..., .validate = FALSE) AzureML::workspace(..., .validate = .validate)

#  ------------------------------------------------------------------------

context("workspace - connect to workspace with no config file")

test_that("Can connect to workspace with no config file", {

  opts <- getOption("AzureML.config")
  options(AzureML.config = tempfile(fileext = ".tmp"))
  on.exit(options(AzureML.config = opts))
  
  expect_error(
    ws <- workspace(), 
    "config file is missing"
  )
  expect_is(workspace("x", "y"), "Workspace")
  expect_equal({ws <- workspace("x", "y"); ls(ws)}, 
               c("datasets", "experiments", "id", "services"))
})
