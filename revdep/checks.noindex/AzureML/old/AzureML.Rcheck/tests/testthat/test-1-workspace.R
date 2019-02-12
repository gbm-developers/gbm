if(interactive()) library("testthat")


settingsFile <- AzureML.config.default
workspace <- function(..., .validate = FALSE) AzureML::workspace(..., .validate = .validate)

#  ------------------------------------------------------------------------

context("workspace - connect to workspace")

test_that("Can connect to workspace with supplied id and auth", {
  AzureML:::skip_if_missing_config(settingsFile)
  
  js <- read.AzureML.config(settingsFile)
  id <- js$workspace$id
  auth <- js$workspace$authorization_token
  
  expect_true(!is.null(id))
  expect_true(!is.null(auth))
  
  ws <- workspace(id, auth)
  
  expect_is(ws, c("Workspace"))
  expect_equal(ls(ws), c("datasets", "experiments", "id", "services"))
  expect_equal(ws$id, id)
})

test_that("Can connect to workspace with config file", {
  AzureML:::skip_if_missing_config(settingsFile)
  
  ws <- workspace()
  
  expect_is(ws, c("Workspace"))
  expect_equal(ls(ws), c("datasets", "experiments", "id", "services"))
})


test_that("workspace() throws helpful 401 error with invalid id", {
  # AzureML:::skip_if_missing_config(settingsFile)
  
  .catchError <- function(expr){
    tryCatch(expr, error = function(e)e)$message
  }
  .expect_error_in <- function(object, msgs){
    if(missing(object) || is.null(object)) return(FALSE)
    ptn <- sprintf("[%s]", paste(sprintf("(%s)", msgs), collapse = "|"))
    expect_true(grepl(ptn, object))
  }
  
  m <- .catchError(workspace(id = "x", auth = "y", .validate = TRUE))
  msg <- c("invalid workspaceId", 
           "401 (Unauthorised).  Please check your workspace ID and auth codes."
  )
  
  .expect_error_in(m, msg = msg)
  
})




#  ------------------------------------------------------------------------

context("workspace - reading from settings.json file")

test_that("workspace() adds api_endpoint and management_endpoint if missing from config", {
  tf <- tempfile(fileext = ".json")
  on.exit(unlink(tf))
  write.AzureML.config("x", "y", file = tf)
  ws <- workspace(config = tf)
  expect_equal(ws$id, "x")
  expect_equal(
    ws$.api_endpoint, 
    default_api(ws$.api_endpoint)[["api_endpoint"]]
  )
  expect_equal(
    ws$.management_endpoint, 
    default_api(ws$.api_endpoint)[["management_endpoint"]]
  )
})

test_that("workspace() throws helpful error if config file does not exist", {
  expect_error(
    workspace(config = "file_does_not_exist"),
    "config file is missing: 'file_does_not_exist'"
  )
})

test_that("workspace() throws helpful error if config is invalid json", {
  tf <- tempfile(fileext = ".json")
  on.exit(unlink(tf))
  writeLines("garbage", con = tf)
  msg <- tryCatch(workspace(config = tf), error = function(e)e)$message
  expect_true(
    grepl("Your config file contains invalid json", msg)
  )
})

