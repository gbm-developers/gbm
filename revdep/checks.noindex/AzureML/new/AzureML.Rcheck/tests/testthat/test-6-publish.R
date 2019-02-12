if(interactive()) library(testthat)

context("Publish API")
settingsFile <- AzureML.config.default


test_that(".getexports finds function and creates zip string", {
  AzureML:::skip_if_missing_config(settingsFile)
  AzureML:::skip_if_offline()
  
  ws <<- workspace()
  endpoint <<- NA
  
  
  funEnv <- new.env()
  assign("add", function(x, y) x + y, envir = funEnv)
  
  exportEnv = new.env()
  AzureML:::.getexports(substitute(add), e = exportEnv, env = funEnv)
  
  expect_equal(
    ls(exportEnv),
    "add"
  )
  
  za <- AzureML:::zipAvailable()
  if(!za) skip(AzureML:::zipNotAvailableMessage)
  expect_true(za)
  
  z <- AzureML:::packageEnv(exportEnv)
  expect_is(z, "character")
  expect_true(nchar(z) > 1)
  
})



test_that("publishWebService throws error if fun is not a function", {
  AzureML:::skip_if_missing_config(settingsFile)
  AzureML:::skip_if_offline()

  add <- function(x,y) x + y
  
  timestamped_name <- paste0("webservice-test-publish-",
                             format(Sys.time(), format="%Y-%m-%d--%H-%M-%S"))
  
  expect_error({
    endpoint <- publishWebService(ws,
                                  fun = "add",
                                  name = timestamped_name,
                                  inputSchema = list(x="numeric",
                                                     y="numeric"),
                                  outputSchema = list(ans="numeric")
    )
    if(is.Endpoint(endpoint)) deleteWebService(ws, timestamped_name)
  },
  "You must specify 'fun' as a function, not a character"
  )
})

timestamped_name <- paste0("webservice-test-publish-",
                           format(Sys.time(), format="%Y-%m-%d--%H-%M-%S"))



test_that("publishWebService works with simple function", {
  AzureML:::skip_if_missing_config(settingsFile)
  AzureML:::skip_if_offline()
  
  add <- function(x,y) x + y
  
  endpoint <- publishWebService(ws,
                                fun = add,
                                name = timestamped_name,
                                inputSchema = list(x="numeric",
                                                   y="numeric"),
                                outputSchema = list(ans="numeric")
  )
  
  endpoint <<- endpoint # Used to test updateWebservice in next test
  
  
  expect_is(endpoint, "data.frame")
  expect_is(endpoint, "Endpoint")
  expect_is(endpoint$WorkspaceId, "character")
  expect_is(endpoint$WebServiceId, "character")
  expect_equal(ws$id, endpoint$WorkspaceId)
  
  # Now test if we can consume the service we just published
  res <- consume(endpoint, x=pi, y=2)
  expect_is(res, "data.frame")
  expect_equal(res$ans, pi + 2, tolerance = 1e-8)
})


test_that("updateWebService works with simple function", {
  # Now test updateWebService
  AzureML:::skip_if_missing_config(settingsFile)
  AzureML:::skip_if_offline()
  
  endpoint <- updateWebService(ws,
                               serviceId = endpoint$WebServiceId,
                               fun = function(x, y) x - y,
                               inputSchema = list(x="numeric",
                                                  y="numeric"),
                               outputSchema = list(ans="numeric"))
  
  # Now test if we can consume the service we just updated
  for(i in 1:10){
    Sys.sleep(3) # Allow some time for the service to update and refresh
    res <- consume(endpoint, x=pi, y=2)
    if(isTRUE(all.equal(res$ans, pi - 2, tolerance = 1e-8))) break
  }
  expect_is(res, "data.frame")
  expect_equal(res$ans, pi - 2, tolerance = 1e-8)
  
  deleteWebService(ws, timestamped_name)
})


test_that("publishWebService works with data frame input", {
  AzureML:::skip_if_missing_config(settingsFile)
  AzureML:::skip_if_offline()
  
  timestamped_name <- paste0("webservice-test-publish-",
                             format(Sys.time(), format="%Y-%m-%d--%H-%M-%S"))
  
  if(!require("lme4")) skip("You need to install lme4 to run this test")
  
  set.seed(1)
  train <- sleepstudy[sample(nrow(sleepstudy), 120),]
  m <- lm(Reaction ~ Days + Subject, data = train)
  
  # Deine a prediction function to publish based on the model:
  sleepyPredict <- function(newdata){
    predict(m, newdata=newdata)
  }
  
  endpoint <- publishWebService(ws, fun = sleepyPredict, 
                                name = timestamped_name,
                                inputSchema = sleepstudy)
  
  expect_is(endpoint, "data.frame")
  expect_is(endpoint, "Endpoint")
  expect_is(endpoint$WorkspaceId, "character")
  expect_is(endpoint$WebServiceId, "character")
  expect_equal(ws$id, endpoint$WorkspaceId)
  
  
  # Now test if we can consume the service we just published
  res <- consume(endpoint, sleepstudy)$ans
  expect_is(res, "numeric")
  expect_equal(length(res), nrow(sleepstudy))
  
  deleteWebService(ws, timestamped_name)
})

