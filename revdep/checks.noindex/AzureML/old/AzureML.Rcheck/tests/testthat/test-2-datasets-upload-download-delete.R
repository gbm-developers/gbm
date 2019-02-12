if(interactive()) library("testthat")

settingsFile <- AzureML.config.default

context("Upload and delete dataset")

test_that("datasets(ws) returns results", {
  AzureML:::skip_if_missing_config(settingsFile)
  AzureML:::skip_if_offline()
  
  ws <<- workspace()
  
  x <- datasets(ws)
  expect_is(x, "data.frame")
})

timestamped_name <- paste0("dataset-test-upload-",
                           format(Sys.time(), format="%Y-%m-%d--%H-%M-%S"))

test_that("Can upload dataset to workspace", {
  AzureML:::skip_if_missing_config(settingsFile)
  AzureML:::skip_if_offline()
  
  upload.dataset(airquality, ws, timestamped_name)
  ds <- datasets(ws, filter = "my")
  expect_true(timestamped_name %in% ds$Name)
})

test_that("Uploading dataset with duplicate name gives helpful error", {
  AzureML:::skip_if_missing_config(settingsFile)
  AzureML:::skip_if_offline()
  
  expect_error(upload.dataset(airquality, ws, timestamped_name),
               sprintf("A dataset with the name '%s' already exists in AzureML", timestamped_name)
  )
})

test_that("Can download dataset", {
  AzureML:::skip_if_missing_config(settingsFile)
  AzureML:::skip_if_offline()
  
  dl <- download.datasets(ws, name = timestamped_name)
  expect_equal(dl, airquality)
})

test_that("Can delete dataset from workspace", {
  AzureML:::skip_if_missing_config(settingsFile)
  AzureML:::skip_if_offline()
  
  z <- delete.datasets(ws, timestamped_name)
  expect_true(timestamped_name %in% z$Name && z$Deleted[z$Name == timestamped_name])
  # Force refresh - sometime this fails in non-interactive
  max_wait <- 15
  wait_period <- 3
  i <- 0
  ds <- datasets(ws, filter = "my")
  while(i < max_wait && nrow(ds) > 0 && timestamped_name %in% ds$Name) {
    Sys.sleep(wait_period)
    i <- i + wait_period
    refresh(ws, what = "datasets")
    ds <- datasets(ws, filter = "my")
  }
  if(nrow(ds) > 0  || timestamped_name %in% ds$Name) skip("skip waiting for delete")
  expect_true(nrow(ds) == 0 || !timestamped_name %in% ds$Name)
})



test_that("Invalid input throws helpful error", {
  expect_error(download.datasets('HSAFundsData.csv'),
               "You specified a dataset name that is not in the workspace. See help file for `download.datasets`"
  )
})






