if(interactive()) library("testthat")

settingsFile <- AzureML.config.default

context("Download multiple datasets")

test_that("datasets(ws) returns results", {
  AzureML:::skip_if_missing_config(settingsFile)
  AzureML:::skip_if_offline()

  ws <<- workspace()
  
  x <- datasets(ws)
  expect_is(x, "data.frame")
})

### Additional tests of download.datasets(.).  We could expand the same dataset formats
### to other tests (upload, delete, etc).

## csv and .tsv files:
test_that("Can download multiple .csv and .tsv files", {
  AzureML:::skip_if_missing_config(settingsFile)
  AzureML:::skip_if_offline()
  
  # ds <- datasets(ws, filter = "samples")
  # ds[grep("[CT]SV", ds$DataTypeId), ]
  names <- c("Time Series Dataset", 
             "Sample Named Entity Recognition Articles")
  
  res <- suppressWarnings(download.datasets(ws, names))
  expect_equal(names(res), names)
  
  res <- suppressWarnings(download.datasets(datasets(ws), names))
  expect_equal(names(res), names)
  
})

test_that("Can download .zip files", {
  AzureML:::skip_if_missing_config(settingsFile)
  AzureML:::skip_if_offline()
  
  # ds <- datasets(ws, filter = "samples")
  # ds[ds$DataTypeId == "Zip", ]
  names <- c("text.preprocessing.zip", "fraudTemplateUtil.zip")
  
  res <- download.datasets(ws, names)
  expect_equal(names(res), names)
  
  res <- download.datasets(datasets(ws), names)
  expect_equal(names(res), names)
})


test_that("Can download .arff files", {
  AzureML:::skip_if_missing_config(settingsFile)
  AzureML:::skip_if_offline()
  
  # ds <- datasets(ws, filter = "samples")
  # ds[ds$DataTypeId == "ARFF", ]
  names <- c("Breast cancer data", "Forest fires data", "Iris Two Class Data")
  
  res <- download.datasets(ws, names)
  expect_equal(names(res), names)
  
  res <- download.datasets(datasets(ws), names)
  expect_equal(names(res), names)
})

