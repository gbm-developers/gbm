if(interactive()) library(testthat)

context("Discover API schema")


test_that("discoverSchema() returns help page information", {
  AzureML:::skip_if_offline()

  schemaUrl <- "https://studio.azureml.net/apihelp/workspaces/xxxxx/webservices/yyyyy/endpoints/zzzzz"
  expect_equal(
    getDetailsFromUrl(schemaUrl),
    c("xxxxx", 
      "yyyyy",
      "zzzzz")
  )
  
  schemaUrl <- "https://studio.azureml.net/apihelp/workspaces/f5e8e9bc4eed4034b78567449cfca779/webservices/d42667a354e34a3f98888ba86300fc2f/endpoints/b4caf0f0ebfd451bbc187741894e213b/score"
  
  expect_equal(
    getDetailsFromUrl(schemaUrl),
    c("f5e8e9bc4eed4034b78567449cfca779", 
      "d42667a354e34a3f98888ba86300fc2f",
      "b4caf0f0ebfd451bbc187741894e213b")
  )
  
  url <- "https://ussouthcentral.services.azureml.net/workspaces/f5e8e9bc4eed4034b78567449cfca779/services/b4caf0f0ebfd451bbc187741894e213b/execute?api-version=2.0&format=swagger"
  expect_error(
    getDetailsFromUrl(url)
  )
  
  schema <- discoverSchema(schemaUrl)
  schema$sampleInput$Gender <- "male"
  schema$sampleInput$PortEmbarkation <- "C"
  
  expect_equal(length(schema), 4)
  expect_equivalent(schema$requestUrl, url)
  expect_equivalent(schema$columnNames, 
                    list("Survived", 
                         "PassengerClass", 
                         "Gender", 
                         "Age", 
                         "SiblingSpouse", 
                         "ParentChild", 
                         "FarePrice", 
                         "PortEmbarkation")
  )
  expect_equivalent(schema$sampleInput, 
                    list(Survived = 1, 
                         PassengerClass = 1, 
                         Gender = "male", 
                         Age = 1, 
                         SiblingSpouse = 1, 
                         ParentChild = 1, 
                         FarePrice = 1, 
                         PortEmbarkation = "C"))
})

