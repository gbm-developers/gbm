if(!file.exists("~/.azureml/settings.json")){
  message("To run tests, add a file ~/.azureml/settings.json containing AzureML keys")
  message("Some tests skipped. See ?workspace for help")
}