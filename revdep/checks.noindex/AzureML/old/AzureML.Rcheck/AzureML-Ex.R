pkgname <- "AzureML"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('AzureML')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("consume")
### * consume

flush(stderr()); flush(stdout())

### Name: consume
### Title: Use a web service to score data in list (key=value) format.
### Aliases: consume

### ** Examples

## Not run: 
##D # Use a default configuration in ~/.azureml, alternatively
##D # see help for `?workspace`.
##D 
##D ws <- workspace()
##D   
##D # Publish a simple model using the lme4::sleepdata ---------------------------
##D 
##D library(lme4)
##D set.seed(1)
##D train <- sleepstudy[sample(nrow(sleepstudy), 120),]
##D m <- lm(Reaction ~ Days + Subject, data = train)
##D 
##D # Deine a prediction function to publish based on the model:
##D sleepyPredict <- function(newdata){
##D   predict(m, newdata=newdata)
##D }
##D 
##D ep <- publishWebService(ws, fun = sleepyPredict, name="sleepy lm",
##D                         inputSchema = sleepstudy,
##D                         data.frame=TRUE)
##D 
##D # OK, try this out, and compare with raw data
##D ans <- consume(ep, sleepstudy)$ans
##D plot(ans, sleepstudy$Reaction)
##D 
##D # Remove the service
##D deleteWebService(ws, "sleepy lm")
##D 
##D 
##D 
##D # Another data frame example -------------------------------------------------
##D 
##D # If your function can consume a whole data frame at once, you can also
##D # supply data in that form, resulting in more efficient computation.
##D # The following example builds a simple linear model on a subset of the
##D # airquality data and publishes a prediction function based on the model.
##D set.seed(1)
##D m <- lm(Ozone ~ ., data=airquality[sample(nrow(airquality), 100),])
##D # Define a prediction function based on the model:
##D fun <- function(newdata)
##D {
##D   predict(m, newdata=newdata)
##D }
##D # Note the definition of inputSchema and use of the data.frame argument.
##D ep <- publishWebService(ws, fun=fun, name="Ozone",
##D                         inputSchema = airquality,
##D                         data.frame=TRUE)
##D ans <- consume(ep, airquality)$ans
##D plot(ans, airquality$Ozone)
##D deleteWebService(ws, "Ozone")
##D 
##D 
##D 
##D # Train a model using diamonds in ggplot2 ------------------------------------
##D # This example also demonstrates how to deal with factor in the data
##D 
##D data(diamonds, package="ggplot2")
##D set.seed(1)
##D train_idx = sample.int(nrow(diamonds), 30000)
##D test_idx = sample(setdiff(seq(1, nrow(diamonds)), train_idx), 500)
##D train <- diamonds[train_idx, ]
##D test  <- diamonds[test_idx, ]
##D 
##D model <- glm(price ~ carat + clarity + color + cut - 1, data = train, 
##D              family = Gamma(link = "log"))
##D 
##D diamondLevels <- diamonds[1, ]
##D 
##D # The model works reasonably well, except for some outliers
##D plot(exp(predict(model, test)) ~ test$price)
##D 
##D # Create a prediction function that converts characters correctly to factors
##D 
##D predictDiamonds <- function(x){
##D   x$cut     <- factor(x$cut,     
##D                       levels = levels(diamondLevels$cut), ordered = TRUE)
##D   x$clarity <- factor(x$clarity, 
##D                       levels = levels(diamondLevels$clarity), ordered = TRUE)
##D   x$color   <- factor(x$color,   
##D                       levels = levels(diamondLevels$color), ordered = TRUE)
##D   exp(predict(model, newdata = x))
##D }
##D 
##D 
##D # Publish the service
##D 
##D ws <- workspace()
##D ep <- publishWebService(ws, fun = predictDiamonds, name = "diamonds",
##D                         inputSchema = test,
##D                         data.frame = TRUE
##D )
##D 
##D # Consume the service
##D results <- consume(ep, test)$ans
##D plot(results ~ test$price)
##D 
##D deleteWebService(ws, "diamonds")
##D 
##D 
##D 
##D # Simple example using scalar input ------------------------------------------
##D 
##D ws <- workspace()
##D 
##D # Really simple example:
##D add <- function(x,y) x + y
##D endpoint <- publishWebService(ws, 
##D                               fun = add, 
##D                               name = "addme", 
##D                               inputSchema = list(x="numeric", 
##D                                                  y="numeric"), 
##D                               outputSchema = list(ans="numeric"))
##D consume(endpoint, list(x=pi, y=2))
##D 
##D # Now remove the web service named "addme" that we just published
##D deleteWebService(ws, "addme")
##D 
##D 
##D 
##D # Send a custom R function for evaluation in AzureML -------------------------
##D 
##D # A neat trick to evaluate any expression in the Azure ML virtual
##D # machine R session and view its output:
##D ep <- publishWebService(ws, 
##D                         fun =  function(expr) {
##D                           paste(capture.output(
##D                             eval(parse(text=expr))), collapse="\n")
##D                         },
##D                         name="commander", 
##D                         inputSchema = list(x = "character"),
##D                         outputSchema = list(ans = "character"))
##D cat(consume(ep, list(x = "getwd()"))$ans)
##D cat(consume(ep, list(x = ".packages(all=TRUE)"))$ans)
##D cat(consume(ep, list(x = "R.Version()"))$ans)
##D 
##D # Remove the service we just published
##D deleteWebService(ws, "commander")
##D 
##D 
##D 
##D # Understanding the scoping rules --------------------------------------------
##D 
##D # The following example illustrates scoping rules. Note that the function
##D # refers to the variable y defined outside the function body. That value
##D # will be exported with the service.
##D y <- pi
##D ep <- publishWebService(ws, 
##D                         fun = function(x) x + y, 
##D                         name = "lexical scope",
##D                         inputSchema = list(x = "numeric"), 
##D                         outputSchema = list(ans = "numeric"))
##D cat(consume(ep, list(x=2))$ans)
##D 
##D # Remove the service we just published
##D deleteWebService(ws, "lexical scope")
##D 
##D 
##D # Demonstrate scalar inputs but sending a data frame for scoring -------------
##D 
##D # Example showing the use of consume to score all the rows of a data frame
##D # at once, and other invocations for evaluating multiple sets of input
##D # values. The columns of the data frame correspond to the input parameters
##D # of the web service in this example:
##D f <- function(a,b,c,d) list(sum = a+b+c+d, prod = a*b*c*d)
##D ep <-  publishWebService(ws, 
##D                          f, 
##D                          name = "rowSums",
##D                          inputSchema = list(
##D                            a = "numeric", 
##D                            b = "numeric", 
##D                            c = "numeric", 
##D                            d = "numeric"
##D                          ),
##D                          outputSchema = list(
##D                            sum ="numeric", 
##D                            prod = "numeric")
##D )
##D x <- head(iris[,1:4])  # First four columns of iris
##D 
##D # Note the following will FAIL because of a name mismatch in the arguments
##D # (with an informative error):
##D consume(ep, x, retryDelay=1)
##D # We need the columns of the data frame to match the inputSchema:
##D names(x) <- letters[1:4]
##D # Now we can evaluate all the rows of the data frame in one call:
##D consume(ep, x)
##D # output should look like:
##D #    sum    prod
##D # 1 10.2   4.998
##D # 2  9.5   4.116
##D # 3  9.4  3.9104
##D # 4  9.4   4.278
##D # 5 10.2    5.04
##D # 6 11.4 14.3208
##D 
##D # You can use consume to evaluate just a single set of input values with this
##D # form:
##D consume(ep, a=1, b=2, c=3, d=4)
##D 
##D # or, equivalently,
##D consume(ep, list(a=1, b=2, c=3, d=4))
##D 
##D # You can evaluate multiple sets of input values with a data frame input:
##D consume(ep, data.frame(a=1:2, b=3:4, c=5:6, d=7:8))
##D 
##D # or, equivalently, with multiple lists:
##D consume(ep, list(a=1, b=3, c=5, d=7), list(a=2, b=4, c=6, d=8))
##D 
##D # Remove the service we just published
##D deleteWebService(ws, "rowSums")
##D 
##D # A more efficient way to do the same thing using data frame input/output:
##D f <- function(df) with(df, list(sum = a+b+c+d, prod = a*b*c*d))
##D ep = publishWebService(ws, f, name="rowSums2", 
##D                        inputSchema = data.frame(a = 0, b = 0, c = 0, d = 0))
##D consume(ep, data.frame(a=1:2, b=3:4, c=5:6, d=7:8))
##D deleteWebService(ws, "rowSums2")
##D 
##D 
##D 
##D # Automatically discover dependencies ----------------------------------------
##D 
##D # The publishWebService function uses `miniCRAN` to include dependencies on
##D # packages required by your function. The next example uses the `lmer`
##D # function from the lme4 package, and also shows how to publish a function
##D # that consumes a data frame by setting data.frame=TRUE.  Note! This example
##D # depends on a lot of packages and may take some time to upload to Azure.
##D library(lme4)
##D # Build a sample mixed effects model on just a subset of the sleepstudy data...
##D set.seed(1)
##D m <- lmer(Reaction ~ Days + (Days | Subject), 
##D           data=sleepstudy[sample(nrow(sleepstudy), 120),])
##D # Deine a prediction function to publish based on the model:
##D fun <- function(newdata)
##D {
##D   predict(m, newdata=newdata)
##D }
##D ep <- publishWebService(ws, fun=fun, name="sleepy lmer",
##D                         inputSchema= sleepstudy,
##D                         packages="lme4",
##D                         data.frame=TRUE)
##D 
##D # OK, try this out, and compare with raw data
##D ans = consume(ep, sleepstudy)$ans
##D plot(ans, sleepstudy$Reaction)
##D 
##D # Remove the service
##D deleteWebService(ws, "sleepy lmer")
## End(Not run)



cleanEx()
nameEx("datasets")
### * datasets

flush(stderr()); flush(stdout())

### Name: datasets
### Title: List datasets in an AzureML workspace.
### Aliases: datasets

### ** Examples

## Not run: 
##D   library(AzureML)
##D   
##D   # Use the default config file ~/azureml/settings.json with format:
##D   #   {"workspace":{
##D   #     "id":"test_id",
##D   #     "authorization_token": "test_token",
##D   #     "api_endpoint":"api_endpoint",
##D   #     "management_endpoint":"management_endpoint"
##D   #    }}
##D   # or, optionally set the `id` and `auth` parameters in the workspace
##D   # function.
##D   ws <- workspace()
##D   
##D   # List datasets
##D   ws$datasets
##D   datasets(ws)
##D   
##D   dataset <- "New York weather"
##D   ds <- match(dataset, ws$datasets$Name)
##D   frame <- download.datasets(ws$datasets[ds, ])
##D   head(frame)
##D 
##D   # Alternative approach:
##D   frame <- download.datasets(ws, name=dataset)
##D   head(frame)
## End(Not run)



cleanEx()
nameEx("deleteWebService")
### * deleteWebService

flush(stderr()); flush(stdout())

### Name: deleteWebService
### Title: Delete a Microsoft Azure Web Service
### Aliases: deleteWebService

### ** Examples

## Not run: 
##D # Use a default configuration in ~/.azureml, alternatively
##D # see help for `?workspace`.
##D 
##D ws <- workspace()
##D   
##D # Publish a simple model using the lme4::sleepdata ---------------------------
##D 
##D library(lme4)
##D set.seed(1)
##D train <- sleepstudy[sample(nrow(sleepstudy), 120),]
##D m <- lm(Reaction ~ Days + Subject, data = train)
##D 
##D # Deine a prediction function to publish based on the model:
##D sleepyPredict <- function(newdata){
##D   predict(m, newdata=newdata)
##D }
##D 
##D ep <- publishWebService(ws, fun = sleepyPredict, name="sleepy lm",
##D                         inputSchema = sleepstudy,
##D                         data.frame=TRUE)
##D 
##D # OK, try this out, and compare with raw data
##D ans <- consume(ep, sleepstudy)$ans
##D plot(ans, sleepstudy$Reaction)
##D 
##D # Remove the service
##D deleteWebService(ws, "sleepy lm")
##D 
##D 
##D 
##D # Another data frame example -------------------------------------------------
##D 
##D # If your function can consume a whole data frame at once, you can also
##D # supply data in that form, resulting in more efficient computation.
##D # The following example builds a simple linear model on a subset of the
##D # airquality data and publishes a prediction function based on the model.
##D set.seed(1)
##D m <- lm(Ozone ~ ., data=airquality[sample(nrow(airquality), 100),])
##D # Define a prediction function based on the model:
##D fun <- function(newdata)
##D {
##D   predict(m, newdata=newdata)
##D }
##D # Note the definition of inputSchema and use of the data.frame argument.
##D ep <- publishWebService(ws, fun=fun, name="Ozone",
##D                         inputSchema = airquality,
##D                         data.frame=TRUE)
##D ans <- consume(ep, airquality)$ans
##D plot(ans, airquality$Ozone)
##D deleteWebService(ws, "Ozone")
##D 
##D 
##D 
##D # Train a model using diamonds in ggplot2 ------------------------------------
##D # This example also demonstrates how to deal with factor in the data
##D 
##D data(diamonds, package="ggplot2")
##D set.seed(1)
##D train_idx = sample.int(nrow(diamonds), 30000)
##D test_idx = sample(setdiff(seq(1, nrow(diamonds)), train_idx), 500)
##D train <- diamonds[train_idx, ]
##D test  <- diamonds[test_idx, ]
##D 
##D model <- glm(price ~ carat + clarity + color + cut - 1, data = train, 
##D              family = Gamma(link = "log"))
##D 
##D diamondLevels <- diamonds[1, ]
##D 
##D # The model works reasonably well, except for some outliers
##D plot(exp(predict(model, test)) ~ test$price)
##D 
##D # Create a prediction function that converts characters correctly to factors
##D 
##D predictDiamonds <- function(x){
##D   x$cut     <- factor(x$cut,     
##D                       levels = levels(diamondLevels$cut), ordered = TRUE)
##D   x$clarity <- factor(x$clarity, 
##D                       levels = levels(diamondLevels$clarity), ordered = TRUE)
##D   x$color   <- factor(x$color,   
##D                       levels = levels(diamondLevels$color), ordered = TRUE)
##D   exp(predict(model, newdata = x))
##D }
##D 
##D 
##D # Publish the service
##D 
##D ws <- workspace()
##D ep <- publishWebService(ws, fun = predictDiamonds, name = "diamonds",
##D                         inputSchema = test,
##D                         data.frame = TRUE
##D )
##D 
##D # Consume the service
##D results <- consume(ep, test)$ans
##D plot(results ~ test$price)
##D 
##D deleteWebService(ws, "diamonds")
##D 
##D 
##D 
##D # Simple example using scalar input ------------------------------------------
##D 
##D ws <- workspace()
##D 
##D # Really simple example:
##D add <- function(x,y) x + y
##D endpoint <- publishWebService(ws, 
##D                               fun = add, 
##D                               name = "addme", 
##D                               inputSchema = list(x="numeric", 
##D                                                  y="numeric"), 
##D                               outputSchema = list(ans="numeric"))
##D consume(endpoint, list(x=pi, y=2))
##D 
##D # Now remove the web service named "addme" that we just published
##D deleteWebService(ws, "addme")
##D 
##D 
##D 
##D # Send a custom R function for evaluation in AzureML -------------------------
##D 
##D # A neat trick to evaluate any expression in the Azure ML virtual
##D # machine R session and view its output:
##D ep <- publishWebService(ws, 
##D                         fun =  function(expr) {
##D                           paste(capture.output(
##D                             eval(parse(text=expr))), collapse="\n")
##D                         },
##D                         name="commander", 
##D                         inputSchema = list(x = "character"),
##D                         outputSchema = list(ans = "character"))
##D cat(consume(ep, list(x = "getwd()"))$ans)
##D cat(consume(ep, list(x = ".packages(all=TRUE)"))$ans)
##D cat(consume(ep, list(x = "R.Version()"))$ans)
##D 
##D # Remove the service we just published
##D deleteWebService(ws, "commander")
##D 
##D 
##D 
##D # Understanding the scoping rules --------------------------------------------
##D 
##D # The following example illustrates scoping rules. Note that the function
##D # refers to the variable y defined outside the function body. That value
##D # will be exported with the service.
##D y <- pi
##D ep <- publishWebService(ws, 
##D                         fun = function(x) x + y, 
##D                         name = "lexical scope",
##D                         inputSchema = list(x = "numeric"), 
##D                         outputSchema = list(ans = "numeric"))
##D cat(consume(ep, list(x=2))$ans)
##D 
##D # Remove the service we just published
##D deleteWebService(ws, "lexical scope")
##D 
##D 
##D # Demonstrate scalar inputs but sending a data frame for scoring -------------
##D 
##D # Example showing the use of consume to score all the rows of a data frame
##D # at once, and other invocations for evaluating multiple sets of input
##D # values. The columns of the data frame correspond to the input parameters
##D # of the web service in this example:
##D f <- function(a,b,c,d) list(sum = a+b+c+d, prod = a*b*c*d)
##D ep <-  publishWebService(ws, 
##D                          f, 
##D                          name = "rowSums",
##D                          inputSchema = list(
##D                            a = "numeric", 
##D                            b = "numeric", 
##D                            c = "numeric", 
##D                            d = "numeric"
##D                          ),
##D                          outputSchema = list(
##D                            sum ="numeric", 
##D                            prod = "numeric")
##D )
##D x <- head(iris[,1:4])  # First four columns of iris
##D 
##D # Note the following will FAIL because of a name mismatch in the arguments
##D # (with an informative error):
##D consume(ep, x, retryDelay=1)
##D # We need the columns of the data frame to match the inputSchema:
##D names(x) <- letters[1:4]
##D # Now we can evaluate all the rows of the data frame in one call:
##D consume(ep, x)
##D # output should look like:
##D #    sum    prod
##D # 1 10.2   4.998
##D # 2  9.5   4.116
##D # 3  9.4  3.9104
##D # 4  9.4   4.278
##D # 5 10.2    5.04
##D # 6 11.4 14.3208
##D 
##D # You can use consume to evaluate just a single set of input values with this
##D # form:
##D consume(ep, a=1, b=2, c=3, d=4)
##D 
##D # or, equivalently,
##D consume(ep, list(a=1, b=2, c=3, d=4))
##D 
##D # You can evaluate multiple sets of input values with a data frame input:
##D consume(ep, data.frame(a=1:2, b=3:4, c=5:6, d=7:8))
##D 
##D # or, equivalently, with multiple lists:
##D consume(ep, list(a=1, b=3, c=5, d=7), list(a=2, b=4, c=6, d=8))
##D 
##D # Remove the service we just published
##D deleteWebService(ws, "rowSums")
##D 
##D # A more efficient way to do the same thing using data frame input/output:
##D f <- function(df) with(df, list(sum = a+b+c+d, prod = a*b*c*d))
##D ep = publishWebService(ws, f, name="rowSums2", 
##D                        inputSchema = data.frame(a = 0, b = 0, c = 0, d = 0))
##D consume(ep, data.frame(a=1:2, b=3:4, c=5:6, d=7:8))
##D deleteWebService(ws, "rowSums2")
##D 
##D 
##D 
##D # Automatically discover dependencies ----------------------------------------
##D 
##D # The publishWebService function uses `miniCRAN` to include dependencies on
##D # packages required by your function. The next example uses the `lmer`
##D # function from the lme4 package, and also shows how to publish a function
##D # that consumes a data frame by setting data.frame=TRUE.  Note! This example
##D # depends on a lot of packages and may take some time to upload to Azure.
##D library(lme4)
##D # Build a sample mixed effects model on just a subset of the sleepstudy data...
##D set.seed(1)
##D m <- lmer(Reaction ~ Days + (Days | Subject), 
##D           data=sleepstudy[sample(nrow(sleepstudy), 120),])
##D # Deine a prediction function to publish based on the model:
##D fun <- function(newdata)
##D {
##D   predict(m, newdata=newdata)
##D }
##D ep <- publishWebService(ws, fun=fun, name="sleepy lmer",
##D                         inputSchema= sleepstudy,
##D                         packages="lme4",
##D                         data.frame=TRUE)
##D 
##D # OK, try this out, and compare with raw data
##D ans = consume(ep, sleepstudy)$ans
##D plot(ans, sleepstudy$Reaction)
##D 
##D # Remove the service
##D deleteWebService(ws, "sleepy lmer")
## End(Not run)



cleanEx()
nameEx("download.datasets")
### * download.datasets

flush(stderr()); flush(stdout())

### Name: download.datasets
### Title: Download one or more datasets from an AzureML workspace.
### Aliases: download.datasets

### ** Examples

## Not run: 
##D   library(AzureML)
##D   
##D   name <- "Blood donation data"
##D   
##D   ws <- workspace()
##D   
##D   # The following three alternatives produce the same output:
##D   frame1 <- download.datasets(ws, name)
##D   frame2 <- download.datasets(datasets(ws), name)
##D 
##D   # Note that one can examine all the names, sizes, etc. of the datasets
##D   # in ws by examining d:
##D   d <- datasets(ws)
##D   frame3 <- download.datasets(subset(d, Name == name))
##D 
##D   head(frame1)
## End(Not run)



cleanEx()
nameEx("endpointHelp")
### * endpointHelp

flush(stderr()); flush(stdout())

### Name: endpointHelp
### Title: Display AzureML Web Service Endpoint Help Screens.
### Aliases: endpointHelp

### ** Examples

## Not run: 
##D ws <- workspace()
##D 
##D s <- services(ws)
##D e <- endpoints(ws, s[1,])
##D endpointHelp(e)
##D 
##D Particularly useful way to see expected service input and output:
##D endpointHelp(e)$definitions
##D 
## End(Not run)



cleanEx()
nameEx("endpoints")
### * endpoints

flush(stderr()); flush(stdout())

### Name: endpoints
### Title: List AzureML Web Service Endpoints
### Aliases: endpoints getEndpoints

### ** Examples

## Not run: 
##D workspace_id <- ""          # Your AzureML workspace id
##D authorization_token <- ""   # Your AsureML authorization token
##D 
##D ws <- workspace(
##D   id = workspace_id,
##D   auth = authorization_token
##D )
##D 
##D s <- services(ws)
##D endpoints(ws, s$Id[1])
##D 
##D # Note that you can alternatively just use the entire row that
##D # describes the service.
##D endpoints(ws, s[1,])
##D 
##D # Equivalent:
##D getEndpoints(ws, s$Id[1])
## End(Not run)



cleanEx()
nameEx("experiments")
### * experiments

flush(stderr()); flush(stdout())

### Name: experiments
### Title: List experiments in an AzureML workspace.
### Aliases: experiments

### ** Examples

## Not run: 
##D   library(AzureML)
##D   
##D   experiment <- "dd01c7e4a424432c9a9f83142d5cfec4.f-id.d2f351dd4cec4c06a4592ac83f7af55a"
##D   node_id <- '2a472ae1-ecb1-4f40-ae4e-cd3cecb1003f-268'
##D   
##D   ws <- workspace()
##D   
##D   ws$experiments
##D   experiments(ws)
##D   frame <- download.intermediate.dataset(ws, experiment, node_id,
##D                                          port_name = "Results dataset", 
##D                                          data_type_id = "GenericCSV")
##D   head(frame)
## End(Not run)



cleanEx()
nameEx("publishWebService")
### * publishWebService

flush(stderr()); flush(stdout())

### Name: publishWebService
### Title: Publish a function as a Microsoft Azure Web Service.
### Aliases: publishWebService updateWebService

### ** Examples

## Not run: 
##D # Use a default configuration in ~/.azureml, alternatively
##D # see help for `?workspace`.
##D 
##D ws <- workspace()
##D   
##D # Publish a simple model using the lme4::sleepdata ---------------------------
##D 
##D library(lme4)
##D set.seed(1)
##D train <- sleepstudy[sample(nrow(sleepstudy), 120),]
##D m <- lm(Reaction ~ Days + Subject, data = train)
##D 
##D # Deine a prediction function to publish based on the model:
##D sleepyPredict <- function(newdata){
##D   predict(m, newdata=newdata)
##D }
##D 
##D ep <- publishWebService(ws, fun = sleepyPredict, name="sleepy lm",
##D                         inputSchema = sleepstudy,
##D                         data.frame=TRUE)
##D 
##D # OK, try this out, and compare with raw data
##D ans <- consume(ep, sleepstudy)$ans
##D plot(ans, sleepstudy$Reaction)
##D 
##D # Remove the service
##D deleteWebService(ws, "sleepy lm")
##D 
##D 
##D 
##D # Another data frame example -------------------------------------------------
##D 
##D # If your function can consume a whole data frame at once, you can also
##D # supply data in that form, resulting in more efficient computation.
##D # The following example builds a simple linear model on a subset of the
##D # airquality data and publishes a prediction function based on the model.
##D set.seed(1)
##D m <- lm(Ozone ~ ., data=airquality[sample(nrow(airquality), 100),])
##D # Define a prediction function based on the model:
##D fun <- function(newdata)
##D {
##D   predict(m, newdata=newdata)
##D }
##D # Note the definition of inputSchema and use of the data.frame argument.
##D ep <- publishWebService(ws, fun=fun, name="Ozone",
##D                         inputSchema = airquality,
##D                         data.frame=TRUE)
##D ans <- consume(ep, airquality)$ans
##D plot(ans, airquality$Ozone)
##D deleteWebService(ws, "Ozone")
##D 
##D 
##D 
##D # Train a model using diamonds in ggplot2 ------------------------------------
##D # This example also demonstrates how to deal with factor in the data
##D 
##D data(diamonds, package="ggplot2")
##D set.seed(1)
##D train_idx = sample.int(nrow(diamonds), 30000)
##D test_idx = sample(setdiff(seq(1, nrow(diamonds)), train_idx), 500)
##D train <- diamonds[train_idx, ]
##D test  <- diamonds[test_idx, ]
##D 
##D model <- glm(price ~ carat + clarity + color + cut - 1, data = train, 
##D              family = Gamma(link = "log"))
##D 
##D diamondLevels <- diamonds[1, ]
##D 
##D # The model works reasonably well, except for some outliers
##D plot(exp(predict(model, test)) ~ test$price)
##D 
##D # Create a prediction function that converts characters correctly to factors
##D 
##D predictDiamonds <- function(x){
##D   x$cut     <- factor(x$cut,     
##D                       levels = levels(diamondLevels$cut), ordered = TRUE)
##D   x$clarity <- factor(x$clarity, 
##D                       levels = levels(diamondLevels$clarity), ordered = TRUE)
##D   x$color   <- factor(x$color,   
##D                       levels = levels(diamondLevels$color), ordered = TRUE)
##D   exp(predict(model, newdata = x))
##D }
##D 
##D 
##D # Publish the service
##D 
##D ws <- workspace()
##D ep <- publishWebService(ws, fun = predictDiamonds, name = "diamonds",
##D                         inputSchema = test,
##D                         data.frame = TRUE
##D )
##D 
##D # Consume the service
##D results <- consume(ep, test)$ans
##D plot(results ~ test$price)
##D 
##D deleteWebService(ws, "diamonds")
##D 
##D 
##D 
##D # Simple example using scalar input ------------------------------------------
##D 
##D ws <- workspace()
##D 
##D # Really simple example:
##D add <- function(x,y) x + y
##D endpoint <- publishWebService(ws, 
##D                               fun = add, 
##D                               name = "addme", 
##D                               inputSchema = list(x="numeric", 
##D                                                  y="numeric"), 
##D                               outputSchema = list(ans="numeric"))
##D consume(endpoint, list(x=pi, y=2))
##D 
##D # Now remove the web service named "addme" that we just published
##D deleteWebService(ws, "addme")
##D 
##D 
##D 
##D # Send a custom R function for evaluation in AzureML -------------------------
##D 
##D # A neat trick to evaluate any expression in the Azure ML virtual
##D # machine R session and view its output:
##D ep <- publishWebService(ws, 
##D                         fun =  function(expr) {
##D                           paste(capture.output(
##D                             eval(parse(text=expr))), collapse="\n")
##D                         },
##D                         name="commander", 
##D                         inputSchema = list(x = "character"),
##D                         outputSchema = list(ans = "character"))
##D cat(consume(ep, list(x = "getwd()"))$ans)
##D cat(consume(ep, list(x = ".packages(all=TRUE)"))$ans)
##D cat(consume(ep, list(x = "R.Version()"))$ans)
##D 
##D # Remove the service we just published
##D deleteWebService(ws, "commander")
##D 
##D 
##D 
##D # Understanding the scoping rules --------------------------------------------
##D 
##D # The following example illustrates scoping rules. Note that the function
##D # refers to the variable y defined outside the function body. That value
##D # will be exported with the service.
##D y <- pi
##D ep <- publishWebService(ws, 
##D                         fun = function(x) x + y, 
##D                         name = "lexical scope",
##D                         inputSchema = list(x = "numeric"), 
##D                         outputSchema = list(ans = "numeric"))
##D cat(consume(ep, list(x=2))$ans)
##D 
##D # Remove the service we just published
##D deleteWebService(ws, "lexical scope")
##D 
##D 
##D # Demonstrate scalar inputs but sending a data frame for scoring -------------
##D 
##D # Example showing the use of consume to score all the rows of a data frame
##D # at once, and other invocations for evaluating multiple sets of input
##D # values. The columns of the data frame correspond to the input parameters
##D # of the web service in this example:
##D f <- function(a,b,c,d) list(sum = a+b+c+d, prod = a*b*c*d)
##D ep <-  publishWebService(ws, 
##D                          f, 
##D                          name = "rowSums",
##D                          inputSchema = list(
##D                            a = "numeric", 
##D                            b = "numeric", 
##D                            c = "numeric", 
##D                            d = "numeric"
##D                          ),
##D                          outputSchema = list(
##D                            sum ="numeric", 
##D                            prod = "numeric")
##D )
##D x <- head(iris[,1:4])  # First four columns of iris
##D 
##D # Note the following will FAIL because of a name mismatch in the arguments
##D # (with an informative error):
##D consume(ep, x, retryDelay=1)
##D # We need the columns of the data frame to match the inputSchema:
##D names(x) <- letters[1:4]
##D # Now we can evaluate all the rows of the data frame in one call:
##D consume(ep, x)
##D # output should look like:
##D #    sum    prod
##D # 1 10.2   4.998
##D # 2  9.5   4.116
##D # 3  9.4  3.9104
##D # 4  9.4   4.278
##D # 5 10.2    5.04
##D # 6 11.4 14.3208
##D 
##D # You can use consume to evaluate just a single set of input values with this
##D # form:
##D consume(ep, a=1, b=2, c=3, d=4)
##D 
##D # or, equivalently,
##D consume(ep, list(a=1, b=2, c=3, d=4))
##D 
##D # You can evaluate multiple sets of input values with a data frame input:
##D consume(ep, data.frame(a=1:2, b=3:4, c=5:6, d=7:8))
##D 
##D # or, equivalently, with multiple lists:
##D consume(ep, list(a=1, b=3, c=5, d=7), list(a=2, b=4, c=6, d=8))
##D 
##D # Remove the service we just published
##D deleteWebService(ws, "rowSums")
##D 
##D # A more efficient way to do the same thing using data frame input/output:
##D f <- function(df) with(df, list(sum = a+b+c+d, prod = a*b*c*d))
##D ep = publishWebService(ws, f, name="rowSums2", 
##D                        inputSchema = data.frame(a = 0, b = 0, c = 0, d = 0))
##D consume(ep, data.frame(a=1:2, b=3:4, c=5:6, d=7:8))
##D deleteWebService(ws, "rowSums2")
##D 
##D 
##D 
##D # Automatically discover dependencies ----------------------------------------
##D 
##D # The publishWebService function uses `miniCRAN` to include dependencies on
##D # packages required by your function. The next example uses the `lmer`
##D # function from the lme4 package, and also shows how to publish a function
##D # that consumes a data frame by setting data.frame=TRUE.  Note! This example
##D # depends on a lot of packages and may take some time to upload to Azure.
##D library(lme4)
##D # Build a sample mixed effects model on just a subset of the sleepstudy data...
##D set.seed(1)
##D m <- lmer(Reaction ~ Days + (Days | Subject), 
##D           data=sleepstudy[sample(nrow(sleepstudy), 120),])
##D # Deine a prediction function to publish based on the model:
##D fun <- function(newdata)
##D {
##D   predict(m, newdata=newdata)
##D }
##D ep <- publishWebService(ws, fun=fun, name="sleepy lmer",
##D                         inputSchema= sleepstudy,
##D                         packages="lme4",
##D                         data.frame=TRUE)
##D 
##D # OK, try this out, and compare with raw data
##D ans = consume(ep, sleepstudy)$ans
##D plot(ans, sleepstudy$Reaction)
##D 
##D # Remove the service
##D deleteWebService(ws, "sleepy lmer")
## End(Not run)



cleanEx()
nameEx("services")
### * services

flush(stderr()); flush(stdout())

### Name: services
### Title: List Available Web Services.
### Aliases: services getWebServices

### ** Examples

## Not run: 
##D workspace_id <- ""          # Your AzureML workspace id
##D authorization_token <- ""   # Your AzureML authorization token
##D 
##D ws <- workspace(
##D   id = workspace_id,
##D   auth = authorization_token
##D )
##D 
##D # Equivalent:
##D services(ws)
##D getWebServices(ws)
## End(Not run)



cleanEx()
nameEx("upload.dataset")
### * upload.dataset

flush(stderr()); flush(stdout())

### Name: upload.dataset
### Title: Upload an R data frame to an AzureML workspace.
### Aliases: upload.dataset

### ** Examples

## Not run: 
##D   library(AzureML)
##D   
##D   ws <- workspace()
##D   
##D   # Upload the R airquality data.frame to the workspace.
##D   upload.dataset(airquality, ws, "airquality")
##D 
##D   # Example datasets (airquality should be among them now)
##D   head(datasets(ws))
##D 
##D   # Now delete what we've just uploaded
##D   delete.datasets(ws, "airquality")
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
