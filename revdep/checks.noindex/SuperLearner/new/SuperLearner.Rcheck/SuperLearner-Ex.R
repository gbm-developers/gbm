pkgname <- "SuperLearner"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SuperLearner')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CV.SuperLearner")
### * CV.SuperLearner

flush(stderr()); flush(stdout())

### Name: CV.SuperLearner
### Title: Function to get V-fold cross-validated risk estimate for super
###   learner
### Aliases: CV.SuperLearner print.CV.SuperLearner coef.CV.SuperLearner
### Keywords: models

### ** Examples

## Not run: 
##D set.seed(23432)
##D ## training set
##D n <- 500
##D p <- 50
##D X <- matrix(rnorm(n*p), nrow = n, ncol = p)
##D colnames(X) <- paste("X", 1:p, sep="")
##D X <- data.frame(X)
##D Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)
##D 
##D ## build Library and run Super Learner
##D SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.polymars", "SL.mean")
##D 
##D test <- CV.SuperLearner(Y = Y, X = X, V = 10, SL.library = SL.library,
##D   verbose = TRUE, method = "method.NNLS")
##D test
##D summary(test)
##D ## Look at the coefficients across folds
##D coef(test)
##D 
##D # Example with specifying cross-validation options for both 
##D # CV.SuperLearner (cvControl) and the internal SuperLearners (innerCvControl)
##D test <- CV.SuperLearner(Y = Y, X = X, SL.library = SL.library,
##D   cvControl = list(V = 10, shuffle = FALSE),
##D   innerCvControl = list(list(V = 5)),
##D   verbose = TRUE, method = "method.NNLS")
##D 
##D ## examples with snow
##D library(parallel)
##D cl <- makeCluster(2, type = "PSOCK") # can use different types here
##D clusterSetRNGStream(cl, iseed = 2343)
##D testSNOW <- CV.SuperLearner(Y = Y, X = X, SL.library = SL.library, method = "method.NNLS",
##D   parallel = cl)
##D summary(testSNOW)
##D stopCluster(cl)
## End(Not run)



cleanEx()
nameEx("SL.biglasso")
### * SL.biglasso

flush(stderr()); flush(stdout())

### Name: SL.biglasso
### Title: SL wrapper for biglasso
### Aliases: SL.biglasso

### ** Examples


data(Boston, package = "MASS")
Y = Boston$medv
# Remove outcome from covariate dataframe.
X = Boston[, -14]

set.seed(1)

# Sample rows to speed up example.
row_subset = sample(nrow(X), 30)

# Subset rows and columns & use only 2 folds to speed up example.
sl = SuperLearner(Y[row_subset], X[row_subset, 1:2, drop = FALSE],
                  family = gaussian(), cvControl = list(V = 2),
                  SL.library = "SL.biglasso")
sl

pred = predict(sl, X)
summary(pred$pred)




cleanEx()
nameEx("SL.extraTrees")
### * SL.extraTrees

flush(stderr()); flush(stdout())

### Name: SL.extraTrees
### Title: extraTrees SuperLearner wrapper
### Aliases: SL.extraTrees

### ** Examples


data(Boston, package = "MASS")
Y = Boston$medv
# Remove outcome from covariate dataframe.
X = Boston[, -14]

set.seed(1)

# Sample rows to speed up example.
row_subset = sample(nrow(X), 30)

sl = SuperLearner(Y[row_subset], X[row_subset, ], family = gaussian(),
cvControl = list(V = 2), SL.library = c("SL.mean", "SL.extraTrees"))

print(sl)




cleanEx()
nameEx("SL.glm")
### * SL.glm

flush(stderr()); flush(stdout())

### Name: SL.glm
### Title: Wrapper for glm
### Aliases: SL.glm

### ** Examples


data(Boston, package = "MASS")
Y = Boston$medv
# Remove outcome from covariate dataframe.
X = Boston[, -14]

set.seed(1)

sl = SuperLearner(Y, X, family = gaussian(),
                  SL.library = c("SL.mean", "SL.glm"))

print(sl)




cleanEx()
nameEx("SL.glmnet")
### * SL.glmnet

flush(stderr()); flush(stdout())

### Name: SL.glmnet
### Title: Elastic net regression, including lasso and ridge
### Aliases: SL.glmnet

### ** Examples


# Load a test dataset.
data(PimaIndiansDiabetes2, package = "mlbench")
data = PimaIndiansDiabetes2

# Omit observations with missing data.
data = na.omit(data)

Y = as.numeric(data$diabetes == "pos")
X = subset(data, select = -diabetes)

set.seed(1, "L'Ecuyer-CMRG")

sl = SuperLearner(Y, X, family = binomial(),
                  SL.library = c("SL.mean", "SL.glm", "SL.glmnet"))
sl




cleanEx()
nameEx("SL.kernelKnn")
### * SL.kernelKnn

flush(stderr()); flush(stdout())

### Name: SL.kernelKnn
### Title: SL wrapper for KernelKNN
### Aliases: SL.kernelKnn

### ** Examples


# Load a test dataset.
data(PimaIndiansDiabetes2, package = "mlbench")

data = PimaIndiansDiabetes2

# Omit observations with missing data.
data = na.omit(data)

Y_bin = as.numeric(data$diabetes)
X = subset(data, select = -diabetes)

set.seed(1)

sl = SuperLearner(Y_bin, X, family = binomial(),
                 SL.library = c("SL.mean", "SL.kernelKnn"))
sl




cleanEx()
nameEx("SL.ksvm")
### * SL.ksvm

flush(stderr()); flush(stdout())

### Name: SL.ksvm
### Title: Wrapper for Kernlab's SVM algorithm
### Aliases: SL.ksvm

### ** Examples


data(Boston, package = "MASS")
Y = Boston$medv
# Remove outcome from covariate dataframe.
X = Boston[, -14]

set.seed(1)

sl = SuperLearner(Y, X, family = gaussian(),
                 SL.library = c("SL.mean", "SL.ksvm"))
sl

pred = predict(sl, X)
summary(pred$pred)




cleanEx()
nameEx("SL.lda")
### * SL.lda

flush(stderr()); flush(stdout())

### Name: SL.lda
### Title: SL wrapper for MASS:lda
### Aliases: SL.lda

### ** Examples


data(Boston, package = "MASS")
Y = as.numeric(Boston$medv > 23)
# Remove outcome from covariate dataframe.
X = Boston[, -14]

set.seed(1)

# Use only 2 CV folds to speed up example.
sl = SuperLearner(Y, X, family = binomial(), cvControl = list(V = 2),
                 SL.library = c("SL.mean", "SL.lda"))
sl

pred = predict(sl, X)
summary(pred$pred)




cleanEx()
nameEx("SL.lm")
### * SL.lm

flush(stderr()); flush(stdout())

### Name: SL.lm
### Title: Wrapper for lm
### Aliases: SL.lm

### ** Examples


data(Boston, package = "MASS")
Y = Boston$medv
# Remove outcome from covariate dataframe.
X = Boston[, -14]

set.seed(1)

sl = SuperLearner(Y, X, family = gaussian(),
                  SL.library = c("SL.mean", "SL.lm"))

print(sl)




cleanEx()
nameEx("SL.qda")
### * SL.qda

flush(stderr()); flush(stdout())

### Name: SL.qda
### Title: SL wrapper for MASS:qda
### Aliases: SL.qda

### ** Examples


data(Boston, package = "MASS")
Y = as.numeric(Boston$medv > 23)
# Remove outcome from covariate dataframe.
X = Boston[, -14]

set.seed(1)

# Use only 2 CV folds to speed up example.
sl = SuperLearner(Y, X, family = binomial(), cvControl = list(V = 2),
                 SL.library = c("SL.mean", "SL.qda"))
sl

pred = predict(sl, X)
summary(pred$pred)





cleanEx()
nameEx("SL.ranger")
### * SL.ranger

flush(stderr()); flush(stdout())

### Name: SL.ranger
### Title: SL wrapper for ranger
### Aliases: SL.ranger

### ** Examples


data(Boston, package = "MASS")
Y = Boston$medv
# Remove outcome from covariate dataframe.
X = Boston[, -14]

set.seed(1)

# Use only 2 CV folds to speed up example.
sl = SuperLearner(Y, X, family = gaussian(), cvControl = list(V = 2),
                 SL.library = c("SL.mean", "SL.ranger"))
sl

pred = predict(sl, X)
summary(pred$pred)




cleanEx()
nameEx("SampleSplitSuperLearner")
### * SampleSplitSuperLearner

flush(stderr()); flush(stdout())

### Name: SampleSplitSuperLearner
### Title: Super Learner Prediction Function
### Aliases: SampleSplitSuperLearner
### Keywords: models

### ** Examples

## Not run: 
##D ## simulate data
##D set.seed(23432)
##D ## training set
##D n <- 500
##D p <- 50
##D X <- matrix(rnorm(n*p), nrow = n, ncol = p)
##D colnames(X) <- paste("X", 1:p, sep="")
##D X <- data.frame(X)
##D Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)
##D 
##D ## test set
##D m <- 1000
##D newX <- matrix(rnorm(m*p), nrow = m, ncol = p)
##D colnames(newX) <- paste("X", 1:p, sep="")
##D newX <- data.frame(newX)
##D newY <- newX[, 1] + sqrt(abs(newX[, 2] * newX[, 3])) + newX[, 2] -
##D   newX[, 3] + rnorm(m)
##D 
##D # generate Library and run Super Learner
##D SL.library <- c("SL.glm", "SL.randomForest", "SL.gam",
##D   "SL.polymars", "SL.mean")
##D test <- SampleSplitSuperLearner(Y = Y, X = X, newX = newX, SL.library = SL.library,
##D   verbose = TRUE, method = "method.NNLS")
##D test
##D 
##D # library with screening
##D SL.library <- list(c("SL.glmnet", "All"), c("SL.glm", "screen.randomForest",
##D   "All", "screen.SIS"), "SL.randomForest", c("SL.polymars", "All"), "SL.mean")
##D test <- SuperLearner(Y = Y, X = X, newX = newX, SL.library = SL.library,
##D   verbose = TRUE, method = "method.NNLS")
##D test
##D 
##D # binary outcome
##D set.seed(1)
##D N <- 200
##D X <- matrix(rnorm(N*10), N, 10)
##D X <- as.data.frame(X)
##D Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + 
##D   .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
##D 
##D SL.library <- c("SL.glmnet", "SL.glm", "SL.knn", "SL.gam", "SL.mean")
##D 
##D # least squares loss function
##D test.NNLS <- SampleSplitSuperLearner(Y = Y, X = X, SL.library = SL.library, 
##D   verbose = TRUE, method = "method.NNLS", family = binomial())
##D test.NNLS
## End(Not run)



cleanEx()
nameEx("SuperLearner")
### * SuperLearner

flush(stderr()); flush(stdout())

### Name: SuperLearner
### Title: Super Learner Prediction Function
### Aliases: SuperLearner mcSuperLearner snowSuperLearner
###   print.SuperLearner coef.SuperLearner
### Keywords: models

### ** Examples

## Not run: 
##D ## simulate data
##D set.seed(23432)
##D ## training set
##D n <- 500
##D p <- 50
##D X <- matrix(rnorm(n*p), nrow = n, ncol = p)
##D colnames(X) <- paste("X", 1:p, sep="")
##D X <- data.frame(X)
##D Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)
##D 
##D ## test set
##D m <- 1000
##D newX <- matrix(rnorm(m*p), nrow = m, ncol = p)
##D colnames(newX) <- paste("X", 1:p, sep="")
##D newX <- data.frame(newX)
##D newY <- newX[, 1] + sqrt(abs(newX[, 2] * newX[, 3])) + newX[, 2] -
##D   newX[, 3] + rnorm(m)
##D 
##D # generate Library and run Super Learner
##D SL.library <- c("SL.glm", "SL.randomForest", "SL.gam",
##D   "SL.polymars", "SL.mean")
##D test <- SuperLearner(Y = Y, X = X, newX = newX, SL.library = SL.library,
##D   verbose = TRUE, method = "method.NNLS")
##D test
##D 
##D # library with screening
##D SL.library <- list(c("SL.glmnet", "All"), c("SL.glm", "screen.randomForest",
##D   "All", "screen.SIS"), "SL.randomForest", c("SL.polymars", "All"), "SL.mean")
##D test <- SuperLearner(Y = Y, X = X, newX = newX, SL.library = SL.library,
##D   verbose = TRUE, method = "method.NNLS")
##D test
##D 
##D # binary outcome
##D set.seed(1)
##D N <- 200
##D X <- matrix(rnorm(N*10), N, 10)
##D X <- as.data.frame(X)
##D Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] +
##D   .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
##D 
##D SL.library <- c("SL.glmnet", "SL.glm", "SL.knn", "SL.gam", "SL.mean")
##D 
##D # least squares loss function
##D test.NNLS <- SuperLearner(Y = Y, X = X, SL.library = SL.library,
##D   verbose = TRUE, method = "method.NNLS", family = binomial())
##D test.NNLS
##D 
##D # negative log binomial likelihood loss function
##D test.NNloglik <- SuperLearner(Y = Y, X = X, SL.library = SL.library,
##D   verbose = TRUE, method = "method.NNloglik", family = binomial())
##D test.NNloglik
##D 
##D # 1 - AUC loss function
##D test.AUC <- SuperLearner(Y = Y, X = X, SL.library = SL.library,
##D   verbose = TRUE, method = "method.AUC", family = binomial())
##D test.AUC
##D 
##D # 2
##D # adapted from library(SIS)
##D set.seed(1)
##D # training
##D b <- c(2, 2, 2, -3*sqrt(2))
##D n <- 150
##D p <- 200
##D truerho <- 0.5
##D corrmat <- diag(rep(1-truerho, p)) + matrix(truerho, p, p)
##D corrmat[, 4] = sqrt(truerho)
##D corrmat[4, ] = sqrt(truerho)
##D corrmat[4, 4] = 1
##D cholmat <- chol(corrmat)
##D x <- matrix(rnorm(n*p, mean=0, sd=1), n, p)
##D x <- x ##D 
##D feta <- x[, 1:4] ##D 
##D fprob <- exp(feta) / (1 + exp(feta))
##D y <- rbinom(n, 1, fprob)
##D 
##D # test
##D m <- 10000
##D newx <- matrix(rnorm(m*p, mean=0, sd=1), m, p)
##D newx <- newx ##D 
##D newfeta <- newx[, 1:4] ##D 
##D newfprob <- exp(newfeta) / (1 + exp(newfeta))
##D newy <- rbinom(m, 1, newfprob)
##D 
##D DATA2 <- data.frame(Y = y, X = x)
##D newDATA2 <- data.frame(Y = newy, X=newx)
##D 
##D create.SL.knn <- function(k = c(20, 30)) {
##D   for(mm in seq(length(k))){
##D     eval(parse(text = paste('SL.knn.', k[mm], '<- function(..., k = ', k[mm],
##D       ') SL.knn(..., k = k)', sep = '')), envir = .GlobalEnv)
##D   }
##D   invisible(TRUE)
##D }
##D create.SL.knn(c(20, 30, 40, 50, 60, 70))
##D 
##D # library with screening
##D SL.library <- list(c("SL.glmnet", "All"), c("SL.glm", "screen.randomForest"),
##D   "SL.randomForest", "SL.knn", "SL.knn.20", "SL.knn.30", "SL.knn.40",
##D   "SL.knn.50", "SL.knn.60", "SL.knn.70",
##D   c("SL.polymars", "screen.randomForest"))
##D test <- SuperLearner(Y = DATA2$Y, X = DATA2[, -1], newX = newDATA2[, -1],
##D   SL.library = SL.library, verbose = TRUE, family = binomial())
##D test
##D 
##D ## examples with multicore
##D set.seed(23432, "L'Ecuyer-CMRG")  # use L'Ecuyer for multicore seeds. see ?set.seed for details
##D ## training set
##D n <- 500
##D p <- 50
##D X <- matrix(rnorm(n*p), nrow = n, ncol = p)
##D colnames(X) <- paste("X", 1:p, sep="")
##D X <- data.frame(X)
##D Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)
##D 
##D ## test set
##D m <- 1000
##D newX <- matrix(rnorm(m*p), nrow = m, ncol = p)
##D colnames(newX) <- paste("X", 1:p, sep="")
##D newX <- data.frame(newX)
##D newY <- newX[, 1] + sqrt(abs(newX[, 2] * newX[, 3])) + newX[, 2] - newX[, 3] + rnorm(m)
##D 
##D # generate Library and run Super Learner
##D SL.library <- c("SL.glm", "SL.randomForest", "SL.gam",
##D   "SL.polymars", "SL.mean")
##D 
##D testMC <- mcSuperLearner(Y = Y, X = X, newX = newX, SL.library = SL.library,
##D   method = "method.NNLS")
##D testMC
##D 
##D ## examples with snow
##D library(parallel)
##D cl <- makeCluster(2, type = "PSOCK") # can use different types here
##D clusterSetRNGStream(cl, iseed = 2343)
##D # make SL functions available on the clusters, use assignment to avoid printing
##D foo <- clusterEvalQ(cl, library(SuperLearner))  
##D testSNOW <- snowSuperLearner(cluster = cl, Y = Y, X = X, newX = newX,
##D   SL.library = SL.library, method = "method.NNLS")
##D testSNOW
##D stopCluster(cl)
##D 
##D ## snow example with user-generated wrappers
##D # If you write your own wrappers and are using snowSuperLearner()
##D # These new wrappers need to be added to the SuperLearner namespace and exported to the clusters
##D # Using a simple example here, but can define any new SuperLearner wrapper
##D my.SL.wrapper <- function(...) SL.glm(...)
##D # assign function into SuperLearner namespace
##D environment(my.SL.wrapper) <-asNamespace("SuperLearner")
##D 
##D cl <- makeCluster(2, type = "PSOCK") # can use different types here
##D clusterSetRNGStream(cl, iseed = 2343)
##D # make SL functions available on the clusters, use assignment to avoid printing	
##D foo <- clusterEvalQ(cl, library(SuperLearner))  
##D clusterExport(cl, c("my.SL.wrapper"))  # copy the function to all clusters
##D testSNOW <- snowSuperLearner(cluster = cl, Y = Y, X = X, newX = newX,
##D   SL.library = c("SL.glm", "SL.mean", "my.SL.wrapper"), method = "method.NNLS")
##D testSNOW
##D stopCluster(cl)
##D 
##D ## timing
##D replicate(5, system.time(SuperLearner(Y = Y, X = X, newX = newX,
##D   SL.library = SL.library, method = "method.NNLS")))
##D 
##D replicate(5, system.time(mcSuperLearner(Y = Y, X = X, newX = newX,
##D   SL.library = SL.library, method = "method.NNLS")))
##D 
##D cl <- makeCluster(2, type = 'PSOCK')
##D # make SL functions available on the clusters, use assignment to avoid printing	
##D foo <- clusterEvalQ(cl, library(SuperLearner))  
##D replicate(5, system.time(snowSuperLearner(cl, Y = Y, X = X, newX = newX,
##D   SL.library = SL.library, method = "method.NNLS")))
##D stopCluster(cl)
##D 
## End(Not run)



cleanEx()
nameEx("create.Learner")
### * create.Learner

flush(stderr()); flush(stdout())

### Name: create.Learner
### Title: Factory for learner wrappers
### Aliases: create.Learner

### ** Examples

## Not run: 
##D # Create a randomForest learner with ntree set to 1000 rather than the
##D # default of 500.
##D create_rf = create.Learner("SL.randomForest", list(ntree = 1000))
##D create_rf
##D sl = SuperLearner(Y = Y, X = X, SL.library = create_rf$names, family = binomial())
##D sl
##D # Clean up global environment.
##D rm(list = create_rf$names)
##D # Create a randomForest learner that optimizes over mtry
##D create_rf = create.Learner("SL.randomForest",
##D                      tune = list(mtry = round(c(1, sqrt(ncol(X)), ncol(X)))))
##D create_rf
##D sl = SuperLearner(Y = Y, X = X, SL.library = create_rf$names, family = binomial())
##D sl
##D # Clean up global environment.
##D rm(list = create_rf$names)
##D 
##D # Optimize elastic net over alpha, with a custom environment and detailed names.
##D learners = new.env()
##D create_enet = create.Learner("SL.glmnet", env = learners, detailed_names = T,
##D                            tune = list(alpha = seq(0, 1, length.out=5)))
##D create_enet
##D # List the environment to review what functions were created.
##D ls(learners)
##D # We can simply list the environment to specify the library.
##D sl = SuperLearner(Y = Y, X = X, SL.library = ls(learners), family = binomial(), env = learners)
##D sl
## End(Not run)




cleanEx()
nameEx("create.SL.xgboost")
### * create.SL.xgboost

flush(stderr()); flush(stdout())

### Name: create.SL.xgboost
### Title: Factory for XGBoost SL wrappers
### Aliases: create.SL.xgboost

### ** Examples


# Create a new environment to store the learner functions.
# This keeps the global environment organized.
sl_env = new.env()
# Create 2 * 2 * 1 * 3 = 12 combinations of hyperparameters.
tune = list(ntrees = c(100, 500), max_depth = c(1, 2), minobspernode = 10,
            shrinkage = c(0.1, 0.01, 0.001))
# Generate a separate learner for each combination.
xgb_grid = create.SL.xgboost(tune = tune, env = sl_env)
# Review the function configurations.
xgb_grid
# Attach the environment so that the custom learner functions can be accessed.
attach(sl_env)
## Not run: 
##D sl = SuperLearner(Y = Y, X = X, SL.library = xgb_grid$names)
## End(Not run)
detach(sl_env)



cleanEx()
nameEx("listWrappers")
### * listWrappers

flush(stderr()); flush(stdout())

### Name: listWrappers
### Title: list all wrapper functions in SuperLearner
### Aliases: listWrappers
### Keywords: utilities

### ** Examples

listWrappers(what = "SL")
listWrappers(what = "screen")



cleanEx()
nameEx("recombineCVSL")
### * recombineCVSL

flush(stderr()); flush(stdout())

### Name: recombineCVSL
### Title: Recombine a CV.SuperLearner fit using a new metalearning method
### Aliases: recombineCVSL
### Keywords: models

### ** Examples

## Not run: 
##D 
##D # Binary outcome example adapted from SuperLearner examples
##D 
##D set.seed(1)
##D N <- 200
##D X <- matrix(rnorm(N*10), N, 10)
##D X <- as.data.frame(X)
##D Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + 
##D   .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
##D 
##D SL.library <- c("SL.glmnet", "SL.glm", "SL.knn", "SL.gam", "SL.mean")
##D 
##D # least squares loss function
##D set.seed(1) # for reproducibility
##D cvfit_nnls <- CV.SuperLearner(Y = Y, X = X, V = 10, SL.library = SL.library, 
##D   verbose = TRUE, method = "method.NNLS", family = binomial())
##D cvfit_nnls$coef
##D #    SL.glmnet_All SL.glm_All  SL.knn_All SL.gam_All SL.mean_All
##D # 1      0.0000000 0.00000000 0.000000000  0.4143862   0.5856138
##D # 2      0.0000000 0.00000000 0.304802397  0.3047478   0.3904498
##D # 3      0.0000000 0.00000000 0.002897533  0.5544075   0.4426950
##D # 4      0.0000000 0.20322642 0.000000000  0.1121891   0.6845845
##D # 5      0.1743973 0.00000000 0.032471026  0.3580624   0.4350693
##D # 6      0.0000000 0.00000000 0.099881535  0.3662309   0.5338876
##D # 7      0.0000000 0.00000000 0.234876082  0.2942472   0.4708767
##D # 8      0.0000000 0.06424676 0.113988158  0.5600208   0.2617443
##D # 9      0.0000000 0.00000000 0.338030342  0.2762604   0.3857093
##D # 10     0.3022442 0.00000000 0.294226204  0.1394534   0.2640762
##D 
##D 
##D # negative log binomial likelihood loss function
##D cvfit_nnloglik <- recombineCVSL(cvfit_nnls, method = "method.NNloglik")
##D cvfit_nnloglik$coef
##D #    SL.glmnet_All SL.glm_All SL.knn_All SL.gam_All SL.mean_All
##D # 1      0.0000000  0.0000000 0.00000000  0.5974799  0.40252010
##D # 2      0.0000000  0.0000000 0.31177345  0.6882266  0.00000000
##D # 3      0.0000000  0.0000000 0.01377469  0.8544238  0.13180152
##D # 4      0.0000000  0.1644188 0.00000000  0.2387919  0.59678930
##D # 5      0.2142254  0.0000000 0.00000000  0.3729426  0.41283197
##D # 6      0.0000000  0.0000000 0.00000000  0.5847150  0.41528502
##D # 7      0.0000000  0.0000000 0.47538172  0.5080311  0.01658722
##D # 8      0.0000000  0.0000000 0.00000000  1.0000000  0.00000000
##D # 9      0.0000000  0.0000000 0.45384961  0.2923480  0.25380243
##D # 10     0.3977816  0.0000000 0.27927906  0.1606384  0.16230097
##D 
##D # If we use the same seed as the original `cvfit_nnls`, then
##D # the recombineCVSL and CV.SuperLearner results will be identical
##D # however, the recombineCVSL version will be much faster since
##D # it doesn't have to re-fit all the base learners, V times each.
##D set.seed(1)
##D cvfit_nnloglik2 <- CV.SuperLearner(Y = Y, X = X, V = 10, SL.library = SL.library,
##D   verbose = TRUE, method = "method.NNloglik", family = binomial())
##D cvfit_nnloglik2$coef
##D #    SL.glmnet_All SL.glm_All SL.knn_All SL.gam_All SL.mean_All
##D # 1      0.0000000  0.0000000 0.00000000  0.5974799  0.40252010
##D # 2      0.0000000  0.0000000 0.31177345  0.6882266  0.00000000
##D # 3      0.0000000  0.0000000 0.01377469  0.8544238  0.13180152
##D # 4      0.0000000  0.1644188 0.00000000  0.2387919  0.59678930
##D # 5      0.2142254  0.0000000 0.00000000  0.3729426  0.41283197
##D # 6      0.0000000  0.0000000 0.00000000  0.5847150  0.41528502
##D # 7      0.0000000  0.0000000 0.47538172  0.5080311  0.01658722
##D # 8      0.0000000  0.0000000 0.00000000  1.0000000  0.00000000
##D # 9      0.0000000  0.0000000 0.45384961  0.2923480  0.25380243
##D # 10     0.3977816  0.0000000 0.27927906  0.1606384  0.16230097
##D 
## End(Not run)



cleanEx()
nameEx("recombineSL")
### * recombineSL

flush(stderr()); flush(stdout())

### Name: recombineSL
### Title: Recombine a SuperLearner fit using a new metalearning method
### Aliases: recombineSL
### Keywords: models

### ** Examples

## Not run: 
##D 
##D # Binary outcome example adapted from SuperLearner examples
##D 
##D set.seed(1)
##D N <- 200
##D X <- matrix(rnorm(N*10), N, 10)
##D X <- as.data.frame(X)
##D Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + 
##D   .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
##D 
##D SL.library <- c("SL.glmnet", "SL.glm", "SL.knn", "SL.gam", "SL.mean")
##D 
##D # least squares loss function
##D set.seed(1) # for reproducibility
##D fit_nnls <- SuperLearner(Y = Y, X = X, SL.library = SL.library, 
##D   verbose = TRUE, method = "method.NNLS", family = binomial())
##D fit_nnls
##D #                    Risk       Coef
##D # SL.glmnet_All 0.2439433 0.01293059
##D # SL.glm_All    0.2461245 0.08408060
##D # SL.knn_All    0.2604000 0.09600353
##D # SL.gam_All    0.2471651 0.40761918
##D # SL.mean_All   0.2486049 0.39936611
##D 
##D 
##D # negative log binomial likelihood loss function
##D fit_nnloglik <- recombineSL(fit_nnls, Y = Y, method = "method.NNloglik")
##D fit_nnloglik
##D #                    Risk      Coef
##D # SL.glmnet_All 0.6815911 0.1577228
##D # SL.glm_All    0.6918926 0.0000000
##D # SL.knn_All          Inf 0.0000000
##D # SL.gam_All    0.6935383 0.6292881
##D # SL.mean_All   0.6904050 0.2129891
##D 
##D # If we use the same seed as the original `fit_nnls`, then
##D # the recombineSL and SuperLearner results will be identical
##D # however, the recombineSL version will be much faster since
##D # it doesn't have to re-fit all the base learners.
##D set.seed(1)
##D fit_nnloglik2 <- SuperLearner(Y = Y, X = X, SL.library = SL.library,
##D   verbose = TRUE, method = "method.NNloglik", family = binomial())
##D fit_nnloglik2
##D #                    Risk      Coef
##D # SL.glmnet_All 0.6815911 0.1577228
##D # SL.glm_All    0.6918926 0.0000000
##D # SL.knn_All          Inf 0.0000000
##D # SL.gam_All    0.6935383 0.6292881
##D # SL.mean_All   0.6904050 0.2129891
##D 
## End(Not run)



cleanEx()
nameEx("trimLogit")
### * trimLogit

flush(stderr()); flush(stdout())

### Name: trimLogit
### Title: truncated-probabilities logit transformation
### Aliases: trimLogit
### Keywords: models

### ** Examples

x <- c(0.00000001, 0.0001, 0.001, 0.01, 0.1, 0.3, 0.7, 0.9, 0.99, 
  0.999, 0.9999, 0.99999999)
trimLogit(x, trim = 0.001)
data.frame(Prob = x, Logit = qlogis(x), trimLogit = trimLogit(x, 0.001))



cleanEx()
nameEx("write.SL.template")
### * write.SL.template

flush(stderr()); flush(stdout())

### Name: write.SL.template
### Title: Wrapper functions for prediction algorithms in SuperLearner
### Aliases: write.SL.template SL.template predict.SL.template SL.bayesglm
###   predict.SL.bayesglm SL.caret predict.SL.caret SL.caret.rpart
###   predict.SL.cforest SL.earth predict.SL.earth SL.gam predict.SL.gam
###   SL.gbm predict.SL.gbm SL.glm.interaction SL.ipredbagg
###   predict.SL.ipredbagg SL.knn predict.SL.knn SL.loess predict.SL.loess
###   SL.logreg predict.SL.logreg SL.mean predict.SL.mean SL.nnet
###   predict.SL.nnet SL.polymars predict.SL.polymars SL.randomForest
###   predict.SL.randomForest SL.rpart SL.rpartPrune predict.SL.rpart
###   SL.step predict.SL.step SL.step.forward SL.step.interaction
###   SL.stepAIC predict.SL.stepAIC SL.svm predict.SL.svm SL.ridge
###   predict.SL.ridge SL.leekasso predict.SL.leekasso SL.nnls
###   predict.SL.nnls
### Keywords: utilities

### ** Examples

write.SL.template(file = '')



cleanEx()
nameEx("write.method.template")
### * write.method.template

flush(stderr()); flush(stdout())

### Name: write.method.template
### Title: Method to estimate the coefficients for the super learner
### Aliases: write.method.template method.template method.NNLS method.NNLS2
###   method.NNloglik method.CC_LS method.CC_nloglik method.AUC
### Keywords: utilities

### ** Examples

write.method.template(file = '')



cleanEx()
nameEx("write.screen.template")
### * write.screen.template

flush(stderr()); flush(stdout())

### Name: write.screen.template
### Title: screening algorithms for SuperLearner
### Aliases: write.screen.template screen.template All screen.randomForest
###   screen.SIS screen.ttest screen.corP screen.corRank screen.glmnet
### Keywords: utilities

### ** Examples

write.screen.template(file = '')



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
