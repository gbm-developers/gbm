pkgname <- "imputeR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('imputeR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CubistR")
### * CubistR

flush(stderr()); flush(stdout())

### Name: CubistR
### Title: Cubist method for imputation
### Aliases: CubistR

### ** Examples

data(parkinson)
missdata <- SimIm(parkinson, 0.1)



cleanEx()
nameEx("Detect")
### * Detect

flush(stderr()); flush(stdout())

### Name: Detect
### Title: Detect variable type in a data matrix
### Aliases: Detect

### ** Examples

data(parkinson)
Detect(parkinson)
data(spect)
Detect(spect)
data(tic)
table(Detect(tic))



cleanEx()
nameEx("Rmse")
### * Rmse

flush(stderr()); flush(stdout())

### Name: Rmse
### Title: calculate the RMSE or NRMSE
### Aliases: Rmse

### ** Examples

data(parkinson)
# introduce 10% random missing values into the parkinson data
missdata <- SimIm(parkinson, 0.1)

# impute the missing values by LASSO



cleanEx()
nameEx("SimEval")
### * SimEval

flush(stderr()); flush(stdout())

### Name: SimEval
### Title: Evaluate imputation performance by simulation
### Aliases: SimEval

### ** Examples

data(parkinson)
# WARNING: simulation may take considerable time.



cleanEx()
nameEx("SimIm")
### * SimIm

flush(stderr()); flush(stdout())

### Name: SimIm
### Title: Introduce some missing values into a data matrix
### Aliases: SimIm

### ** Examples

# Create data without missing values as example
simdata <- matrix(rnorm(100), 10, 10)

# Now let's introduce some missing values into the dataset
missingdata <- SimIm(simdata, p = 0.15)

# count the number of missing values afterwards
sum(is.na(missingdata))

#------------------

# There is no missing values in the original parkinson data
data(parkinson)

# Let's introduce some missing values into the dataset
missdata <- SimIm(parkinson, 0.1)

# count the number of missing values afterwards
sum(is.na(missdata))



cleanEx()
nameEx("gbmC")
### * gbmC

flush(stderr()); flush(stdout())

### Name: gbmC
### Title: boosting tree for imputation
### Aliases: gbmC

### ** Examples

data(spect)
missdata <- SimIm(spect, 0.1)



cleanEx()
nameEx("glmboostR")
### * glmboostR

flush(stderr()); flush(stdout())

### Name: glmboostR
### Title: Boosting for regression
### Aliases: glmboostR

### ** Examples

data(parkinson)
missdata <- SimIm(parkinson, 0.1)



cleanEx()
nameEx("guess")
### * guess

flush(stderr()); flush(stdout())

### Name: guess
### Title: Impute by (educated) guessing
### Aliases: guess

### ** Examples

data(parkinson)
# introduce some random missing values
missdata <- SimIm(parkinson, 0.1)
# impute by mean imputation
impdata <- guess(missdata)
# caculate the NRMSE
Rmse(impdata, missdata, parkinson, norm = TRUE)
# by random guessing, the NRMSE should be much bigger
impdata2 <- guess(missdata, "random")
Rmse(impdata2, missdata, parkinson, norm = TRUE)



cleanEx()
nameEx("impute")
### * impute

flush(stderr()); flush(stdout())

### Name: impute
### Title: General Imputation Framework in R
### Aliases: impute

### ** Examples

data(parkinson)
# introduce 10% random missing values into the parkinson data
missdata <- SimIm(parkinson, 0.1)
# impute the missing values by LASSO



cleanEx()
nameEx("lassoC")
### * lassoC

flush(stderr()); flush(stdout())

### Name: lassoC
### Title: logistic regression with lasso for imputation
### Aliases: lassoC

### ** Examples

data(spect)
missdata <- SimIm(spect, 0.1)



cleanEx()
nameEx("lassoR")
### * lassoR

flush(stderr()); flush(stdout())

### Name: lassoR
### Title: LASSO for regression
### Aliases: lassoR

### ** Examples

data(parkinson)
missdata <- SimIm(parkinson, 0.1)



cleanEx()
nameEx("major")
### * major

flush(stderr()); flush(stdout())

### Name: major
### Title: Majority imputation for a vector
### Aliases: major

### ** Examples

a <- c(rep(0, 10), rep(1, 15), rep(2, 5))
a[sample(seq_along(a), 5)] <- NA
a
b <- major(a)
b



cleanEx()
nameEx("mixError")
### * mixError

flush(stderr()); flush(stdout())

### Name: mixError
### Title: Calculate mixed error when the imputed matrix is mixed type
### Aliases: mixError

### ** Examples

data(tic)
Detect(tic)
missdata <- SimIm(tic, 0.3)



cleanEx()
nameEx("mixGuess")
### * mixGuess

flush(stderr()); flush(stdout())

### Name: mixGuess
### Title: Naive imputation for mixed type data
### Aliases: mixGuess

### ** Examples

data(tic)
missdata <- SimIm(tic, 0.1)
sum(is.na(missdata))
impdata <- mixGuess(missdata)
sum(is.na(impdata))



cleanEx()
nameEx("mr")
### * mr

flush(stderr()); flush(stdout())

### Name: mr
### Title: calculate miss-classification error
### Aliases: mr

### ** Examples

data(spect)
Detect(spect)
missdata <- SimIm(spect, 0.1)



cleanEx()
nameEx("orderbox")
### * orderbox

flush(stderr()); flush(stdout())

### Name: orderbox
### Title: Ordered boxplot for a data matrix
### Aliases: orderbox

### ** Examples

data(parkinson)



cleanEx()
nameEx("pcrR")
### * pcrR

flush(stderr()); flush(stdout())

### Name: pcrR
### Title: Principle component regression for imputation
### Aliases: pcrR

### ** Examples

data(parkinson)
missdata <- SimIm(parkinson, 0.1)



cleanEx()
nameEx("plotIm")
### * plotIm

flush(stderr()); flush(stdout())

### Name: plotIm
### Title: Plot function for imputation
### Aliases: plotIm

### ** Examples

data(parkinson)
# introduce 10% random missing values into the parkinson data
missdata <- SimIm(parkinson, 0.1)

# impute the missing values by LASSO



cleanEx()
nameEx("plsR")
### * plsR

flush(stderr()); flush(stdout())

### Name: plsR
### Title: Partial Least Square regression for imputation
### Aliases: plsR

### ** Examples

data(parkinson)
missdata <- SimIm(parkinson, 0.1)



cleanEx()
nameEx("ridgeC")
### * ridgeC

flush(stderr()); flush(stdout())

### Name: ridgeC
### Title: Ridge regression with lasso for imputation
### Aliases: ridgeC

### ** Examples

data(spect)
missdata <- SimIm(spect, 0.1)



cleanEx()
nameEx("ridgeR")
### * ridgeR

flush(stderr()); flush(stdout())

### Name: ridgeR
### Title: Ridge shrinkage for regression
### Aliases: ridgeR

### ** Examples

data(parkinson)
missdata <- SimIm(parkinson, 0.1)



cleanEx()
nameEx("rpartC")
### * rpartC

flush(stderr()); flush(stdout())

### Name: rpartC
### Title: classification tree for imputation
### Aliases: rpartC

### ** Examples

data(spect)
missdata <- SimIm(spect, 0.1)



cleanEx()
nameEx("stepBackC")
### * stepBackC

flush(stderr()); flush(stdout())

### Name: stepBackC
### Title: Best subset for classification (backward)
### Aliases: stepBackC

### ** Examples

data(spect)
missdata <- SimIm(spect, 0.1)



cleanEx()
nameEx("stepBackR")
### * stepBackR

flush(stderr()); flush(stdout())

### Name: stepBackR
### Title: Best subset (backward direction) for regression
### Aliases: stepBackR

### ** Examples

data(parkinson)
missdata <- SimIm(parkinson, 0.1)



cleanEx()
nameEx("stepBothC")
### * stepBothC

flush(stderr()); flush(stdout())

### Name: stepBothC
### Title: Best subset for classification (both direction)
### Aliases: stepBothC

### ** Examples

data(spect)
missdata <- SimIm(spect, 0.1)



cleanEx()
nameEx("stepBothR")
### * stepBothR

flush(stderr()); flush(stdout())

### Name: stepBothR
### Title: Best subset for regression (both direction)
### Aliases: stepBothR

### ** Examples

data(parkinson)
missdata <- SimIm(parkinson, 0.1)



cleanEx()
nameEx("stepForC")
### * stepForC

flush(stderr()); flush(stdout())

### Name: stepForC
### Title: Best subset for classification (forward direction)
### Aliases: stepForC

### ** Examples

data(spect)
missdata <- SimIm(spect, 0.1)



cleanEx()
nameEx("stepForR")
### * stepForR

flush(stderr()); flush(stdout())

### Name: stepForR
### Title: Best subset (forward direction) for regression
### Aliases: stepForR

### ** Examples

data(parkinson)
missdata <- SimIm(parkinson, 0.1)



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
