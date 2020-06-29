pkgname <- "EZtune"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('EZtune')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("eztune")
### * eztune

flush(stderr()); flush(stdout())

### Name: eztune
### Title: Supervised Learning Function
### Aliases: eztune

### ** Examples

library(mlbench)
data(Sonar)
sonar <- Sonar[sample(1:nrow(Sonar), 100), ]

y <- sonar[, 61]
x <- sonar[, 1:10]

# Optimize an SVM using the default fast setting and Hooke-Jeeves
eztune(x, y)

# Optimize an SVM with 3-fold cross validation and Hooke-Jeeves
eztune(x, y, fast = FALSE, cross = 3)

# Optimize GBM using training set of 50 observations and Hooke-Jeeves
eztune(x, y, method = "gbm", fast = 50)

# Optimize SVM with 25% of the observations as a training dataset
# using a genetic algorithm
eztune(x, y, method = "svm", optimizer = "ga", fast = 0.25)




cleanEx()
nameEx("eztune_cv")
### * eztune_cv

flush(stderr()); flush(stdout())

### Name: eztune_cv
### Title: Cross Validated Accuracy for Supervised Learning Model
### Aliases: eztune_cv

### ** Examples

library(mlbench)
data(Sonar)
sonar <- Sonar[sample(1:nrow(Sonar), 100), ]

y <- sonar[, 61]
x <- sonar[, 1:10]

sonar_default <- eztune(x, y)
eztune_cv(x, y, sonar_default)

sonar_svm <- eztune(x, y, fast = FALSE, cross = 3)
eztune_cv(x, y, sonar_svm)

sonar_gbm <- eztune(x, y, method = "gbm", fast = 50)
eztune_cv(x, y, sonar_gbm)





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
