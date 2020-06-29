pkgname <- "EnsembleBase"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('EnsembleBase')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ALL.Regression.FitObj-class")
### * ALL.Regression.FitObj-class

flush(stderr()); flush(stdout())

### Name: ALL.Regression.FitObj-class
### Title: Classes '"KNN.Regression.FitObj"', '"NNET.Regression.FitObj"',
###   '"RF.Regression.FitObj"', '"SVM.Regression.FitObj"',
###   '"GBM.Regression.FitObj"', '"PENREG.Regression.FitObj"',
###   '"BART.Regression.FitObj"'
### Aliases: KNN.Regression.FitObj-class NNET.Regression.FitObj-class
###   RF.Regression.FitObj-class SVM.Regression.FitObj-class
###   GBM.Regression.FitObj-class PENREG.Regression.FitObj-class
###   BART.Regression.FitObj-class
### Keywords: classes

### ** Examples

showClass("KNN.Regression.FitObj")



cleanEx()
nameEx("BaseLearner.Batch.FitObj-class")
### * BaseLearner.Batch.FitObj-class

flush(stderr()); flush(stdout())

### Name: BaseLearner.Batch.FitObj-class
### Title: Classes '"BaseLearner.Batch.FitObj"' and
###   '"Regression.Batch.FitObj"'
### Aliases: BaseLearner.Batch.FitObj-class Regression.Batch.FitObj-class
### Keywords: classes

### ** Examples

showClass("BaseLearner.Batch.FitObj")



cleanEx()
nameEx("BaseLearner.Config-class")
### * BaseLearner.Config-class

flush(stderr()); flush(stdout())

### Name: BaseLearner.Config-class
### Title: Classes '"BaseLearner.Config"', '"Regression.Config"'
### Aliases: BaseLearner.Config-class Regression.Config-class
### Keywords: classes

### ** Examples

showClass("BaseLearner.Config")



cleanEx()
nameEx("Instance-class")
### * Instance-class

flush(stderr()); flush(stdout())

### Name: Instance-class
### Title: Classes '"Instance"' and '"Instance.List"'
### Aliases: Instance-class Instance.List-class
### Keywords: classes

### ** Examples

showClass("Instance")



cleanEx()
nameEx("Regression.Batch.Fit")
### * Regression.Batch.Fit

flush(stderr()); flush(stdout())

### Name: Regression.Batch.Fit
### Title: Batch Training, Prediction and Diagnostics of Regression Base
###   Learners
### Aliases: Regression.Batch.Fit predict.Regression.Batch.FitObj
###   plot.Regression.Batch.FitObj

### ** Examples

data(servo)
myformula <- class~motor+screw+pgain+vgain
myconfigs <- make.configs("knn")
perc.train <- 0.7
index.train <- sample(1:nrow(servo), size = round(perc.train*nrow(servo)))
data.train <- servo[index.train,]
data.predict <- servo[-index.train,]
ret <- Regression.Batch.Fit(myconfigs, myformula, data.train, ncores=2)
newpred <- predict(ret, data.predict)



cleanEx()
nameEx("Regression.CV.Batch.Fit")
### * Regression.CV.Batch.Fit

flush(stderr()); flush(stdout())

### Name: Regression.CV.Batch.Fit
### Title: CV Batch Training and Diagnostics of Regression Base Learners
### Aliases: Regression.CV.Batch.Fit predict.Regression.CV.Batch.FitObj
###   plot.Regression.CV.Batch.FitObj

### ** Examples

data(servo)
myformula <- class~motor+screw+pgain+vgain

perc.train <- 0.7
index.train <- sample(1:nrow(servo)
  , size = round(perc.train*nrow(servo)))
data.train <- servo[index.train,]
data.predict <- servo[-index.train,]

parts <- generate.partitions(1, nrow(data.train))
myconfigs <- make.configs("knn"
  , config.df = expand.grid(kernel = "rectangular", k = c(5, 10)))
instances <- make.instances(myconfigs, parts)

ret <- Regression.CV.Batch.Fit(instances, myformula, data.train)
newpred <- predict(ret, data.predict)



cleanEx()
nameEx("Regression.CV.Fit")
### * Regression.CV.Fit

flush(stderr()); flush(stdout())

### Name: Regression.CV.Fit
### Title: Cross-Validated Training and Prediction of Regression Base
###   Learners
### Aliases: Regression.CV.Fit predict.Regression.CV.FitObj

### ** Examples

data(servo)
myformula <- class~motor+screw+pgain+vgain
myconfig <- make.configs("knn", config.df=data.frame(kernel="rectangular", k=10))
perc.train <- 0.7
index.train <- sample(1:nrow(servo), size = round(perc.train*nrow(servo)))
data.train <- servo[index.train,]
data.predict <- servo[-index.train,]
mypartition <- generate.partition(nrow(data.train),nfold=3)
ret <- Regression.CV.Fit(myconfig[[1]], myformula, data.train, mypartition)
newpred <- predict(ret, data.predict)



cleanEx()
nameEx("servo")
### * servo

flush(stderr()); flush(stdout())

### Name: servo
### Title: Servo Data Set
### Aliases: servo
### Keywords: datasets

### ** Examples

data(servo)
lm(class~motor+screw+pgain+vgain, servo)



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
