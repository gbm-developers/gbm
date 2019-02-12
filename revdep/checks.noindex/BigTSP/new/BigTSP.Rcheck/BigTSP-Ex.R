pkgname <- "BigTSP"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('BigTSP')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("LDCA")
### * LDCA

flush(stderr()); flush(stdout())

### Name: LDCA
### Title: Linear Discriminant Analysis based on Top Scoring Pair
### Aliases: LDCA
### Keywords: ~kwd1 ~kwd2

### ** Examples

library(glmnet)
x=matrix(rnorm(100*20),100,20)
y=rbinom(100,1,0.5)
fit=LDCA(x,y)
print(fit)
predict(fit,newx=x[1:10,]) # make predictions



cleanEx()
nameEx("cv.LDCA")
### * cv.LDCA

flush(stderr()); flush(stdout())

### Name: cv.LDCA
### Title: Cross validation for LDCA
### Aliases: cv.LDCA
### Keywords: ~kwd1 ~kwd2

### ** Examples

library(glmnet)
x=matrix(rnorm(50*20),50,20)
y=rbinom(50,1,0.5)
cvfit=cv.LDCA(x,y,nfolds=5)
predict(cvfit,x[1:10,],s="lambda.min")



cleanEx()
nameEx("predict.LDCA")
### * predict.LDCA

flush(stderr()); flush(stdout())

### Name: predict.LDCA
### Title: predict function for LDCA
### Aliases: predict.LDCA
### Keywords: ~kwd1 ~kwd2

### ** Examples

library(glmnet)
x=matrix(rnorm(50*20),50,20)
y=rbinom(50,1,0.5)
cvfit=cv.LDCA(x,y,nfolds=5)
predict(cvfit,x[1:10,],s="lambda.min")



cleanEx()
nameEx("predict.cv.LDCA")
### * predict.cv.LDCA

flush(stderr()); flush(stdout())

### Name: predict.cv.LDCA
### Title: prediction function for cv.LDCA
### Aliases: predict.cv.LDCA
### Keywords: ~kwd1 ~kwd2

### ** Examples

library(glmnet)
x=matrix(rnorm(50*20),50,20)
y=rbinom(50,1,0.5)
cvfit=cv.LDCA(x,y,nfolds=5)
predict(cvfit,x[1:10,],s="lambda.min")



cleanEx()
nameEx("predict.tsp.gbm")
### * predict.tsp.gbm

flush(stderr()); flush(stdout())

### Name: predict.tsp.gbm
### Title: prediction function for tsp.gbm
### Aliases: predict.tsp.gbm
### Keywords: ~kwd1 ~kwd2

### ** Examples

library(gbm)
x=matrix(rnorm(100*20),100,20)
y=rbinom(100,1,0.5)
fit=tsp.gbm(x,y)
predict(fit,x[1:10,],n.trees=5)



cleanEx()
nameEx("predict.tsp.randomForest")
### * predict.tsp.randomForest

flush(stderr()); flush(stdout())

### Name: predict.tsp.randomForest
### Title: prediction function for tsp.randomForest
### Aliases: predict.tsp.randomForest
### Keywords: ~kwd1 ~kwd2

### ** Examples

library(randomForest)
x=matrix(rnorm(100*20),100,20)
y=rbinom(100,1,0.5)
y=as.factor(y)
fit=tsp.randomForest(x,y)
predict(fit,x[1:10,])



cleanEx()
nameEx("predict.tsp.tree")
### * predict.tsp.tree

flush(stderr()); flush(stdout())

### Name: predict.tsp.tree
### Title: prediction function for tsp.tree
### Aliases: predict.tsp.tree
### Keywords: ~kwd1 ~kwd2

### ** Examples

library(tree)
x=matrix(rnorm(100*20),100,20)
y=rbinom(100,1,0.5)
y=as.factor(y)
data=data.frame(y,x)
tr=tsp.tree(x,y)
predict(tr,data[1:10,])



cleanEx()
nameEx("print.LDCA")
### * print.LDCA

flush(stderr()); flush(stdout())

### Name: print.LDCA
### Title: print the LDCA object
### Aliases: print.LDCA
### Keywords: ~kwd1 ~kwd2

### ** Examples

library(glmnet)
x=matrix(rnorm(100*20),100,20)
y=rbinom(100,1,0.5)
fit=LDCA(x,y)
print(fit)



cleanEx()
nameEx("print.cv.LDCA")
### * print.cv.LDCA

flush(stderr()); flush(stdout())

### Name: print.cv.LDCA
### Title: print function for cv.LDCA
### Aliases: print.cv.LDCA
### Keywords: ~kwd1 ~kwd2

### ** Examples

library(glmnet)
x=matrix(rnorm(50*20),50,20)
y=rbinom(50,1,0.5)
cvfit=cv.LDCA(x,y,nfolds=5)
print(cvfit)



cleanEx()
nameEx("tsp.gbm")
### * tsp.gbm

flush(stderr()); flush(stdout())

### Name: tsp.gbm
### Title: Fits generalized boosted logistic regression models based on Top
###   Scoring Pairs.
### Aliases: tsp.gbm
### Keywords: ~kwd1 ~kwd2

### ** Examples

library(gbm)
x=matrix(rnorm(100*20),100,20)
y=rbinom(100,1,0.5)
fit=tsp.gbm(x,y)
predict(fit,x[1:10,],n.trees=5)



cleanEx()
nameEx("tsp.randomForest")
### * tsp.randomForest

flush(stderr()); flush(stdout())

### Name: tsp.randomForest
### Title: Classification with Random Forest based on Top Scoring Pairs
### Aliases: tsp.randomForest
### Keywords: ~kwd1 ~kwd2

### ** Examples

library(randomForest)
x=matrix(rnorm(100*20),100,20)
y=rbinom(100,1,0.5)
y=as.factor(y)
fit=tsp.randomForest(x,y)
predict(fit,x[1:10,])
plot(fit)



cleanEx()
nameEx("tsp.tree")
### * tsp.tree

flush(stderr()); flush(stdout())

### Name: tsp.tree
### Title: Fit a Classification Tree based on Top Scoring Pairs.
### Aliases: tsp.tree

### ** Examples

library(tree)
x=matrix(rnorm(100*20),100,20)
y=rbinom(100,1,0.5)
y=as.factor(y)
data=data.frame(y,x)
tr=tsp.tree(x,y)
predict(tr,data[1:10,])
plot(tr)
text(tr)



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
