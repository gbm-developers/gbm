pkgname <- "bst"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('bst')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("bst")
### * bst

flush(stderr()); flush(stdout())

### Name: bst
### Title: Boosting for Classification and Regression
### Aliases: bst print.bst predict.bst plot.bst coef.bst fpartial.bst
### Keywords: classification

### ** Examples

x <- matrix(rnorm(100*5),ncol=5)
c <- 2*x[,1]
p <- exp(c)/(exp(c)+exp(-c))
y <- rbinom(100,1,p)
y[y != 1] <- -1
x <- as.data.frame(x)
dat.m <- bst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls")
predict(dat.m)
dat.m1 <- bst(x, y, ctrl = bst_control(twinboost=TRUE, 
coefir=coef(dat.m), xselect.init = dat.m$xselect, mstop=50))
dat.m2 <- rbst(x, y, ctrl = bst_control(mstop=50, s=0, trace=TRUE), 
rfamily = "thinge", learner = "ls")
predict(dat.m2)



cleanEx()
nameEx("bst.sel")
### * bst.sel

flush(stderr()); flush(stdout())

### Name: bst.sel
### Title: Function to select number of predictors
### Aliases: bst.sel
### Keywords: models regression

### ** Examples

## Not run: 
##D x <- matrix(rnorm(100*100), nrow = 100, ncol = 100)
##D y <- x[,1] * 2 + x[,2] * 2.5 + rnorm(100)
##D sel <- bst.sel(x, y, q=10)
##D library("hdi")
##D fit.multi <- hdi(x, y, method = "multi.split",
##D model.selector =bst.sel,
##D args.model.selector=list(type="firstq", q=10))
##D fit.multi
##D fit.multi$pval[1:10] ## the first 10 p-values
##D fit.multi <- hdi(x, y, method = "multi.split",
##D model.selector =bst.sel,
##D args.model.selector=list(type="cv"))
##D fit.multi
##D fit.multi$pval[1:10] ## the first 10 p-values
## End(Not run)



cleanEx()
nameEx("cv.bst")
### * cv.bst

flush(stderr()); flush(stdout())

### Name: cv.bst
### Title: Cross-Validation for Boosting
### Aliases: cv.bst

### ** Examples

## Not run: 
##D x <- matrix(rnorm(100*5),ncol=5)
##D c <- 2*x[,1]
##D p <- exp(c)/(exp(c)+exp(-c))
##D y <- rbinom(100,1,p)
##D y[y != 1] <- -1
##D x <- as.data.frame(x)
##D cv.bst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls", type="loss")
##D cv.bst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls", type="error")
##D dat.m <- bst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls")
##D dat.m1 <- cv.bst(x, y, ctrl = bst_control(twinboost=TRUE, coefir=coef(dat.m), 
##D xselect.init = dat.m$xselect, mstop=50), family = "hinge", learner="ls")
## End(Not run)



cleanEx()
nameEx("cv.rbst")
### * cv.rbst

flush(stderr()); flush(stdout())

### Name: cv.rbst
### Title: Cross-Validation for Nonconvex Loss Boosting
### Aliases: cv.rbst

### ** Examples

## Not run: 
##D x <- matrix(rnorm(100*5),ncol=5)
##D c <- 2*x[,1]
##D p <- exp(c)/(exp(c)+exp(-c))
##D y <- rbinom(100,1,p)
##D y[y != 1] <- -1
##D x <- as.data.frame(x)
##D cv.rbst(x, y, ctrl = bst_control(mstop=50), rfamily = "thinge", learner = "ls", type="lose")
##D cv.rbst(x, y, ctrl = bst_control(mstop=50), rfamily = "thinge", learner = "ls", type="error")
##D dat.m <- rbst(x, y, ctrl = bst_control(mstop=50), rfamily = "thinge", learner = "ls")
##D dat.m1 <- cv.rbst(x, y, ctrl = bst_control(twinboost=TRUE, coefir=coef(dat.m), 
##D xselect.init = dat.m$xselect, mstop=50), family = "thinge", learner="ls")
## End(Not run)



cleanEx()
nameEx("ex1data")
### * ex1data

flush(stderr()); flush(stdout())

### Name: ex1data
### Title: Generating Three-class Data with 50 Predictors
### Aliases: ex1data
### Keywords: classification

### ** Examples

## Not run: 
##D dat <- ex1data(100, p=5)
##D mhingebst(x=dat$x, y=dat$y)
## End(Not run)



cleanEx()
nameEx("mada")
### * mada

flush(stderr()); flush(stdout())

### Name: mada
### Title: Multi-class AdaBoost
### Aliases: mada
### Keywords: classification

### ** Examples

data(iris)
mada(xtr=iris[,-5], ytr=iris[,5])



cleanEx()
nameEx("mbst")
### * mbst

flush(stderr()); flush(stdout())

### Name: mbst
### Title: Boosting for Multi-Classification
### Aliases: mbst print.mbst predict.mbst fpartial.mbst
### Keywords: classification

### ** Examples

x <- matrix(rnorm(100*5),ncol=5)
c <- quantile(x[,1], prob=c(0.33, 0.67))
y <- rep(1, 100)
y[x[,1] > c[1] & x[,1] < c[2] ] <- 2
y[x[,1] > c[2]] <- 3
x <- as.data.frame(x)
dat.m <- mbst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls")
predict(dat.m)
dat.m1 <- mbst(x, y, ctrl = bst_control(twinboost=TRUE, 
f.init=predict(dat.m), xselect.init = dat.m$xselect, mstop=50))
dat.m2 <- rmbst(x, y, ctrl = bst_control(mstop=50, s=1, trace=TRUE), 
rfamily = "thinge", learner = "ls")
predict(dat.m2)



cleanEx()
nameEx("mhingebst")
### * mhingebst

flush(stderr()); flush(stdout())

### Name: mhingebst
### Title: Boosting for Multi-class Classification
### Aliases: mhingebst print.mhingebst predict.mhingebst fpartial.mhingebst
### Keywords: classification

### ** Examples

## Not run: 
##D dat <- ex1data(100, p=5)
##D res <- mhingebst(x=dat$x, y=dat$y)
## End(Not run)



cleanEx()
nameEx("mhingeova")
### * mhingeova

flush(stderr()); flush(stdout())

### Name: mhingeova
### Title: Multi-class HingeBoost
### Aliases: mhingeova print.mhingeova
### Keywords: classification

### ** Examples

## Not run: 
##D dat1 <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/
##D thyroid-disease/ann-train.data")
##D dat2 <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/
##D thyroid-disease/ann-test.data")
##D res <- mhingeova(xtr=dat1[,-22], ytr=dat1[,22], xte=dat2[,-22], yte=dat2[,22], 
##D cost=c(2/3, 0.5, 0.5), nu=0.5, learner="ls", m1=100, K=5, cv1=FALSE, 
##D twinboost=TRUE, m2= 200, cv2=FALSE)
##D res <- mhingeova(xtr=dat1[,-22], ytr=dat1[,22], xte=dat2[,-22], yte=dat2[,22], 
##D cost=c(2/3, 0.5, 0.5), nu=0.5, learner="ls", m1=100, K=5, cv1=FALSE, 
##D twinboost=TRUE, m2= 200, cv2=TRUE)
## End(Not run)



cleanEx()
nameEx("rbst")
### * rbst

flush(stderr()); flush(stdout())

### Name: rbst
### Title: Robust Boosting for Robust Loss Functions
### Aliases: rbst
### Keywords: classification

### ** Examples

x <- matrix(rnorm(100*5),ncol=5)
c <- 2*x[,1]
p <- exp(c)/(exp(c)+exp(-c))
y <- rbinom(100,1,p)
y[y != 1] <- -1
y[1:10] <- -y[1:10]
x <- as.data.frame(x)
dat.m <- bst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls")
predict(dat.m)
dat.m1 <- bst(x, y, ctrl = bst_control(twinboost=TRUE, 
coefir=coef(dat.m), xselect.init = dat.m$xselect, mstop=50))
dat.m2 <- rbst(x, y, ctrl = bst_control(mstop=50, s=0, trace=TRUE), 
rfamily = "thinge", learner = "ls")
predict(dat.m2)



cleanEx()
nameEx("rbstpath")
### * rbstpath

flush(stderr()); flush(stdout())

### Name: rbstpath
### Title: Robust Boosting Path for Nonconvex Loss Functions
### Aliases: rbstpath
### Keywords: classification

### ** Examples

x <- matrix(rnorm(100*5),ncol=5)
c <- 2*x[,1]
p <- exp(c)/(exp(c)+exp(-c))
y <- rbinom(100,1,p)
y[y != 1] <- -1
y[1:10] <- -y[1:10]
x <- as.data.frame(x)
dat.m <- bst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls")
predict(dat.m)
dat.m1 <- bst(x, y, ctrl = bst_control(twinboost=TRUE, 
coefir=coef(dat.m), xselect.init = dat.m$xselect, mstop=50))
dat.m2 <- rbst(x, y, ctrl = bst_control(mstop=50, s=0, trace=TRUE), 
rfamily = "thinge", learner = "ls")
predict(dat.m2)
rmstop <- seq(10, 40, by=10)
dat.m3 <- rbstpath(x, y, rmstop, ctrl=bst_control(s=0), rfamily = "thinge", learner = "ls")



cleanEx()
nameEx("rmbst")
### * rmbst

flush(stderr()); flush(stdout())

### Name: rmbst
### Title: Robust Boosting for Multi-class Robust Loss Functions
### Aliases: rmbst
### Keywords: classification

### ** Examples

x <- matrix(rnorm(100*5),ncol=5)
c <- quantile(x[,1], prob=c(0.33, 0.67))
y <- rep(1, 100)
y[x[,1] > c[1] & x[,1] < c[2] ] <- 2
y[x[,1] > c[2]] <- 3
x <- as.data.frame(x)
x <- as.data.frame(x)
dat.m <- mbst(x, y, ctrl = bst_control(mstop=50), family = "hinge", learner = "ls")
predict(dat.m)
dat.m1 <- mbst(x, y, ctrl = bst_control(twinboost=TRUE, 
f.init=predict(dat.m), xselect.init = dat.m$xselect, mstop=50))
dat.m2 <- rmbst(x, y, ctrl = bst_control(mstop=50, s=1, trace=TRUE), 
rfamily = "thinge", learner = "ls")
predict(dat.m2)



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
