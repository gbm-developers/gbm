pkgname <- "mboost"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('mboost')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("FP")
### * FP

flush(stderr()); flush(stdout())

### Name: FP
### Title: Fractional Polynomials
### Aliases: FP
### Keywords: datagen

### ** Examples


    data("bodyfat", package = "TH.data")
    tbodyfat <- bodyfat

    ### map covariates into [1, 2]
    indep <- names(tbodyfat)[-2]
    tbodyfat[indep] <- lapply(bodyfat[indep], function(x) {
        x <- x - min(x)
        x / max(x) + 1
    })

    ### generate formula
    fpfm <- as.formula(paste("DEXfat ~ ",
        paste("FP(", indep, ", scaling = FALSE)", collapse = "+")))
    fpfm

    ### fit linear model
    bf_fp <- glmboost(fpfm, data = tbodyfat,
                      control = boost_control(mstop = 3000))

    ### when to stop
    mstop(aic <- AIC(bf_fp))
    plot(aic)

    ### coefficients
    cf <- coef(bf_fp[mstop(aic)])
    length(cf)
    cf[abs(cf) > 0]




cleanEx()
nameEx("Family")
### * Family

flush(stderr()); flush(stdout())

### Name: Family
### Title: Gradient Boosting Families
### Aliases: Family AdaExp Binomial GaussClass GaussReg Gaussian Huber
###   Laplace Poisson GammaReg CoxPH QuantReg ExpectReg NBinomial PropOdds
###   Weibull Loglog Lognormal AUC Gehan Hurdle Multinomial Cindex RCG
### Keywords: models

### ** Examples

### Define a new family
MyGaussian <- function(){
       Family(ngradient = function(y, f, w = 1) y - f,
       loss = function(y, f) (y - f)^2,
       name = "My Gauss Variant")
}
# Now use the new family
data(bodyfat, package = "TH.data")
mod <- mboost(DEXfat ~ ., data = bodyfat, family = MyGaussian())
# N.B. that the family needs to be called with empty brackets



### Multinomial logit model via a linear array model
## One needs to convert the data to a list
myiris <- as.list(iris)
## ... and define a dummy vector with one factor level less
## than the outcome, which is used as reference category.
myiris$class <- factor(levels(iris$Species)[-nlevels(iris$Species)])
## Now fit the linear array model
mlm <- mboost(Species ~ bols(Sepal.Length, df = 2) %O%
                        bols(class, df = 2, contrasts.arg = "contr.dummy"),
              data = myiris,
              family = Multinomial())
coef(mlm) ## one should use more boosting iterations.
head(round(pred <- predict(mlm, type = "response"), 2))

## Prediction with new data:
newdata <- as.list(iris[1,])
## One always needs to keep the dummy vector class as above!
newdata$class <- factor(levels(iris$Species)[-nlevels(iris$Species)])
pred2 <- predict(mlm, type = "response", newdata = newdata)
## check results
pred[1, ]
pred2









cleanEx()
nameEx("baselearners")
### * baselearners

flush(stderr()); flush(stdout())

### Name: baselearners
### Title: Base-learners for Gradient Boosting
### Aliases: baselearners baselearner base-learner bols bbs bspatial brad
###   bkernel brandom btree bmono bmrf buser bns bss %+% %X% %O%
### Keywords: models

### ** Examples


  set.seed(290875)

  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n) + 0.25 * x1
  x3 <- as.factor(sample(0:1, 100, replace = TRUE))
  x4 <- gl(4, 25)
  y <- 3 * sin(x1) + x2^2 + rnorm(n)
  weights <- drop(rmultinom(1, n, rep.int(1, n) / n))

  ### set up base-learners
  spline1 <- bbs(x1, knots = 20, df = 4)
  extract(spline1, "design")[1:10, 1:10]
  extract(spline1, "penalty")
  knots.x2 <- quantile(x2, c(0.25, 0.5, 0.75))
  spline2 <- bbs(x2, knots = knots.x2, df = 5)
  ols3 <- bols(x3)
  extract(ols3)
  ols4 <- bols(x4)

  ### compute base-models
  drop(ols3$dpp(weights)$fit(y)$model) ## same as:
  coef(lm(y ~ x3, weights = weights))

  drop(ols4$dpp(weights)$fit(y)$model) ## same as:
  coef(lm(y ~ x4, weights = weights))

  ### fit model, component-wise
  mod1 <- mboost_fit(list(spline1, spline2, ols3, ols4), y, weights)

  ### more convenient formula interface
  mod2 <- mboost(y ~ bbs(x1, knots = 20, df = 4) +
                     bbs(x2, knots = knots.x2, df = 5) +
                     bols(x3) + bols(x4), weights = weights)
  all.equal(coef(mod1), coef(mod2))


  ### grouped linear effects
  # center x1 and x2 first
  x1 <- scale(x1, center = TRUE, scale = FALSE)
  x2 <- scale(x2, center = TRUE, scale = FALSE)
  model <- gamboost(y ~ bols(x1, x2, intercept = FALSE) +
                        bols(x1, intercept = FALSE) +
                        bols(x2, intercept = FALSE),
                        control = boost_control(mstop = 50))
  coef(model, which = 1)   # one base-learner for x1 and x2
  coef(model, which = 2:3) # two separate base-learners for x1 and x2
                           # zero because they were (not yet) selected.

  ### example for bspatial
  x1 <- runif(250,-pi,pi)
  x2 <- runif(250,-pi,pi)

  y <- sin(x1) * sin(x2) + rnorm(250, sd = 0.4)

  spline3 <- bspatial(x1, x2, knots = 12)
  Xmat <- extract(spline3, "design")
  ## 12 inner knots + 4 boundary knots = 16 knots per direction
  ## THUS: 16 * 16 = 256 columns
  dim(Xmat)
  extract(spline3, "penalty")[1:10, 1:10]

  ## specify number of knots separately
  form1 <- y ~ bspatial(x1, x2, knots = list(x1 = 12, x2 = 14))

  ## decompose spatial effect into parametric part and
  ## deviation with one df
  form2 <- y ~ bols(x1) + bols(x2) + bols(x1, by = x2, intercept = FALSE) +
               bspatial(x1, x2, knots = 12, center = TRUE, df = 1)


  ## specify radial basis function base-learner for spatial effect
  ## and use data-adaptive effective range (theta = NULL, see 'args')
  form3 <- y ~ brad(x1, x2)
  ## Now use different settings, e.g. 50 knots and theta fixed to 0.4
  ## (not really a good setting)
  form4 <- y ~ brad(x1, x2, knots = 50, args = list(theta = 0.4))


  ### random intercept
  id <- factor(rep(1:10, each = 5))
  raneff <- brandom(id)
  extract(raneff, "design")
  extract(raneff, "penalty")

  ## random intercept with non-observed category
  set.seed(1907)
  y <- rnorm(50, mean = rep(rnorm(10), each = 5), sd = 0.1)
  plot(y ~ id)
  # category 10 not observed
  obs <- c(rep(1, 45), rep(0, 5))
  model <- gamboost(y ~ brandom(id), weights = obs)
  coef(model)
  fitted(model)[46:50] # just the grand mean as usual for
                       # random effects models


  ### random slope
  z <- runif(50)
  raneff <- brandom(id, by = z)
  extract(raneff, "design")
  extract(raneff, "penalty")

  ### specify simple interaction model (with main effect)
  n <- 210
  x <- rnorm(n)
  X <- model.matrix(~ x)
  z <- gl(3, n/3)
  Z <- model.matrix(~z)
  beta <- list(c(0,1), c(-3,4), c(2, -4))
  y <- rnorm(length(x), mean = (X * Z[,1]) %*% beta[[1]] +
                               (X * Z[,2]) %*% beta[[2]] +
                               (X * Z[,3]) %*% beta[[3]])
  plot(y ~ x, col = z)
  ## specify main effect and interaction
  mod_glm <- gamboost(y ~ bols(x) + bols(x, by = z),
                  control = boost_control(mstop = 100))
  nd <- data.frame(x, z)
  nd <- nd[order(x),]
  nd$pred_glm <- predict(mod_glm, newdata = nd)
  for (i in seq(along = levels(z)))
      with(nd[nd$z == i,], lines(x, pred_glm, col = z))
  mod_gam <- gamboost(y ~ bbs(x) + bbs(x, by = z, df = 8),
                      control = boost_control(mstop = 100))
  nd$pred_gam <- predict(mod_gam, newdata = nd)
  for (i in seq(along = levels(z)))
      with(nd[nd$z == i,], lines(x, pred_gam, col = z, lty = "dashed"))
  ### convenience function for plotting
  par(mfrow = c(1,3))
  plot(mod_gam)


  ### remove intercept from base-learner
  ### and add explicit intercept to the model
  tmpdata <- data.frame(x = 1:100, y = rnorm(1:100), int = rep(1, 100))
  mod <- gamboost(y ~ bols(int, intercept = FALSE) +
                      bols(x, intercept = FALSE),
                  data = tmpdata,
                  control = boost_control(mstop = 1000))
  cf <- unlist(coef(mod))
  ## add offset
  cf[1] <- cf[1] + mod$offset
  signif(cf, 3)
  signif(coef(lm(y ~ x, data = tmpdata)), 3)

  ### much quicker and better with (mean-) centering
  tmpdata$x_center <- tmpdata$x - mean(tmpdata$x)
  mod_center <- gamboost(y ~ bols(int, intercept = FALSE) +
                             bols(x_center, intercept = FALSE),
                         data = tmpdata,
                         control = boost_control(mstop = 100))
  cf_center <- unlist(coef(mod_center, which=1:2))
  ## due to the shift in x direction we need to subtract
  ## beta_1 * mean(x) to get the correct intercept
  cf_center[1] <- cf_center[1] + mod_center$offset -
                  cf_center[2] * mean(tmpdata$x)
  signif(cf_center, 3)
  signif(coef(lm(y ~ x, data = tmpdata)), 3)


  ### cyclic P-splines
  set.seed(781)
  x <- runif(200, 0,(2*pi))
  y <- rnorm(200, mean=sin(x), sd=0.2)
  newX <- seq(0,2*pi, length=100)
  ### model without cyclic constraints
  mod <- gamboost(y ~ bbs(x, knots = 20))
  ### model with cyclic constraints
  mod_cyclic <- gamboost(y ~ bbs(x, cyclic=TRUE, knots = 20,
                                 boundary.knots=c(0, 2*pi)))
  par(mfrow = c(1,2))
  plot(x,y, main="bbs (non-cyclic)", cex=0.5)
  lines(newX, sin(newX), lty="dotted")
  lines(newX + 2 * pi, sin(newX), lty="dashed")
  lines(newX, predict(mod, data.frame(x = newX)),
        col="red", lwd = 1.5)
  lines(newX + 2 * pi, predict(mod, data.frame(x = newX)),
        col="blue", lwd=1.5)
  plot(x,y, main="bbs (cyclic)", cex=0.5)
  lines(newX, sin(newX), lty="dotted")
  lines(newX + 2 * pi, sin(newX), lty="dashed")
  lines(newX, predict(mod_cyclic, data.frame(x = newX)),
        col="red", lwd = 1.5)
  lines(newX + 2 * pi, predict(mod_cyclic, data.frame(x = newX)),
        col="blue", lwd = 1.5)

  ### use buser() to mimic p-spline base-learner:
  set.seed(1907)
  x <- rnorm(100)
  y <- rnorm(100, mean = x^2, sd = 0.1)
  mod1 <- gamboost(y ~ bbs(x))
  ## now extract design and penalty matrix
  X <- extract(bbs(x), "design")
  K <- extract(bbs(x), "penalty")
  ## use X and K in buser()
  mod2 <- gamboost(y ~ buser(X, K))
  max(abs(predict(mod1) - predict(mod2)))  # same results

  ### use buser() to mimic penalized ordinal base-learner:
  z <- as.ordered(sample(1:3, 100, replace=TRUE))
  y <- rnorm(100, mean = as.numeric(z), sd = 0.1)
  X <- extract(bols(z))
  K <- extract(bols(z), "penalty")
  index <- extract(bols(z), "index")
  mod1 <- gamboost(y ~  buser(X, K, df = 1, index = index))
  mod2 <- gamboost(y ~  bols(z, df = 1))
  max(abs(predict(mod1) - predict(mod2)))  # same results

  ### kronecker product for matrix-valued responses
  data("volcano", package = "datasets")
  layout(matrix(1:2, ncol = 2))

  ## estimate mean of image treating image as matrix
  image(volcano, main = "data")
  x1 <- 1:nrow(volcano)
  x2 <- 1:ncol(volcano)

  vol <- as.vector(volcano)
  mod <- mboost(vol ~ bbs(x1, df = 3, knots = 10)%O%
                      bbs(x2, df = 3, knots = 10),
                      control = boost_control(nu = 0.25))
  mod[250]

  volf <- matrix(fitted(mod), nrow = nrow(volcano))
  image(volf, main = "fitted")


  ### setting contrasts via contrasts.arg
  x <- as.factor(sample(1:4, 100, replace = TRUE))

  ## compute base-learners with different reference categories
  BL1 <- bols(x, contrasts.arg = contr.treatment(4, base = 1)) # default
  BL2 <- bols(x, contrasts.arg = contr.treatment(4, base = 2))
  ## compute 'sum to zero contrasts' using character string
  BL3 <- bols(x, contrasts.arg = "contr.sum")

  ## extract model matrices to check if it works
  extract(BL1)
  extract(BL2)
  extract(BL3)

  ### setting contrasts using named lists in contrasts.arg
  x2 <- as.factor(sample(1:4, 100, replace = TRUE))

  BL4 <- bols(x, x2,
              contrasts.arg = list(x = contr.treatment(4, base = 2),
                                   x2 = "contr.helmert"))
  extract(BL4)

  ### using special contrast: "contr.dummy":
  BL5 <- bols(x, contrasts.arg = "contr.dummy")
  extract(BL5)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("blackboost")
### * blackboost

flush(stderr()); flush(stdout())

### Name: blackboost
### Title: Gradient Boosting with Regression Trees
### Aliases: blackboost
### Keywords: models regression

### ** Examples


### a simple two-dimensional example: cars data
cars.gb <- blackboost(dist ~ speed, data = cars,
                      control = boost_control(mstop = 50))
cars.gb

### plot fit
plot(dist ~ speed, data = cars)
lines(cars$speed, predict(cars.gb), col = "red")





cleanEx()
nameEx("boost_family-class")
### * boost_family-class

flush(stderr()); flush(stdout())

### Name: boost_family-class
### Title: Class "boost\_family": Gradient Boosting Family
### Aliases: boost_family-class show,boost_family-method
### Keywords: classes

### ** Examples


    Laplace()




cleanEx()
nameEx("confint")
### * confint

flush(stderr()); flush(stdout())

### Name: confint.mboost
### Title: Pointwise Bootstrap Confidence Intervals
### Aliases: confint.mboost confint.glmboost plot.mboost.ci lines.mboost.ci
###   print.glmboost.ci
### Keywords: methods

### ** Examples




cleanEx()
nameEx("cvrisk")
### * cvrisk

flush(stderr()); flush(stdout())

### Name: cvrisk
### Title: Cross-Validation
### Aliases: cvrisk cvrisk.mboost print.cvrisk plot.cvrisk cv
### Keywords: models regression

### ** Examples


  data("bodyfat", package = "TH.data")

  ### fit linear model to data
  model <- glmboost(DEXfat ~ ., data = bodyfat, center = TRUE)

  ### AIC-based selection of number of boosting iterations
  maic <- AIC(model)
  maic

  ### inspect coefficient path and AIC-based stopping criterion
  par(mai = par("mai") * c(1, 1, 1, 1.8))
  plot(model)
  abline(v = mstop(maic), col = "lightgray")

  ### 10-fold cross-validation
  cv10f <- cv(model.weights(model), type = "kfold")
  cvm <- cvrisk(model, folds = cv10f, papply = lapply)
  print(cvm)
  mstop(cvm)
  plot(cvm)

  ### 25 bootstrap iterations (manually)
  set.seed(290875)
  n <- nrow(bodyfat)
  bs25 <- rmultinom(25, n, rep(1, n)/n)
  cvm <- cvrisk(model, folds = bs25, papply = lapply)
  print(cvm)
  mstop(cvm)
  plot(cvm)

  ### same by default
  set.seed(290875)
  cvrisk(model, papply = lapply)

  ### 25 bootstrap iterations (using cv)
  set.seed(290875)
  bs25_2 <- cv(model.weights(model), type="bootstrap")
  all(bs25 == bs25_2)



### cvrisk in parallel modes:

## Not run: 
##D ## at least not automatically
##D 
##D ## parallel::mclapply() which is used here for parallelization only runs 
##D ## on unix systems (here we use 2 cores)
##D 
##D     cvrisk(model, mc.cores = 2)
##D 
##D ## infrastructure needs to be set up in advance
##D 
##D     cl <- makeCluster(25) # e.g. to run cvrisk on 25 nodes via PVM
##D     myApply <- function(X, FUN, ...) {
##D       myFun <- function(...) {
##D           library("mboost") # load mboost on nodes
##D           FUN(...)
##D       }
##D       ## further set up steps as required
##D       parLapply(cl = cl, X, myFun, ...)
##D     }
##D     cvrisk(model, papply = myApply)
##D     stopCluster(cl)
## End(Not run)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("gamboost")
### * gamboost

flush(stderr()); flush(stdout())

### Name: mboost
### Title: Gradient Boosting for Additive Models
### Aliases: mboost gamboost
### Keywords: models nonlinear

### ** Examples


    ### a simple two-dimensional example: cars data
    cars.gb <- gamboost(dist ~ speed, data = cars, dfbase = 4,
                        control = boost_control(mstop = 50))
    cars.gb
    AIC(cars.gb, method = "corrected")

    ### plot fit for mstop = 1, ..., 50
    plot(dist ~ speed, data = cars)
    tmp <- sapply(1:mstop(AIC(cars.gb)), function(i)
        lines(cars$speed, predict(cars.gb[i]), col = "red"))
    lines(cars$speed, predict(smooth.spline(cars$speed, cars$dist),
                              cars$speed)$y, col = "green")

    ### artificial example: sinus transformation
    x <- sort(runif(100)) * 10
    y <- sin(x) + rnorm(length(x), sd = 0.25)
    plot(x, y)
    ### linear model
    lines(x, fitted(lm(y ~ sin(x) - 1)), col = "red")
    ### GAM
    lines(x, fitted(gamboost(y ~ x,
                    control = boost_control(mstop = 500))),
          col = "green")




cleanEx()
nameEx("glmboost")
### * glmboost

flush(stderr()); flush(stdout())

### Name: glmboost
### Title: Gradient Boosting with Component-wise Linear Models
### Aliases: glmboost glmboost.formula glmboost.matrix glmboost.default
### Keywords: models regression

### ** Examples


    ### a simple two-dimensional example: cars data
    cars.gb <- glmboost(dist ~ speed, data = cars,
                        control = boost_control(mstop = 2000),
                        center = FALSE)
    cars.gb

    ### coefficients should coincide
    cf <- coef(cars.gb, off2int = TRUE)     ## add offset to intercept
    coef(cars.gb) + c(cars.gb$offset, 0)    ## add offset to intercept (by hand)
    signif(cf, 3)
    signif(coef(lm(dist ~ speed, data = cars)), 3)
    ## almost converged. With higher mstop the results get even better

    ### now we center the design matrix for
    ### much quicker "convergence"
    cars.gb_centered <- glmboost(dist ~ speed, data = cars,
                                 control = boost_control(mstop = 2000),
                                 center = TRUE)

    ## plot coefficient paths of glmboost
    par(mfrow=c(1,2), mai = par("mai") * c(1, 1, 1, 2.5))
    plot(cars.gb, main = "without centering")
    plot(cars.gb_centered, main = "with centering")

    ### alternative loss function: absolute loss
    cars.gbl <- glmboost(dist ~ speed, data = cars,
                         control = boost_control(mstop = 1000),
                         family = Laplace())
    cars.gbl
    coef(cars.gbl, off2int = TRUE)

    ### plot fit
    par(mfrow = c(1,1))
    plot(dist ~ speed, data = cars)
    lines(cars$speed, predict(cars.gb), col = "red")     ## quadratic loss
    lines(cars$speed, predict(cars.gbl), col = "green")  ## absolute loss

    ### Huber loss with adaptive choice of delta
    cars.gbh <- glmboost(dist ~ speed, data = cars,
                         control = boost_control(mstop = 1000),
                         family = Huber())

    lines(cars$speed, predict(cars.gbh), col = "blue")   ## Huber loss
    legend("topleft", col = c("red", "green", "blue"), lty = 1,
           legend = c("Gaussian", "Laplace", "Huber"), bty = "n")

    ### sparse high-dimensional example that makes use of the matrix
    ### interface of glmboost and uses the matrix representation from
    ### package Matrix
    library("Matrix")
    n <- 100
    p <- 10000
    ptrue <- 10
    X <- Matrix(0, nrow = n, ncol = p)
    X[sample(1:(n * p), floor(n * p / 20))] <- runif(floor(n * p / 20))
    beta <- numeric(p)
    beta[sample(1:p, ptrue)] <- 10
    y <- drop(X %*% beta + rnorm(n, sd = 0.1))
    mod <- glmboost(y = y, x = X, center = TRUE) ### mstop needs tuning
    coef(mod, which = which(beta > 0))




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("mboost_fit")
### * mboost_fit

flush(stderr()); flush(stdout())

### Name: mboost_fit
### Title: Model-based Gradient Boosting
### Aliases: mboost_fit
### Keywords: models nonlinear

### ** Examples

  data("bodyfat", package = "TH.data")

  ### formula interface: additive Gaussian model with
  ### a non-linear step-function in `age', a linear function in `waistcirc'
  ### and a smooth non-linear smooth function in `hipcirc'
  mod <- mboost(DEXfat ~ btree(age) + bols(waistcirc) + bbs(hipcirc),
                data = bodyfat)
  layout(matrix(1:6, nc = 3, byrow = TRUE))
  plot(mod, main = "formula")

  ### the same
  with(bodyfat,
       mod <- mboost_fit(list(btree(age), bols(waistcirc), bbs(hipcirc)),
                         response = DEXfat))
  plot(mod, main = "base-learner")



cleanEx()
nameEx("mboost_package")
### * mboost_package

flush(stderr()); flush(stdout())

### Name: mboost-package
### Title: mboost: Model-Based Boosting
### Aliases: mboost-package mboost_package package_mboost package-mboost
### Keywords: package smooth nonparametric models

### ** Examples




cleanEx()
nameEx("methods")
### * methods

flush(stderr()); flush(stdout())

### Name: methods
### Title: Methods for Gradient Boosting Objects
### Aliases: mboost_methods print.glmboost print.mboost summary.mboost
###   coef.mboost coef.glmboost [.mboost AIC.mboost mstop mstop.gbAIC
###   mstop.mboost mstop.cvrisk mstop<- predict.mboost predict.gamboost
###   predict.blackboost predict.glmboost fitted.mboost residuals.mboost
###   resid.mboost variable.names.glmboost variable.names.mboost risk
###   risk.mboost extract extract.mboost extract.gamboost extract.glmboost
###   extract.blackboost extract.blg extract.bl_lin extract.bl_tree
###   logLik.mboost hatvalues.gamboost hatvalues.glmboost selected
###   selected.mboost nuisance nuisance.mboost downstream.test
### Keywords: methods

### ** Examples


  ### a simple two-dimensional example: cars data
  cars.gb <- glmboost(dist ~ speed, data = cars,
                      control = boost_control(mstop = 2000),
                      center = FALSE)
  cars.gb

  ### initial number of boosting iterations
  mstop(cars.gb)

  ### AIC criterion
  aic <- AIC(cars.gb, method = "corrected")
  aic

  ### extract coefficients for glmboost
  coef(cars.gb)
  coef(cars.gb, off2int = TRUE)        # offset added to intercept
  coef(lm(dist ~ speed, data = cars))  # directly comparable

  cars.gb_centered <- glmboost(dist ~ speed, data = cars,
                               center = TRUE)
  selected(cars.gb_centered)           # intercept never selected
  coef(cars.gb_centered)               # intercept implicitly estimated
                                       # and thus returned
  ## intercept is internally corrected for mean-centering
  - mean(cars$speed) * coef(cars.gb_centered, which="speed") # = intercept
  # not asked for intercept thus not returned
  coef(cars.gb_centered, which="speed")
  # explicitly asked for intercept
  coef(cars.gb_centered, which=c("Intercept", "speed"))

  ### enhance or restrict model
  cars.gb <- gamboost(dist ~ speed, data = cars,
                      control = boost_control(mstop = 100, trace = TRUE))
  cars.gb[10]
  cars.gb[100, return = FALSE] # no refitting required
  cars.gb[150, return = FALSE] # only iterations 101 to 150
                               # are newly fitted

  ### coefficients for optimal number of boosting iterations
  coef(cars.gb[mstop(aic)])
  plot(cars$dist, predict(cars.gb[mstop(aic)]),
       ylim = range(cars$dist))
  abline(a = 0, b = 1)

  ### example for extraction of coefficients
  set.seed(1907)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  x4 <- rnorm(n)
  int <- rep(1, n)
  y <- 3 * x1^2 - 0.5 * x2 + rnorm(n, sd = 0.1)
  data <- data.frame(y = y, int = int, x1 = x1, x2 = x2, x3 = x3, x4 = x4)

  model <- gamboost(y ~ bols(int, intercept = FALSE) +
                        bbs(x1, center = TRUE, df = 1) +
                        bols(x1, intercept = FALSE) +
                        bols(x2, intercept = FALSE) +
                        bols(x3, intercept = FALSE) +
                        bols(x4, intercept = FALSE),
                    data = data, control = boost_control(mstop = 500))

  coef(model) # standard output (only selected base-learners)
  coef(model,
       which = 1:length(variable.names(model))) # all base-learners
  coef(model, which = "x1") # shows all base-learners for x1

  cf1 <- coef(model, which = c(1,3,4), aggregate = "cumsum")
  tmp <- sapply(cf1, function(x) x)
  matplot(tmp, type = "l", main = "Coefficient Paths")

  cf1_all <- coef(model, aggregate = "cumsum")
  cf1_all <- lapply(cf1_all, function(x) x[, ncol(x)]) # last element
  ## same as coef(model)

  cf2 <- coef(model, aggregate = "none")
  cf2 <- lapply(cf2, rowSums) # same as coef(model)

  ### example continued for extraction of predictions

  yhat <- predict(model) # standard prediction; here same as fitted(model)
  p1 <- predict(model, which = "x1") # marginal effects of x1
  orderX <- order(data$x1)
  ## rowSums needed as p1 is a matrix
  plot(data$x1[orderX], rowSums(p1)[orderX], type = "b")

  ## better: predictions on a equidistant grid
  new_data <- data.frame(x1 = seq(min(data$x1), max(data$x1), length = 100))
  p2 <- predict(model, newdata = new_data, which = "x1")
  lines(new_data$x1, rowSums(p2), col = "red")

  ### extraction of model characteristics
  extract(model, which = "x1")  # design matrices for x1
  extract(model, what = "penalty", which = "x1") # penalty matrices for x1
  extract(model, what = "lambda", which = "x1") # df and corresponding lambda for x1
       ## note that bols(x1, intercept = FALSE) is unpenalized

  extract(model, what = "bnames")  ## name of complete base-learner
  extract(model, what = "variable.names") ## only variable names
  variable.names(model)            ## the same

  ### extract from base-learners
  extract(bbs(x1), what = "design")
  extract(bbs(x1), what = "penalty")
  ## weights and lambda can only be extracted after using dpp
  weights <- rep(1, length(x1))
  extract(bbs(x1)$dpp(weights), what = "lambda")



cleanEx()
nameEx("plot")
### * plot

flush(stderr()); flush(stdout())

### Name: plot
### Title: Plot effect estimates of boosting models
### Aliases: plot plot.glmboost plot.mboost lines.mboost
### Keywords: methods

### ** Examples


### a simple example: cars data with one random variable
set.seed(1234)
cars$z <- rnorm(50)

########################################
## Plot linear models
########################################

## fit a linear model
cars.lm <- glmboost(dist ~ speed + z, data = cars)

## plot coefficient paths of glmboost
par(mfrow = c(3, 1), mar = c(4, 4, 4, 8))
plot(cars.lm,
     main = "Coefficient paths (offset not included)")
plot(cars.lm, off2int = TRUE,
     main = "Coefficient paths (offset included in intercept)")

## plot coefficient paths only for the first 15 steps,
## i.e., bevore z is selected
mstop(cars.lm) <- 15
plot(cars.lm, off2int = TRUE, main = "z is not yet selected")


########################################
## Plot additive models; basics
########################################

## fit an additive model
cars.gam <- gamboost(dist ~ speed + z, data = cars)

## plot effects
par(mfrow = c(1, 2), mar = c(4, 4, 0.1, 0.1))
plot(cars.gam)

## use same y-lims
plot(cars.gam, ylim = c(-50, 50))

## plot only the effect of speed
plot(cars.gam, which = "speed")
## as partial matching is used we could also use
plot(cars.gam, which = "sp")


########################################
## More complex plots
########################################

## Let us use more boosting iterations and compare the effects.

## We change the plot type and plot both effects in one figure:
par(mfrow = c(1, 1), mar = c(4, 4, 4, 0.1))
mstop(cars.gam) <- 100
plot(cars.gam, which = 1, col = "red", type = "l", rug = FALSE,
     main = "Compare effect for various models")

## Now the same model with 1000 iterations
mstop(cars.gam) <- 1000
lines(cars.gam, which = 1, col = "grey", lty = "dotted")

## There are some gaps in the data. Use newdata to get a smoother curve:
newdata <- data.frame(speed = seq(min(cars$speed), max(cars$speed),
                                  length = 200))
lines(cars.gam, which = 1, col = "grey", lty = "dashed",
      newdata = newdata)

## The model with 1000 steps seems to overfit the data.
## Usually one should use e.g. cross-validation to tune the model.

## Finally we refit the model using linear effects as comparison
cars.glm <- gamboost(dist ~ speed + z, baselearner = bols, data = cars)
lines(cars.glm, which = 1, col = "black")
## We see that all effects are more or less linear.

## Add a legend
legend("topleft", title = "Model",
       legend = c("... with mstop = 100", "... with mstop = 1000",
         "... with linear effects"),
       lty = c("solid", "dashed", "solid"),
       col = c("red", "grey", "black"))




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("stabsel")
### * stabsel

flush(stderr()); flush(stdout())

### Name: stabsel
### Title: Stability Selection
### Aliases: stabsel stabsel.mboost stabsel_parameters.mboost
### Keywords: nonparametric

### ** Examples

  ## make data set available
  data("bodyfat", package = "TH.data")
  ## set seed
  set.seed(1234)

  ### low-dimensional example
  mod <- glmboost(DEXfat ~ ., data = bodyfat)

  ## compute cutoff ahead of running stabsel to see if it is a sensible
  ## parameter choice.
  ##   p = ncol(bodyfat) - 1 (= Outcome) + 1 ( = Intercept)
  stabsel_parameters(q = 3, PFER = 1, p = ncol(bodyfat) - 1 + 1,
                     sampling.type = "MB")
  ## the same:
  stabsel(mod, q = 3, PFER = 1, sampling.type = "MB", eval = FALSE)




cleanEx()
nameEx("survFit")
### * survFit

flush(stderr()); flush(stdout())

### Name: survFit
### Title: Survival Curves for a Cox Proportional Hazards Model
### Aliases: survFit survFit.mboost plot.survFit

### ** Examples


library("survival")
data("ovarian", package = "survival")

fm <- Surv(futime,fustat) ~ age + resid.ds + rx + ecog.ps
fit <- glmboost(fm, data = ovarian, family = CoxPH(),
    control=boost_control(mstop = 500))

S1 <- survFit(fit)
S1
newdata <- ovarian[c(1,3,12),]
S2 <- survFit(fit, newdata = newdata)
S2

plot(S1)



cleanEx()
nameEx("varimp")
### * varimp

flush(stderr()); flush(stdout())

### Name: varimp
### Title: Variable Importance
### Aliases: varimp varimp.mboost plot.varimp as.data.frame.varimp

### ** Examples


data(iris)
### glmboost with multiple variables and intercept
iris$setosa <- factor(iris$Species == "setosa")
iris_glm <- glmboost(setosa ~ 1 + Sepal.Width + Sepal.Length + Petal.Width +
                         Petal.Length,
                     data = iris, control = boost_control(mstop = 50), 
                     family = Binomial(link = c("logit")))
varimp(iris_glm)
### importance plot with four bars only
plot(varimp(iris_glm), nbars = 4)

### gamboost with multiple variables
iris_gam <- gamboost(Sepal.Width ~ 
                         bols(Sepal.Length, by = setosa) +
                         bbs(Sepal.Length, by = setosa, center = TRUE) + 
                         bols(Petal.Width) +
                         bbs(Petal.Width, center = TRUE) + 
                         bols(Petal.Length) +
                         bbs(Petal.Length, center = TRUE),
                     data = iris)
varimp(iris_gam)
### stacked importance plot with base-learners in rev. alphabetical order
plot(varimp(iris_gam), blorder = "rev_alphabetical")

### similar ggplot
## Not run: 
##D  
##D library(ggplot2)
##D ggplot(data.frame(varimp(iris_gam)), aes(variable, reduction, fill = blearner)) + 
##D     geom_bar(stat = "identity") + coord_flip() 
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
