pkgname <- "personalized"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('personalized')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("LaLonde")
### * LaLonde

flush(stderr()); flush(stdout())

### Name: LaLonde
### Title: National Supported Work Study Data
### Aliases: LaLonde
### Keywords: datasets

### ** Examples

data(LaLonde)
y <- LaLonde$outcome

trt <- LaLonde$treat

x.varnames <- c("age", "educ", "black", "hisp", "white",
                "marr", "nodegr", "log.re75", "u75")

# covariates
data.x <- LaLonde[, x.varnames]

# construct design matrix (with no intercept)
x <- model.matrix(~ -1 + ., data = data.x)

const.propens <- function(x, trt)
{
    mean.trt <- mean(trt == "Trt")
    rep(mean.trt, length(trt))
}

subgrp_fit_w <- fit.subgroup(x = x, y = y, trt = trt,
    loss = "logistic_loss_lasso",
    propensity.func = const.propens,
    cutpoint = 0,
    type.measure = "auc",
    nfolds = 10)

summary(subgrp_fit_w)



cleanEx()
nameEx("check.overlap")
### * check.overlap

flush(stderr()); flush(stdout())

### Name: check.overlap
### Title: Check propensity score overlap
### Aliases: check.overlap

### ** Examples

library(personalized)

set.seed(123)
n.obs  <- 250
n.vars <- 15
x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


# simulate non-randomized treatment
xbetat   <- 0.25 + 0.5 * x[,11] - 0.5 * x[,12]
trt.prob <- exp(xbetat) / (1 + exp(xbetat))
trt01    <- rbinom(n.obs, 1, prob = trt.prob)

# create function for fitting propensity score model
prop.func <- function(x, trt)
{
    # fit propensity score model
    propens.model <- cv.glmnet(y = trt,
                               x = x, family = "binomial")
    pi.x <- predict(propens.model, s = "lambda.min",
                    newx = x, type = "response")[,1]
    pi.x
}

check.overlap(x = x,
              trt = trt01,
              propensity.func = prop.func)

# now add density plot with histogram
check.overlap(x = x,
              trt = trt01,
              type = "both",
              propensity.func = prop.func)


# simulated non-randomized treatment with multiple levels
xbetat_1   <- 0.15 + 0.5 * x[,9] - 0.25 * x[,12]
xbetat_2   <- 0.15 - 0.5 * x[,11] + 0.25 * x[,15]
trt.1.prob <- exp(xbetat_1) / (1 + exp(xbetat_1) + exp(xbetat_2))
trt.2.prob <- exp(xbetat_2) / (1 + exp(xbetat_1) + exp(xbetat_2))
trt.3.prob <- 1 - (trt.1.prob + trt.2.prob)
prob.mat <- cbind(trt.1.prob, trt.2.prob, trt.3.prob)
trt    <- apply(prob.mat, 1, function(rr) rmultinom(1, 1, prob = rr))
trt    <- apply(trt, 2, function(rr) which(rr == 1))

# use multinomial logistic regression model with lasso penalty for propensity
propensity.multinom.lasso <- function(x, trt)
{
    if (!is.factor(trt)) trt <- as.factor(trt)
    gfit <- cv.glmnet(y = trt, x = x, family = "multinomial")

    # predict returns a matrix of probabilities:
    # one column for each treatment level
    propens <- drop(predict(gfit, newx = x, type = "response", s = "lambda.min",
                            nfolds = 5, alpha = 0))

    # return the probability corresponding to the
    # treatment that was observed
    probs <- propens[,match(levels(trt), colnames(propens))]

    probs
}

check.overlap(x = x,
              trt = trt,
              type = "histogram",
              propensity.func = propensity.multinom.lasso)






cleanEx()
nameEx("create.augmentation.function")
### * create.augmentation.function

flush(stderr()); flush(stdout())

### Name: create.augmentation.function
### Title: Creation of augmentation functions
### Aliases: create.augmentation.function

### ** Examples

library(personalized)

set.seed(123)
n.obs  <- 500
n.vars <- 15
x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


# simulate non-randomized treatment
xbetat   <- 0.5 + 0.5 * x[,7] - 0.5 * x[,9]
trt.prob <- exp(xbetat) / (1 + exp(xbetat))
trt01    <- rbinom(n.obs, 1, prob = trt.prob)

trt      <- 2 * trt01 - 1

# simulate response
# delta below drives treatment effect heterogeneity
delta <- 2 * (0.5 + x[,2] - x[,3] - x[,11] + x[,1] * x[,12] )
xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13] + 0.5 * x[,15] ^ 2
xbeta <- xbeta + delta * trt

# continuous outcomes
y <- drop(xbeta) + rnorm(n.obs, sd = 2)

aug.func <- create.augmentation.function(family = "gaussian",
                                         crossfit = TRUE,
                                         nfolds.crossfit = 10,
                                         cv.glmnet.args = list(type.measure = "mae",
                                                               nfolds = 5))

prop.func <- create.propensity.function(crossfit = TRUE,
                                        nfolds.crossfit = 10,
                                        cv.glmnet.args = list(type.measure = "auc",
                                                              nfolds = 5))

subgrp.model <- fit.subgroup(x = x, y = y,
                             trt = trt01,
                             propensity.func = prop.func,
                             augment.func = aug.func,
                             loss   = "sq_loss_lasso",
                             nfolds = 10)    # option for cv.glmnet (for ITR estimation)

summary(subgrp.model)




cleanEx()
nameEx("create.propensity.function")
### * create.propensity.function

flush(stderr()); flush(stdout())

### Name: create.propensity.function
### Title: Creation of propensity fitting function
### Aliases: create.propensity.function

### ** Examples

library(personalized)

set.seed(123)
n.obs  <- 500
n.vars <- 15
x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


# simulate non-randomized treatment
xbetat   <- 0.5 + 0.5 * x[,7] - 0.5 * x[,9]
trt.prob <- exp(xbetat) / (1 + exp(xbetat))
trt01    <- rbinom(n.obs, 1, prob = trt.prob)

trt      <- 2 * trt01 - 1

# simulate response
# delta below drives treatment effect heterogeneity
delta <- 2 * (0.5 + x[,2] - x[,3] - x[,11] + x[,1] * x[,12] )
xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13] + 0.5 * x[,15] ^ 2
xbeta <- xbeta + delta * trt

# continuous outcomes
y <- drop(xbeta) + rnorm(n.obs, sd = 2)

aug.func <- create.augmentation.function(family = "gaussian",
                                         crossfit = TRUE,
                                         nfolds.crossfit = 10,
                                         cv.glmnet.args = list(type.measure = "mae",
                                                               nfolds = 5))

prop.func <- create.propensity.function(crossfit = TRUE,
                                        nfolds.crossfit = 10,
                                        cv.glmnet.args = list(type.measure = "mae",
                                                              nfolds = 5))

subgrp.model <- fit.subgroup(x = x, y = y,
                             trt = trt01,
                             propensity.func = prop.func,
                             augment.func = aug.func,
                             loss   = "sq_loss_lasso",
                             nfolds = 10)    # option for cv.glmnet (for ITR estimation)

summary(subgrp.model)




cleanEx()
nameEx("fit.subgroup")
### * fit.subgroup

flush(stderr()); flush(stdout())

### Name: fit.subgroup
### Title: Fitting subgroup identification models
### Aliases: fit.subgroup

### ** Examples

library(personalized)

set.seed(123)
n.obs  <- 500
n.vars <- 15
x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


# simulate non-randomized treatment
xbetat   <- 0.5 + 0.5 * x[,7] - 0.5 * x[,9]
trt.prob <- exp(xbetat) / (1 + exp(xbetat))
trt01    <- rbinom(n.obs, 1, prob = trt.prob)

trt      <- 2 * trt01 - 1

# simulate response
# delta below drives treatment effect heterogeneity
delta <- 2 * (0.5 + x[,2] - x[,3] - x[,11] + x[,1] * x[,12] )
xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13] + 0.5 * x[,15] ^ 2
xbeta <- xbeta + delta * trt

# continuous outcomes
y <- drop(xbeta) + rnorm(n.obs, sd = 2)

# binary outcomes
y.binary <- 1 * (xbeta + rnorm(n.obs, sd = 2) > 0 )

# count outcomes
y.count <- round(abs(xbeta + rnorm(n.obs, sd = 2)))

# time-to-event outcomes
surv.time <- exp(-20 - xbeta + rnorm(n.obs, sd = 1))
cens.time <- exp(rnorm(n.obs, sd = 3))
y.time.to.event  <- pmin(surv.time, cens.time)
status           <- 1 * (surv.time <= cens.time)

# create function for fitting propensity score model
prop.func <- function(x, trt)
{
    # fit propensity score model
    propens.model <- cv.glmnet(y = trt,
                               x = x, family = "binomial")
    pi.x <- predict(propens.model, s = "lambda.min",
                    newx = x, type = "response")[,1]
    pi.x
}


####################  Continuous outcomes ################################


subgrp.model <- fit.subgroup(x = x, y = y,
                           trt = trt01,
                           propensity.func = prop.func,
                           loss   = "sq_loss_lasso",
                           nfolds = 10)              # option for cv.glmnet

summary(subgrp.model)

# estimates of the individual-specific
# treatment effect estimates:
subgrp.model$individual.trt.effects

# fit lasso + gam model with REML option for gam


####################  Using an augmentation function #####################
## augmentation funcions involve modeling the conditional mean E[Y|T, X]
## and returning predictions that are averaged over the treatment values
## return <- 1/2 * (hat{E}[Y|T=1, X] + hat{E}[Y|T=-1, X])
##########################################################################

augment.func <- function(x, y, trt) {
    data <- data.frame(x, y, trt)
    xm <- model.matrix(y~trt*x-1, data = data)

    lmod <- cv.glmnet(y = y, x = xm)
    ## get predictions when trt = 1
    data$trt <- 1
    xm <- model.matrix(y~trt*x-1, data = data)
    preds_1  <- predict(lmod, xm, s = "lambda.min")

    ## get predictions when trt = -1
    data$trt <- -1
    xm <- model.matrix(y~trt*x-1, data = data)
    preds_n1  <- predict(lmod, xm, s = "lambda.min")

    ## return predictions averaged over trt
    return(0.5 * (preds_1 + preds_n1))
}


####################  Binary outcomes ####################################

# use logistic loss for binary outcomes
subgrp.model.bin <- fit.subgroup(x = x, y = y.binary,
                           trt = trt01,
                           propensity.func = prop.func,
                           loss   = "logistic_loss_lasso",
                           type.measure = "auc",    # option for cv.glmnet
                           nfolds = 5)              # option for cv.glmnet

subgrp.model.bin


####################  Count outcomes #####################################

# use poisson loss for count/poisson outcomes
subgrp.model.poisson <- fit.subgroup(x = x, y = y.count,
                           trt = trt01,
                           propensity.func = prop.func,
                           loss   = "poisson_loss_lasso",
                           type.measure = "mse",    # option for cv.glmnet
                           nfolds = 5)              # option for cv.glmnet

subgrp.model.poisson


####################  Time-to-event outcomes #############################

library(survival)


####################  Using custom loss functions ########################

## Use custom loss function for binary outcomes

fit.custom.loss.bin <- function(x, y, weights, offset, ...) {
    df <- data.frame(y = y, x)

    # minimize logistic loss with NO lasso penalty
    # with allowance for efficiency augmentation
    glmf <- glm(y ~ x - 1, weights = weights,
                offset = offset, # offset term allows for efficiency augmentation
                family = binomial(), ...)

    # save coefficients
    cfs = coef(glmf)

    # create prediction function.
    prd = function(x, type = "response") {
         dfte <- cbind(1, x)
         colnames(dfte) <- names(cfs)
         ## predictions must be returned on the scale
         ## of the linear predictor
         predict(glmf, data.frame(dfte), type = "link")
    }
    # return lost of required components
    list(predict = prd, model = glmf, coefficients = cfs)
}



## try exponential loss for
## positive outcomes

fit.expo.loss <- function(x, y, weights, ...)
{
    expo.loss <- function(beta, x, y, weights) {
        sum(weights * y * exp(-drop(x %*% beta)))
    }

    # use optim() to minimize loss function
    opt <- optim(rep(0, NCOL(x)), fn = expo.loss, x = x, y = y, weights = weights)

    coefs <- opt$par

    pred <- function(x, type = "response") {
        tcrossprod(cbind(1, x), t(coefs))
    }

    # return list of required components
    list(predict = pred, model = opt, coefficients = coefs)
}






cleanEx()
nameEx("plot")
### * plot

flush(stderr()); flush(stdout())

### Name: plot.subgroup_fitted
### Title: Plotting results for fitted subgroup identification models
### Aliases: plot.subgroup_fitted plot.subgroup_validated

### ** Examples

library(personalized)

set.seed(123)
n.obs  <- 250
n.vars <- 15
x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


# simulate non-randomized treatment
xbetat   <- 0.5 + 0.5 * x[,11] - 0.5 * x[,13]
trt.prob <- exp(xbetat) / (1 + exp(xbetat))
trt01    <- rbinom(n.obs, 1, prob = trt.prob)

trt      <- 2 * trt01 - 1

# simulate response
delta <- 2 * (0.5 + x[,2] - x[,3] - x[,11] + x[,1] * x[,12])
xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13]
xbeta <- xbeta + delta * trt

# continuous outcomes
y <- drop(xbeta) + rnorm(n.obs, sd = 2)

# create function for fitting propensity score model
prop.func <- function(x, trt)
{
    # fit propensity score model
    propens.model <- cv.glmnet(y = trt,
                               x = x, family = "binomial")
    pi.x <- predict(propens.model, s = "lambda.min",
                    newx = x, type = "response")[,1]
    pi.x
}

subgrp.model <- fit.subgroup(x = x, y = y,
                           trt = trt01,
                           propensity.func = prop.func,
                           loss   = "sq_loss_lasso",
                           nfolds = 5)              # option for cv.glmnet

subgrp.model$subgroup.trt.effects

plot(subgrp.model)

plot(subgrp.model, type = "boxplot")

plot(subgrp.model, type = "interaction")

plot(subgrp.model, type = "conditional")

valmod <- validate.subgroup(subgrp.model, B = 3,
                          method = "training_test",
                          benefit.score.quantiles = c(0.25, 0.5, 0.75),
                          train.fraction = 0.75)

plot(valmod)


plot(valmod, type = "interaction")

# see how summary statistics of subgroups change
# when the subgroups are defined based on different cutoffs
# (25th quantile of bene score, 50th, and 75th)
plot(valmod, type = "conditional")

# visualize the frequency of particular variables
# of being selected across the resampling iterations with
# 'type = "stability"'
# not run:
# plot(valmod, type = "stability")




cleanEx()
nameEx("plotCompare")
### * plotCompare

flush(stderr()); flush(stdout())

### Name: plotCompare
### Title: Plot a comparison results for fitted or validated subgroup
###   identification models
### Aliases: plotCompare

### ** Examples

library(personalized)

set.seed(123)
n.obs  <- 100
n.vars <- 15
x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


# simulate non-randomized treatment
xbetat   <- 0.5 + 0.5 * x[,1] - 0.5 * x[,4]
trt.prob <- exp(xbetat) / (1 + exp(xbetat))
trt01    <- rbinom(n.obs, 1, prob = trt.prob)

trt      <- 2 * trt01 - 1

# simulate response
delta <- 2 * (0.5 + x[,2] - x[,3] - x[,11] + x[,1] * x[,12])
xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13]
xbeta <- xbeta + delta * trt

# continuous outcomes
y <- drop(xbeta) + rnorm(n.obs, sd = 2)

# create function for fitting propensity score model
prop.func <- function(x, trt)
{
    # fit propensity score model
    propens.model <- cv.glmnet(y = trt,
                               x = x, family = "binomial")
    pi.x <- predict(propens.model, s = "lambda.min",
                    newx = x, type = "response")[,1]
    pi.x
}

subgrp.model <- fit.subgroup(x = x, y = y,
                           trt = trt01,
                           propensity.func = prop.func,
                           loss   = "sq_loss_lasso",
                           nfolds = 5)              # option for cv.glmnet


subgrp.model.o <- fit.subgroup(x = x, y = y,
                           trt = trt01,
                           propensity.func = prop.func,
                           loss   = "owl_logistic_flip_loss_lasso",
                           nfolds = 5)

plotCompare(subgrp.model, subgrp.model.o)




cleanEx()
nameEx("predict")
### * predict

flush(stderr()); flush(stdout())

### Name: predict.subgroup_fitted
### Title: Function to predict either benefit scores or treatment
###   recommendations
### Aliases: predict.subgroup_fitted predict.wksvm

### ** Examples

library(personalized)

set.seed(123)
n.obs  <- 1000
n.vars <- 50
x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


# simulate non-randomized treatment
xbetat   <- 0.5 + 0.5 * x[,21] - 0.5 * x[,41]
trt.prob <- exp(xbetat) / (1 + exp(xbetat))
trt01    <- rbinom(n.obs, 1, prob = trt.prob)

trt      <- 2 * trt01 - 1

# simulate response
delta <- 2 * (0.5 + x[,2] - x[,3] - x[,11] + x[,1] * x[,12])
xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13]
xbeta <- xbeta + delta * trt

# continuous outcomes
y <- drop(xbeta) + rnorm(n.obs, sd = 2)

# create function for fitting propensity score model
prop.func <- function(x, trt)
{
    # fit propensity score model
    propens.model <- cv.glmnet(y = trt,
                               x = x, family = "binomial")
    pi.x <- predict(propens.model, s = "lambda.min",
                    newx = x, type = "response")[,1]
    pi.x
}

subgrp.model <- fit.subgroup(x = x, y = y,
                            trt = trt01,
                            propensity.func = prop.func,
                            loss   = "sq_loss_lasso",
                            nfolds = 5)              # option for cv.glmnet

subgrp.model$subgroup.trt.effects
benefit.scores <- predict(subgrp.model, newx = x, type = "benefit.score")

rec.trt.grp <- predict(subgrp.model, newx = x, type = "trt.group")



cleanEx()
nameEx("summarize.subgroups")
### * summarize.subgroups

flush(stderr()); flush(stdout())

### Name: summarize.subgroups
### Title: Summarizing covariates within estimated subgroups
### Aliases: summarize.subgroups summarize.subgroups.default
###   summarize.subgroups.subgroup_fitted

### ** Examples

library(personalized)

set.seed(123)
n.obs  <- 1000
n.vars <- 50
x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


# simulate non-randomized treatment
xbetat   <- 0.5 + 0.5 * x[,21] - 0.5 * x[,41]
trt.prob <- exp(xbetat) / (1 + exp(xbetat))
trt01    <- rbinom(n.obs, 1, prob = trt.prob)

trt      <- 2 * trt01 - 1

# simulate response
delta <- 2 * (0.5 + x[,2] - x[,3] - x[,11] + x[,1] * x[,12])
xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13]
xbeta <- xbeta + delta * trt

# continuous outcomes
y <- drop(xbeta) + rnorm(n.obs, sd = 2)

# create function for fitting propensity score model
prop.func <- function(x, trt)
{
    # fit propensity score model
    propens.model <- cv.glmnet(y = trt,
                               x = x, family = "binomial")
    pi.x <- predict(propens.model, s = "lambda.min",
                    newx = x, type = "response")[,1]
    pi.x
}

subgrp.model <- fit.subgroup(x = x, y = y,
                             trt = trt01,
                             propensity.func = prop.func,
                             loss   = "sq_loss_lasso",
                             nfolds = 5)    # option for cv.glmnet

comp <- summarize.subgroups(subgrp.model)
print(comp, p.value = 0.01)

# or we can simply supply the matrix x and the subgroups
comp2 <- summarize.subgroups(x, subgroup = 1 * (subgrp.model$benefit.scores > 0))

print(comp2, p.value = 0.01)




cleanEx()
nameEx("treatment.effects")
### * treatment.effects

flush(stderr()); flush(stdout())

### Name: treatment.effects
### Title: Calculation of covariate-conditional treatment effects
### Aliases: treatment.effects treatment.effects.default treat.effects
###   treatment.effects.subgroup_fitted

### ** Examples

library(personalized)

set.seed(123)
n.obs  <- 1000
n.vars <- 50
x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


# simulate non-randomized treatment
xbetat   <- 0.5 + 0.5 * x[,21] - 0.5 * x[,41]
trt.prob <- exp(xbetat) / (1 + exp(xbetat))
trt01    <- rbinom(n.obs, 1, prob = trt.prob)

trt      <- 2 * trt01 - 1

# simulate response
delta <- 2 * (0.5 + x[,2] - x[,3] - x[,11] + x[,1] * x[,12])
xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13]
xbeta <- xbeta + delta * trt

# continuous outcomes
y <- drop(xbeta) + rnorm(n.obs, sd = 2)

# time-to-event outcomes
surv.time <- exp(-20 - xbeta + rnorm(n.obs, sd = 1))
cens.time <- exp(rnorm(n.obs, sd = 3))
y.time.to.event  <- pmin(surv.time, cens.time)
status           <- 1 * (surv.time <= cens.time)

# create function for fitting propensity score model
prop.func <- function(x, trt)
{
    # fit propensity score model
    propens.model <- cv.glmnet(y = trt,
                               x = x, family = "binomial")
    pi.x <- predict(propens.model, s = "lambda.min",
                    newx = x, type = "response")[,1]
    pi.x
}

subgrp.model <- fit.subgroup(x = x, y = y,
                             trt = trt01,
                             propensity.func = prop.func,
                             loss   = "sq_loss_lasso",
                             nfolds = 5)    # option for cv.glmnet

trt_eff <- treatment.effects(subgrp.model)
str(trt_eff)

trt_eff


library(survival)
subgrp.model.cox <- fit.subgroup(x = x, y = Surv(y.time.to.event, status),
                           trt = trt01,
                           propensity.func = prop.func,
                           loss   = "cox_loss_lasso",
                           nfolds = 5)              # option for cv.glmnet

trt_eff_c <- treatment.effects(subgrp.model.cox)
str(trt_eff_c)

trt_eff_c




cleanEx()
nameEx("validate.subgroup")
### * validate.subgroup

flush(stderr()); flush(stdout())

### Name: validate.subgroup
### Title: Validating fitted subgroup identification models
### Aliases: validate.subgroup

### ** Examples

library(personalized)

set.seed(123)
n.obs  <- 500
n.vars <- 20
x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


# simulate non-randomized treatment
xbetat   <- 0.5 + 0.5 * x[,11] - 0.5 * x[,13]
trt.prob <- exp(xbetat) / (1 + exp(xbetat))
trt01    <- rbinom(n.obs, 1, prob = trt.prob)

trt      <- 2 * trt01 - 1

# simulate response
delta <- 2 * (0.5 + x[,2] - x[,3] - x[,11] + x[,1] * x[,12])
xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13]
xbeta <- xbeta + delta * trt

# continuous outcomes
y <- drop(xbeta) + rnorm(n.obs, sd = 2)

# create function for fitting propensity score model
prop.func <- function(x, trt)
{
    # fit propensity score model
    propens.model <- cv.glmnet(y = trt,
                               x = x, family = "binomial")
    pi.x <- predict(propens.model, s = "lambda.min",
                    newx = x, type = "response")[,1]
    pi.x
}

subgrp.model <- fit.subgroup(x = x, y = y,
                             trt = trt01,
                             propensity.func = prop.func,
                             loss   = "sq_loss_lasso",
                             nfolds = 5)    # option for cv.glmnet


x.test <- matrix(rnorm(10 * n.obs * n.vars, sd = 3), 10 * n.obs, n.vars)


# simulate non-randomized treatment
xbetat.test   <- 0.5 + 0.5 * x.test[,11] - 0.5 * x.test[,13]
trt.prob.test <- exp(xbetat.test) / (1 + exp(xbetat.test))
trt01.test    <- rbinom(10 * n.obs, 1, prob = trt.prob.test)

trt.test      <- 2 * trt01.test - 1

# simulate response
delta.test <- 2 * (0.5 + x.test[,2] - x.test[,3] - x.test[,11] + x.test[,1] * x.test[,12])
xbeta.test <- x.test[,1] + x.test[,11] - 2 * x.test[,12]^2 + x.test[,13]
xbeta.test <- xbeta.test + delta.test * trt.test

y.test <- drop(xbeta.test) + rnorm(10 * n.obs, sd = 2)

valmod <- validate.subgroup(subgrp.model, B = 3,
                            method = "training_test",
                            train.fraction = 0.75)
valmod

print(valmod, which.quant = c(4, 5))




cleanEx()
nameEx("weighted.ksvm")
### * weighted.ksvm

flush(stderr()); flush(stdout())

### Name: weighted.ksvm
### Title: Fit weighted kernel svm model.
### Aliases: weighted.ksvm

### ** Examples


library(kernlab)

x <- matrix(rnorm(200 * 2), ncol = 2)

y <- 2 * (sin(x[,2]) ^ 2 * exp(-x[,2]) - 0.2 > rnorm(200, sd = 0.1)) - 1

weights <- runif(100, max = 1.5, min = 0.5)

wk <- weighted.ksvm(x = x[1:100,], y = y[1:100], C = c(0.1, 0.5, 1, 2, 10),
                    weights = weights[1:100])

pr <- predict(wk, newx = x[101:200,])

mean(pr == y[101:200])




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
