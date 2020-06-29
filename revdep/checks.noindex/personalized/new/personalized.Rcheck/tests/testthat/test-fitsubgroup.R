

context("fit.subgroup")

test_that("test fit.subgroup for continuous outcomes and various losses", {

    library(kernlab)
    set.seed(123)
    n.obs  <- 100
    n.vars <- 5
    x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


    # simulate non-randomized treatment
    xbetat   <- 0.5 + 0.5 * x[,1] - 0.5 * x[,5]
    trt.prob <- exp(xbetat) / (1 + exp(xbetat))
    trt01    <- rbinom(n.obs, 1, prob = trt.prob)

    trt      <- 2 * trt01 - 1

    # simulate response
    delta <- 2 * (0.5 + x[,2] - x[,3]  )
    xbeta <- x[,1]
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

    prop.func2 <- function(x, trt)
    {
        # fit propensity score model
        propens.model <- cv.glmnet(y = trt,
                                   x = x, family = "binomial")
        pi.x <- predict(propens.model, s = "lambda.min",
                        newx = x, type = "response")
        pi.x
    }

    prop.func3 <- function(x, trt)
    {
        # fit propensity score model
        propens.model <- cv.glmnet(y = trt,
                                   x = x, family = "binomial")
        pi.x <- predict(propens.model, s = "lambda.min",
                        newx = x, type = "response")[,1]
        dim(pi.x) <- NROW(pi.x)
        pi.x
    }

    prop.func.bad <- function(x, trt, something)
    {
        return(numeric(NROW(trt)))
    }

    prop.func.bad2 <- function(x, trt, something, somethingelse)
    {
        return(numeric(NROW(trt)))
    }

    prop.func.bad3 <- function(x, something)
    {
        return(numeric(NROW(x)))
    }

    prop.func.bad4 <- function(x, trt)
    {
        ret <- numeric(NROW(x))
        ret[5] <- 100
        ret
    }

    if (Sys.info()[[1]] != "windows")
    {

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func.bad,
                                  loss   = "sq_loss_lasso",
                                  nfolds = 5))

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func.bad2,
                                  loss   = "sq_loss_lasso",
                                  nfolds = 5))
        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func.bad3,
                                  loss   = "sq_loss_lasso",
                                  nfolds = 5))
        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func.bad4,
                                  loss   = "sq_loss_lasso",
                                  nfolds = 5))

        expect_error(fit.subgroup(x = y, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func))

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))

        print(subgrp.model)
    }

    subgrp.model <- fit.subgroup(x = x, y = y,
                                 trt = trt01,
                                 propensity.func = prop.func,
                                 larger.outcome.better = FALSE,
                                 loss   = "sq_loss_lasso",
                                 nfolds = 5)              # option for cv.glmnet

    print(subgrp.model)


    #### TEST CUSTOM LOSS FUNCTIONS

    fit.custom.loss <- function(x, y, weights, ...) {
        df <- data.frame(y = y, x)

        # minimize squared error loss with NO lasso penalty
        lmf <- lm(y ~ x - 1, weights = weights, ...)

        # save coefficients
        cfs = unname(coef(lmf))

        # create prediction function. Notice
        # how a column of 1's is appended
        # to ensure treatment main effects are included
        # in predictions
        prd = function(x, type = "response")
        {
            # roxygen2 removes the below, so
            # using tcrossprod instead
            #cbind(1, x)
            tcrossprod(cbind(1, x), t(cfs))
        }
        # return lost of required components
        list(predict = prd, model = lmf, coefficients = cfs)
    }


    fit.custom.loss2 <- function(x, y, weights, family, match.id,
                                 n.trts, offset, trt, ...) {
        df <- data.frame(y = y, x)

        # minimize squared error loss with NO lasso penalty
        lmf <- lm(y ~ x - 1, weights = weights, ...)

        # save coefficients
        cfs = unname(coef(lmf))

        # create prediction function. Notice
        # how a column of 1's is appended
        # to ensure treatment main effects are included
        # in predictions
        prd = function(x, type = "response")
        {
            # roxygen2 removes the below, so
            # using tcrossprod instead
            #cbind(1, x)
            tcrossprod(cbind(1, x), t(cfs))
        }
        # return lost of required components
        list(predict = prd, model = lmf, coefficients = cfs)
    }


    fit.custom.loss.bad <- function(x, y, weights, wrongarg, ...) {
        df <- data.frame(y = y, x)

        # minimize squared error loss with NO lasso penalty
        lmf <- lm(y ~ x - 1, weights = weights, ...)

        # save coefficients
        cfs = unname(coef(lmf))

        # create prediction function. Notice
        # how a column of 1's is appended
        # to ensure treatment main effects are included
        # in predictions
        prd = function(x, type = "response")
        {
            # roxygen2 removes the below, so
            # using tcrossprod instead
            #cbind(1, x)
            tcrossprod(cbind(1, x), t(cfs))
        }
        # return lost of required components
        list(predict = prd, model = lmf, coefficients = cfs)
    }


    fit.custom.loss.bad2 <- function(x, y, weights, wrongarg, wrongarg2, wrongarg3, wrongarg4, wrongarg5,
                                     just, toooo, many, arguments, ...) {
        df <- data.frame(y = y, x)

        # minimize squared error loss with NO lasso penalty
        lmf <- lm(y ~ x - 1, weights = weights, ...)

        # save coefficients
        cfs = unname(coef(lmf))

        # create prediction function. Notice
        # how a column of 1's is appended
        # to ensure treatment main effects are included
        # in predictions
        prd = function(x, type = "response")
        {
            # roxygen2 removes the below, so
            # using tcrossprod instead
            #cbind(1, x)
            tcrossprod(cbind(1, x), t(cfs))
        }
        # return lost of required components
        list(predict = prd, model = lmf, coefficients = cfs)
    }

    fit.custom.loss.bad3 <- function(x, y, weights) {
        df <- data.frame(y = y, x)

        # minimize squared error loss with NO lasso penalty
        lmf <- lm(y ~ x - 1, weights = weights, ...)

        # save coefficients
        cfs = unname(coef(lmf))

        # create prediction function. Notice
        # how a column of 1's is appended
        # to ensure treatment main effects are included
        # in predictions
        prd = function(x, type = "response")
        {
            # roxygen2 removes the below, so
            # using tcrossprod instead
            #cbind(1, x)
            tcrossprod(cbind(1, x), t(cfs))
        }
        # return lost of required components
        list(predict = prd, model = lmf, coefficients = cfs)
    }


    fit.custom.loss.bad4 <- function(x, y, ...) {
        df <- data.frame(y = y, x)

        # minimize squared error loss with NO lasso penalty
        lmf <- lm(y ~ x - 1, weights = weights, ...)

        # save coefficients
        cfs = unname(coef(lmf))

        # create prediction function. Notice
        # how a column of 1's is appended
        # to ensure treatment main effects are included
        # in predictions
        prd = function(x, type = "response")
        {
            # roxygen2 removes the below, so
            # using tcrossprod instead
            #cbind(1, x)
            tcrossprod(cbind(1, x), t(cfs))
        }
        # return lost of required components
        list(predict = prd, model = lmf, coefficients = cfs)
    }


    fit.custom.loss.bad5 <- function(x, y, weights, ...) {
        df <- data.frame(y = y, x)

        # minimize squared error loss with NO lasso penalty
        lmf <- lm(y ~ x - 1, weights = weights, ...)

        # save coefficients
        cfs = unname(coef(lmf))

        # create prediction function. Notice
        # how a column of 1's is appended
        # to ensure treatment main effects are included
        # in predictions
        prd = function(x, type = "response")
        {
            # roxygen2 removes the below, so
            # using tcrossprod instead
            #cbind(1, x)
            tcrossprod(cbind(1, x), t(cfs))
        }
        # return lost of required components
        list(predict = prd, model = lmf, badelement = 1000)
    }

    fit.custom.loss.bad6 <- function(x, y, weights, ...) {
        df <- data.frame(y = y, x)

        # minimize squared error loss with NO lasso penalty
        lmf <- lm(y ~ x - 1, weights = weights, ...)

        # save coefficients
        cfs = unname(coef(lmf))

        # create prediction function. Notice
        # how a column of 1's is appended
        # to ensure treatment main effects are included
        # in predictions
        prd = function(x, type = "response")
        {
            # roxygen2 removes the below, so
            # using tcrossprod instead
            #cbind(1, x)
            tcrossprod(cbind(1, x), t(cfs))
        }
        # return lost of required components
        list(predict = prd, badelement = 1000)
    }


    ## predict element returned not a function
    fit.custom.loss.bad7 <- function(x, y, weights, ...) {
        df <- data.frame(y = y, x)

        # minimize squared error loss with NO lasso penalty
        lmf <- lm(y ~ x - 1, weights = weights, ...)

        # save coefficients
        cfs = unname(coef(lmf))

        # create prediction function. Notice
        # how a column of 1's is appended
        # to ensure treatment main effects are included
        # in predictions
        prd = function(x, type = "response")
        {
            # roxygen2 removes the below, so
            # using tcrossprod instead
            #cbind(1, x)
            tcrossprod(cbind(1, x), t(cfs))
        }
        # return lost of required components
        list(predict = cfs, model = lmf, coefficients = cfs)
    }

    ## predict function returns something that's not the right length
    fit.custom.loss.bad8 <- function(x, y, weights, ...) {
        df <- data.frame(y = y, x)

        # minimize squared error loss with NO lasso penalty
        lmf <- lm(y ~ x - 1, weights = weights, ...)

        # save coefficients
        cfs = unname(coef(lmf))

        # create prediction function. Notice
        # how a column of 1's is appended
        # to ensure treatment main effects are included
        # in predictions
        prd = function(x, type = "response")
        {
            # roxygen2 removes the below, so
            # using tcrossprod instead
            #cbind(1, x)
            tcrossprod(cbind(1, x)[1:10,], t(cfs))
        }
        # return lost of required components
        list(predict = prd, model = lmf, coefficients = cfs)
    }




    fit.custom.loss.bin <- function(x, y, weights, offset, ...) {
        df <- data.frame(y = y, x)

        # minimize squared error loss with NO lasso penalty
        glmf <- glm(y ~ x - 1, weights = weights,
                    offset = offset, # offset term allows for efficiency augmentation
                    family = binomial(), ...)

        # save coefficients
        cfs = unname(coef(glmf))

        # create prediction function.
        prd = function(x, type = "response")
        {
            #cbind(1, x)
            tcrossprod(cbind(1, x), t(cfs))
        }
        # return lost of required components
        list(predict = prd, model = glmf, coefficients = cfs)
    }

    if (Sys.info()[[1]] != "windows")
    {
        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     fit.custom.loss = fit.custom.loss,
                                     loss   = "custom")

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))

        trt01_bad <- trt01
        trt01_bad[4] <- "THE TRTMENT"
        expect_error(subgrp.model <- fit.subgroup(x = x, y = y,
                                                  trt = trt01_bad,
                                                  propensity.func = prop.func,
                                                  fit.custom.loss = fit.custom.loss,
                                                  loss   = "custom"))

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss = "custom",
                                     fit.custom.loss = fit.custom.loss)

        # use loss with all possible args provided
        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss = "custom",
                                     fit.custom.loss = fit.custom.loss2)


        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  loss = "custom"))

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  fit.custom.loss = fit.custom.loss.bad))

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  fit.custom.loss = fit.custom.loss.bad2))

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  fit.custom.loss = fit.custom.loss.bad3))

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  fit.custom.loss = fit.custom.loss.bad4))

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  fit.custom.loss = fit.custom.loss.bad5))

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  fit.custom.loss = fit.custom.loss.bad6))

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  fit.custom.loss = fit.custom.loss.bad7))

        expect_warning(expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  fit.custom.loss = fit.custom.loss.bad8)))
    }

    ###############################





    if (Sys.info()[[1]] != "windows")
    {
        # test for factor trt
        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = as.factor(trt01),
                                     propensity.func = prop.func,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)

        expect_is(subgrp.model, "subgroup_fitted")


        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = as.factor(trt01),
                                     propensity.func = prop.func,
                                     loss   = "owl_logistic_loss_lasso",
                                     nfolds = 5)

        expect_is(subgrp.model, "subgroup_fitted")


        expect_warning(fit.subgroup(x = x, y = y,
                                  trt = as.factor(trt01),
                                  propensity.func = prop.func,
                                  reference.trt = 2,
                                  loss   = "owl_logistic_loss_lasso",
                                  nfolds = 5))

        ## only 1 trt level
        expect_error(fit.subgroup(x = x, y = y,
                                    trt = as.factor(rep(1, NROW(y))),
                                    propensity.func = prop.func,
                                    loss   = "owl_logistic_loss_lasso",
                                    nfolds = 5))

        ## too many trt levels
        expect_error(fit.subgroup(x = x, y = y,
                                  trt = as.factor(1:NROW(y)),
                                  propensity.func = prop.func,
                                  loss   = "owl_logistic_loss_lasso",
                                  nfolds = 5))

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = as.factor(trt01),
                                     propensity.func = prop.func,
                                     loss   = "owl_logistic_flip_loss_lasso",
                                     nfolds = 5)

        expect_is(subgrp.model, "subgroup_fitted")
    }


    if (Sys.info()[[1]] != "windows")
    {
        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "owl_hinge_flip_loss",
                                     nfolds = 5)

        expect_is(subgrp.model, "subgroup_fitted")


        trt01_bad <- trt01
        trt01_bad[4] <- "third_trt_level"
        expect_error(subgrp.model <- fit.subgroup(x = x, y = y,
                                                  trt = trt01_bad,
                                                  propensity.func = prop.func,
                                                  loss   = "owl_hinge_flip_loss",
                                                  nfolds = 5))


        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "owl_hinge_loss",
                                     nfolds = 5)

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  loss   = "owl_hinge_loss",
                                  nfolds = -5) )

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "owl_hinge_flip_loss",
                                     kernel = "vanilladot")

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))


        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "owl_hinge_flip_loss",
                                     kernel = "polydot",
                                     kpar = list(degree = 3),
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "owl_hinge_flip_loss",
                                     kernel = "laplacedot",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "owl_hinge_flip_loss",
                                     kernel = "besseldot",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "owl_hinge_flip_loss",
                                     kernel = "anovadot",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = as.factor(trt01),
                                     propensity.func = prop.func,
                                     loss   = "owl_hinge_flip_loss",
                                     kernel = "splinedot",
                                     margin = 0.2,
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))
    }


    subgrp.model <- fit.subgroup(x = x, y = y,
                                 trt = trt01,
                                 larger.outcome.better = FALSE,
                                 propensity.func = prop.func,
                                 loss   = "sq_loss_lasso",
                                 nfolds = 5)              # option for cv.glmnet

    expect_is(subgrp.model, "subgroup_fitted")

    if (Sys.info()[[1]] != "windows")
    {
        # test if pi.x is a matrix with 1 column
        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func2,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        # test if pi.x is a matrix with 1 column
        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func3,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        # no prop func
        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "sq_loss_gam")

        expect_is(subgrp.model, "subgroup_fitted")
        invisible(capture.output(print(subgrp.model)))
        invisible(capture.output(summary(subgrp.model)))

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "sq_loss_lasso_gam",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model)))
        invisible(capture.output(summary(subgrp.model)))
    }


    if (Sys.info()[[1]] != "windows")
    {

        ## tests for gam argument specification
        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     method.gam = "REML",
                                     loss   = "sq_loss_gam")

        expect_is(subgrp.model, "subgroup_fitted")

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     method.gam = "REML",
                                     loss   = "sq_loss_lasso_gam")

        expect_is(subgrp.model, "subgroup_fitted")

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     control = gam.control(epsilon = 1e-10),
                                     loss   = "sq_loss_gam")

        expect_is(subgrp.model, "subgroup_fitted")

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     optimizer = "bfgs",
                                     loss   = "sq_loss_lasso_gam")

        expect_is(subgrp.model, "subgroup_fitted")

        # subgrp.model <- fit.subgroup(x = x, y = y,
        #                              trt = trt01,
        #                              propensity.func = prop.func,
        #                              loss   = "sq_loss_gbm",
        #                              n.trees = 5,
        #                              n.cores = 1)
        #
        # invisible(capture.output(print(subgrp.model)))
        # invisible(capture.output(summary(subgrp.model)))
        # expect_is(subgrp.model, "subgroup_fitted")
    }

    # subgrp.model <- fit.subgroup(x = x, y = y,
    #                              trt = trt01,
    #                              propensity.func = prop.func,
    #                              loss   = "abs_loss_gbm",
    #                              n.trees = 5,
    #                              n.cores = 1)
    #
    # invisible(capture.output(print(subgrp.model)))
    # invisible(capture.output(summary(subgrp.model)))
    # expect_is(subgrp.model, "subgroup_fitted")
})



test_that("test fit.subgroup for time-to-event outcomes and various losses", {
    set.seed(123)
    n.obs  <- 100
    n.vars <- 5
    x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


    # simulate non-randomized treatment
    xbetat   <- 0.5 + 0.5 * x[,1] - 0.5 * x[,5]
    trt.prob <- exp(xbetat) / (1 + exp(xbetat))
    trt01    <- rbinom(n.obs, 1, prob = trt.prob)

    trt      <- 2 * trt01 - 1

    # simulate response
    delta <- 2 * (0.5 + x[,2] - x[,3]  )
    xbeta <- x[,1]
    xbeta <- xbeta + delta * trt

    # continuous outcomes
    y <- drop(xbeta) + rnorm(n.obs, sd = 2)

    # binary outcomes
    y.binary <- 1 * (xbeta + rnorm(n.obs, sd = 2) > 0 )

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

    prop.func2 <- function(x, trt)
    {
        # fit propensity score model
        propens.model <- cv.glmnet(y = trt,
                                   x = x, family = "binomial")
        pi.x <- predict(propens.model, s = "lambda.min",
                        newx = x, type = "response")
        pi.x
    }

    subgrp.model <- fit.subgroup(x = x, y = Surv(y.time.to.event, status),
                                 trt = trt01,
                                 propensity.func = prop.func,
                                 loss   = "cox_loss_lasso",
                                 nfolds = 5)              # option for cv.glmnet

    expect_is(subgrp.model, "subgroup_fitted")

    invisible(capture.output(print(subgrp.model, digits = 2)))

    invisible(capture.output(summary(subgrp.model)))


    if (Sys.info()[[1]] != "windows")
    {

        # test for factor trt
        subgrp.model <- fit.subgroup(x = x, y = Surv(y.time.to.event, status),
                                     trt = as.factor(trt01),
                                     propensity.func = prop.func,
                                     loss   = "cox_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")


        subgrp.model <- fit.subgroup(x = x, y = Surv(y.time.to.event, status),
                                     trt = trt01,
                                     larger.outcome.better = FALSE,
                                     propensity.func = prop.func,
                                     loss   = "cox_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        # test if pi.x is a matrix with 1 column
        subgrp.model <- fit.subgroup(x = x, y = Surv(y.time.to.event, status),
                                     trt = trt01,
                                     propensity.func = prop.func2,
                                     loss   = "cox_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")
    }

})







test_that("test fit.subgroup with augment.func for continuous outcomes and various losses", {
    set.seed(123)
    n.obs  <- 100
    n.vars <- 5
    x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


    # simulate non-randomized treatment
    xbetat   <- 0.5 + 0.5 * x[,1] - 0.5 * x[,5]
    trt.prob <- exp(xbetat) / (1 + exp(xbetat))
    trt01    <- rbinom(n.obs, 1, prob = trt.prob)

    trt      <- 2 * trt01 - 1

    # simulate response
    delta <- 2 * (0.5 + x[,2] - x[,3]  )
    xbeta <- x[,1]
    xbeta <- xbeta + delta * trt

    # continuous outcomes
    y <- drop(xbeta) + rnorm(n.obs, sd = 2)

    # count outcomes
    y.count <- round(abs(drop(xbeta) + rnorm(n.obs, sd = 2)))

    # binary outcomes
    y.binary <- 1 * (xbeta + rnorm(n.obs, sd = 2) > 0 )

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

    augment.func <- function(x, y) {lmod <- lm(y ~ x); return(fitted(lmod))}
    augment.func2 <- function(x, y, trt) {lmod <- lm(y ~ x + trt); return(fitted(lmod))}
    augment.func.bad <- function(x, y, something) {lmod <- lm(y ~ x); return(fitted(lmod))}
    augment.func.bad2 <- function(x, y, something, somethingelse) {lmod <- lm(y ~ x); return(fitted(lmod))}
    augment.func.bad3 <- function(x, something) {lmod <- lm(y ~ x); return(fitted(lmod))}
    augment.func.bad4 <- function(x, y) {lmod <- lm(y ~ x); return(fitted(lmod)[1:10])}

    subgrp.model <- fit.subgroup(x = x, y = y,
                                 trt = trt01,
                                 augment.func = augment.func,
                                 propensity.func = prop.func,
                                 loss   = "sq_loss_lasso",
                                 nfolds = 5)              # option for cv.glmnet

    expect_is(subgrp.model, "subgroup_fitted")




    if (Sys.info()[[1]] != "windows")
    {

        ## test for method which uses adjustment directly
        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     augment.func = augment.func,
                                     propensity.func = prop.func,
                                     loss   = "owl_logistic_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     augment.func = augment.func2,
                                     propensity.func = prop.func,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  augment.func = augment.func.bad,
                                  propensity.func = prop.func,
                                  loss   = "sq_loss_lasso",
                                  nfolds = 5))

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  augment.func = augment.func.bad2,
                                  propensity.func = prop.func,
                                  loss   = "sq_loss_lasso",
                                  nfolds = 5))

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  augment.func = augment.func.bad3,
                                  propensity.func = prop.func,
                                  loss   = "sq_loss_lasso",
                                  nfolds = 5))

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  augment.func = augment.func.bad4,
                                  propensity.func = prop.func,
                                  loss   = "sq_loss_lasso",
                                  nfolds = 5))


        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     augment.func = augment.func,
                                     loss   = "sq_loss_gam")

        expect_is(subgrp.model, "subgroup_fitted")
        invisible(capture.output(print(subgrp.model)))
        invisible(capture.output(summary(subgrp.model)))

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     augment.func = augment.func,
                                     loss   = "sq_loss_lasso_gam",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model)))
        invisible(capture.output(summary(subgrp.model)))

        # subgrp.model <- fit.subgroup(x = x, y = y,
        #                              trt = trt01,
        #                              propensity.func = prop.func,
        #                              loss   = "sq_loss_gbm",
        #                              n.trees = 5,
        #                              n.cores = 1)
        #
        # invisible(capture.output(print(subgrp.model)))
        # invisible(capture.output(summary(subgrp.model)))
        # expect_is(subgrp.model, "subgroup_fitted")

        # subgrp.model <- fit.subgroup(x = x, y = y,
        #                              trt = trt01,
        #                              propensity.func = prop.func,
        #                              loss   = "abs_loss_gbm",
        #                              n.trees = 5,
        #                              n.cores = 1)
        #
        # invisible(capture.output(print(subgrp.model)))
        # invisible(capture.output(summary(subgrp.model)))
        # expect_is(subgrp.model, "subgroup_fitted")

        expect_warning(fit.subgroup(x = x, y = y,
                                    trt = trt01,
                                    propensity.func = prop.func,
                                    loss   = "sq_loss_gbm",
                                    n.trees = 5,
                                    cv.folds = 1,
                                    n.cores = 1))

        # expect_warning(fit.subgroup(x = x, y = y,
        #                           trt = trt01,
        #                           propensity.func = prop.func,
        #                           loss   = "abs_loss_gbm",
        #                           n.trees = 5,
        #                           cv.folds = 1,
        #                           n.cores = 1))

        subgrp.model <- fit.subgroup(x = x, y = y.binary,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "logistic_loss_gbm",
                                     n.trees = 5,
                                     n.cores = 1)

        invisible(capture.output(print(subgrp.model)))
        invisible(capture.output(summary(subgrp.model)))
        expect_is(subgrp.model, "subgroup_fitted")

        expect_warning(fit.subgroup(x = x, y = y.binary,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  loss   = "logistic_loss_gbm",
                                  n.trees = 5,
                                  cv.folds = 1,
                                  n.cores = 1))
    }

    if (Sys.info()[[1]] != "windows")
    {
        subgrp.model <- fit.subgroup(x = x, y = Surv(y.time.to.event, status),
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "cox_loss_gbm",
                                     n.trees = 5,
                                     n.cores = 1)

        invisible(capture.output(print(subgrp.model)))
        invisible(capture.output(summary(subgrp.model)))
        expect_is(subgrp.model, "subgroup_fitted")

        expect_warning(fit.subgroup(x = x, y = Surv(y.time.to.event, status),
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  loss   = "cox_loss_gbm",
                                  n.trees = 5,
                                  cv.folds = 1,
                                  n.cores = 1))
    }
})


test_that("test fit.subgroup for continuous outcomes with no propensity function specified", {
    set.seed(123)
    n.obs  <- 100
    n.vars <- 5
    x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


    # simulate non-randomized treatment
    xbetat   <- 0.5 + 0.5 * x[,1] - 0.5 * x[,5]
    trt.prob <- exp(xbetat) / (1 + exp(xbetat))
    trt01    <- rbinom(n.obs, 1, prob = trt.prob)

    trt      <- 2 * trt01 - 1

    # simulate response
    delta <- 2 * (0.5 + x[,2] - x[,3]  )
    xbeta <- x[,1]
    xbeta <- xbeta + delta * trt

    # continuous outcomes
    y <- drop(xbeta) + rnorm(n.obs, sd = 2)

    # binary outcomes
    y.binary <- 1 * (xbeta + rnorm(n.obs, sd = 2) > 0 )

    # time-to-event outcomes
    surv.time <- exp(-20 - xbeta + rnorm(n.obs, sd = 1))
    cens.time <- exp(rnorm(n.obs, sd = 3))
    y.time.to.event  <- pmin(surv.time, cens.time)
    status           <- 1 * (surv.time <= cens.time)

    subgrp.model <- fit.subgroup(x = x, y = y,
                                 trt = trt01,
                                 loss   = "sq_loss_lasso",
                                 nfolds = 5)              # option for cv.glmnet

    expect_is(subgrp.model, "subgroup_fitted")

    invisible(capture.output(print(subgrp.model, digits = 2)))

    invisible(capture.output(summary(subgrp.model)))

})




test_that("test fit.subgroup for binary outcomes and various losses", {
    set.seed(123)
    n.obs  <- 100
    n.vars <- 5
    x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


    # simulate non-randomized treatment
    xbetat   <- 0.5 + 0.5 * x[,1] - 0.5 * x[,5]
    trt.prob <- exp(xbetat) / (1 + exp(xbetat))
    trt01    <- rbinom(n.obs, 1, prob = trt.prob)

    trt      <- 2 * trt01 - 1

    # simulate response
    delta <- 2 * (0.5 + x[,2] - x[,3]  )
    xbeta <- x[,1] + x[,5]
    xbeta <- xbeta + delta * trt

    # continuous outcomes
    y <- drop(xbeta) + rnorm(n.obs, sd = 2)

    # binary outcomes
    y.binary <- 1 * (xbeta + rnorm(n.obs, sd = 2) > 0 )

    # count outcomes
    y.count <- round(abs(drop(xbeta) + rnorm(n.obs, sd = 2)))

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

    subgrp.model <- fit.subgroup(x = x, y = y.binary,
                                 trt = trt01,
                                 propensity.func = prop.func,
                                 loss   = "logistic_loss_lasso",
                                 nfolds = 5)              # option for cv.glmnet

    expect_is(subgrp.model, "subgroup_fitted")

    invisible(capture.output(print(subgrp.model, digits = 2)))

    invisible(capture.output(summary(subgrp.model)))


    if (Sys.info()[[1]] != "windows")
    {
        subgrp.model <- fit.subgroup(x = x, y = y.count,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "poisson_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))


        subgrp.model <- fit.subgroup(x = x, y = y.count,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "poisson_loss_gam",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))


        subgrp.model <- fit.subgroup(x = x, y = y.count,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "poisson_loss_lasso_gam")

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))

        # subgrp.model <- fit.subgroup(x = x, y = y.count,
        #                              trt = trt01,
        #                              propensity.func = prop.func,
        #                              loss   = "poisson_loss_gbm")
        #
        # expect_is(subgrp.model, "subgroup_fitted")
        #
        # invisible(capture.output(print(subgrp.model, digits = 2)))
        #
        # invisible(capture.output(summary(subgrp.model)))


        augment.func <- function(x, y) {
            lmod <- glm(y ~ x, family = binomial());
            return(predict(lmod, type = "link"))
        }

        subgrp.modela <- fit.subgroup(x = x, y = y.binary,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     augment.func = augment.func,
                                     loss   = "logistic_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.modela, "subgroup_fitted")

        invisible(capture.output(print(subgrp.modela, digits = 2)))

        invisible(capture.output(summary(subgrp.modela)))



        subgrp.modela <- fit.subgroup(x = x, y = y.binary,
                                      trt = trt01,
                                      propensity.func = prop.func,
                                      augment.func = augment.func,
                                      loss   = "logistic_loss_gam",
                                      nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.modela, "subgroup_fitted")

        invisible(capture.output(print(subgrp.modela, digits = 2)))

        invisible(capture.output(summary(subgrp.modela)))


        subgrp.modelg <- fit.subgroup(x = x, y = y.binary,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "logistic_loss_lasso_gam",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.modelg, "subgroup_fitted")

        invisible(capture.output(print(subgrp.modelg, digits = 2)))

        invisible(capture.output(summary(subgrp.modelg)))


        subgrp.modelga <- fit.subgroup(x = x, y = y.binary,
                                      trt = trt01,
                                      propensity.func = prop.func,
                                      augment.func = augment.func,
                                      loss   = "logistic_loss_lasso_gam",
                                      nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.modelga, "subgroup_fitted")

        invisible(capture.output(print(subgrp.modelga, digits = 2)))

        invisible(capture.output(summary(subgrp.modelga)))

        # subgrp.model <- fit.subgroup(x = x, y = y.binary,
        #                              trt = trt01,
        #                              propensity.func = prop.func,
        #                              loss   = "logistic_loss_gam")
        #
        # expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))
    }

    # subgrp.model <- fit.subgroup(x = x, y = y.binary,
    #                              trt = trt01,
    #                              propensity.func = prop.func,
    #                              loss   = "logistic_loss_gbm", n.trees = 5,
    #                              n.cores = 1)
    #
    # expect_is(subgrp.model, "subgroup_fitted")
    #
    # invisible(capture.output(print(subgrp.model, digits = 2)))
    #
    # invisible(capture.output(summary(subgrp.model)))
})





test_that("test fit.subgroup for continuous outcomes and match.id provided", {
    set.seed(123)
    n.obs  <- 1000
    n.vars <- 5
    x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


    # simulate non-randomized treatment
    xbetat   <- -3 + 0.5 * x[,1] - 0.5 * x[,5]
    trt.prob <- exp(xbetat) / (1 + exp(xbetat))
    trt01    <- rbinom(n.obs, 1, prob = trt.prob)

    trt      <- 2 * trt01 - 1

    # simulate response
    delta <- 2 * (0.5 + 5 * sin(x[,2] ^ 2) - x[,3]  )
    xbeta <- x[,1]
    xbeta <- xbeta + delta * trt

    # continuous outcomes
    y <- drop(xbeta) + rnorm(n.obs, sd = 2)

    # binary outcomes
    y.binary <- 1 * (xbeta + rnorm(n.obs, sd = 2) > 0 )

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
                                   x = x, family = "binomial",
                                   type.measure = "auc")
        pi.x <- predict(propens.model, s = "lambda.min",
                        newx = x, type = "response")[,1]
        pi.x
    }

    n.matches <- 2
    pscore <- prop.func(x, trt)

    match.mat <- array(NA, dim = c(sum(trt01), n.matches) )

    trt.idx <- which(trt01 == 1)
    ctrl.idx <- which(trt01 == 0)


    for (m in 1:n.matches)
    {
        for (i in 1:length(trt.idx))
        {
            best.match <- which.min(abs(pscore[trt.idx[i]] - pscore[ctrl.idx]))

            # remove match from pool
            ctrl.idx <- ctrl.idx[ctrl.idx != ctrl.idx[best.match]]

            match.mat[i,m] <- ctrl.idx[best.match]
        }
    }




    n.total <- nrow(match.mat) * (ncol(match.mat) + 1)



    match.id.mat <- matrix(nrow = nrow(match.mat), ncol = ncol(match.mat))


    for (cc in 1:ncol(match.mat))
    {
        match.id.mat[,cc] <- 1:nrow(match.mat)
    }


    match.id  <- c(1:nrow(match.mat), as.vector(match.id.mat))
    match.idx <- c(trt.idx, as.vector(match.mat))

    x.m <- x[match.idx,]
    y.m <- y[match.idx]
    trt.m <- trt01[match.idx]
    y.time.to.event.m <- y.time.to.event[match.idx]
    status.m <- status[match.idx]
    y.binary.m <- y.binary[match.idx]


    subgrp.model.m <- fit.subgroup(x = x.m, y = y.m,
                                   trt = trt.m,
                                   match.id = as.factor(match.id),
                                   loss   = "sq_loss_lasso",
                                   nfolds = 5)              # option for cv.glmnet

    expect_is(subgrp.model.m, "subgroup_fitted")

    invisible(capture.output(print(subgrp.model.m, digits = 2)))

    invisible(capture.output(summary(subgrp.model.m)))

    if (Sys.info()[[1]] != "windows")
    {
        subgrp.model.m <- fit.subgroup(x = x.m, y = y.m,
                                       trt = trt.m,
                                       match.id = match.id,
                                       loss   = "sq_loss_lasso",
                                       nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model.m, "subgroup_fitted")



        subgrp.model.m <- fit.subgroup(x = x.m, y = y.m,
                                       trt = trt.m,
                                       match.id = match.id,
                                       loss   = "owl_hinge_loss",
                                       margin = 5,
                                       nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model.m, "subgroup_fitted")

        subgrp.model.m <- fit.subgroup(x = x.m, y = y.m,
                                       trt = trt.m,
                                       match.id = match.id,
                                       loss   = "owl_hinge_flip_loss",
                                       margin = 5,
                                       nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model.m, "subgroup_fitted")

        subgrp.model.m <- fit.subgroup(x = x.m, y = y.m,
                                       trt = trt.m,
                                       match.id = match.id,
                                       penlty.factor = rep(1, ncol(x.m) + 1),
                                       loss   = "sq_loss_lasso_gam",
                                       nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model.m, "subgroup_fitted")

        expect_error(fit.subgroup(x = x.m, y = y.m,
                                  trt = trt.m,
                                  match.id = match.id,
                                  loss   = "sq_loss_lasso_gam",
                                  nfolds = 2))

        expect_warning(fit.subgroup(x = x.m, y = y.m,
                                    trt = trt.m,
                                    match.id = match.id,
                                    loss   = "sq_loss_lasso_gam",
                                    foldid = sample(1:5, nrow(x.m), replace = TRUE)))

        # subgrp.model.m <- fit.subgroup(x = x.m, y = y.m,
        #                                trt = trt.m,
        #                                match.id = match.id,
        #                                loss   = "sq_loss_gbm", n.trees = 5, n.cores = 1)
        #
        # expect_is(subgrp.model.m, "subgroup_fitted")

        # subgrp.model.m <- fit.subgroup(x = x.m, y = y.m,
        #                                trt = trt.m,
        #                                match.id = match.id,
        #                                loss   = "abs_loss_gbm", n.trees = 5, n.cores = 1)
        #
        # expect_is(subgrp.model.m, "subgroup_fitted")

        subgrp.model.m <- fit.subgroup(x = x.m, y = y.binary.m,
                                       trt = trt.m,
                                       match.id = match.id,
                                       loss   = "logistic_loss_gbm", n.trees = 5, n.cores = 1)

        expect_is(subgrp.model.m, "subgroup_fitted")

        subgrp.model.m <- fit.subgroup(x = x.m, y = y.m,
                                       trt = trt.m,
                                       match.id = as.factor(match.id),
                                       loss   = "sq_loss_lasso_gam")

        expect_is(subgrp.model.m, "subgroup_fitted")

        subgrp.model.m <- fit.subgroup(x = x.m, y = y.m,
                                       trt = trt.m,
                                       propensity.func = prop.func,
                                       match.id = as.factor(match.id),
                                       loss   = "sq_loss_lasso",
                                       nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model.m, "subgroup_fitted")

        subgrp.val <- validate.subgroup(subgrp.model.m, B = 10,
                                        method = "train")

        expect_is(subgrp.val, "subgroup_validated")

        invisible(capture.output(print(subgrp.val, digits = 2)))

        subgrp.val <- validate.subgroup(subgrp.model.m, B = 10,
                                        method = "boot")

        expect_is(subgrp.val, "subgroup_validated")

        invisible(capture.output(print(subgrp.val, digits = 2)))


        expect_warning(subgrp.model.m <- fit.subgroup(x = x.m, y = y.m,
                                       trt = trt.m,
                                       match.id = as.factor(match.id),
                                       loss   = "sq_loss_lasso",
                                       foldid = sample(1:5, nrow(x.m), replace = TRUE)))

        expect_is(subgrp.model.m, "subgroup_fitted")


        expect_warning(subgrp.model.m <- fit.subgroup(x = x.m, y = Surv(y.time.to.event.m, status.m),
                                                      trt = trt.m,
                                                      match.id = as.factor(match.id),
                                                      loss   = "cox_loss_lasso",
                                                      foldid = sample(1:5, nrow(x.m), replace = TRUE)))

        expect_is(subgrp.model.m, "subgroup_fitted")




        expect_error(fit.subgroup(x = x.m, y = y.m,
                                       trt = trt.m,
                                       match.id = as.factor(match.id),
                                       loss   = "sq_loss_lasso",
                                  nfolds = 2) )

        expect_error(fit.subgroup(x = x.m, y = y.m,
                                  trt = trt.m,
                                  propensity.func = prop.func,
                                  match.id = rep(1, length(y.m)), # provide bad match.id (only 1 level)
                                  loss   = "sq_loss_lasso",
                                  nfolds = 5) )

        expect_error(fit.subgroup(x = x.m, y = Surv(y.time.to.event.m, status.m),
                                  trt = trt.m,
                                  match.id = as.factor(match.id),
                                  loss   = "cox_loss_lasso",
                                  nfolds = 2) )

        subgrp.model.m <- fit.subgroup(x = x.m, y = y.m,
                                                      trt = trt.m,
                                                      match.id = as.factor(match.id),
                                                      loss   = "sq_loss_lasso")

        expect_is(subgrp.model.m, "subgroup_fitted")

        subgrp.model.m <- fit.subgroup(x = x.m, y = y.m,
                                       trt = trt.m,
                                       match.id = as.factor(match.id),
                                       penalty.factor = rep(1, ncol(x.m) + 1),
                                       loss   = "sq_loss_lasso")

        expect_is(subgrp.model.m, "subgroup_fitted")

        subgrp.model.m <- fit.subgroup(x = x.m, y = Surv(y.time.to.event.m, status.m),
                                       trt = trt.m,
                                       match.id = as.factor(match.id),
                                       loss   = "cox_loss_lasso")

        expect_is(subgrp.model.m, "subgroup_fitted")

        subgrp.model.m <- fit.subgroup(x = x.m, y = Surv(y.time.to.event.m, status.m),
                                       trt = trt.m,
                                       match.id = as.factor(match.id),
                                       loss   = "cox_loss_gbm", n.trees = 5, n.cores = 1)

        expect_is(subgrp.model.m, "subgroup_fitted")



        subgrp.model.m <- fit.subgroup(x = x.m, y = Surv(y.time.to.event.m, status.m),
                                       trt = trt.m,
                                       augment.func = function(x,y,trt) {return(rep(1, NROW(x)))},
                                       match.id = as.factor(match.id),
                                       penalty.factor = rep(1, ncol(x.m) + 1),
                                       loss   = "cox_loss_lasso")

        expect_is(subgrp.model.m, "subgroup_fitted")



        ## parallel

        subgrp.val <- validate.subgroup(subgrp.model.m, B = 10,
                                        parallel = TRUE,
                                        method = "training")


        subgrp.val <- validate.subgroup(subgrp.model.m, B = 10,
                                        parallel = TRUE,
                                        method = "boot")


    }

})












test_that("test fit.subgroup for continuous outcomes and multiple trts and various losses", {
    set.seed(123)
    n.obs  <- 100
    n.vars <- 5
    x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)

    # simulated non-randomized treatment with multiple levels
    xbetat_1   <- 0.15 + 0.5 * x[,1] - 0.25 * x[,5]
    xbetat_2   <- 0.15 - 0.5 * x[,2] + 0.25 * x[,3]
    trt.1.prob <- exp(xbetat_1) / (1 + exp(xbetat_1) + exp(xbetat_2))
    trt.2.prob <- exp(xbetat_2) / (1 + exp(xbetat_1) + exp(xbetat_2))
    trt.3.prob <- 1 - (trt.1.prob + trt.2.prob)
    prob.mat <- cbind(trt.1.prob, trt.2.prob, trt.3.prob)
    trt    <- apply(prob.mat, 1, function(rr) rmultinom(1, 1, prob = rr))
    trt    <- apply(trt, 2, function(rr) which(rr == 1))


    # simulate response
    delta1 <- 8 * (0.5 + x[,2] - x[,3]  )
    delta2 <- 4 * (-0.5 + x[,1] - x[,5] + x[,4]  )
    xbeta <- x[,1] - 2 * x[,2] - 3 * x[,5]
    xbeta <- xbeta + delta1 * ((trt == 1) - (trt == 3) ) + delta2 * ((trt == 2) - (trt == 3) )


    # continuous outcomes
    y <- drop(xbeta) + rnorm(n.obs, sd = 2)

    # binary outcomes
    y.binary <- 1 * (xbeta + rnorm(n.obs, sd = 2) > 0 )

    # time-to-event outcomes
    surv.time <- exp(-20 - xbeta + rnorm(n.obs, sd = 1))
    cens.time <- exp(rnorm(n.obs, sd = 3))
    y.time.to.event  <- pmin(surv.time, cens.time)
    status           <- 1 * (surv.time <= cens.time)

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

    # use multinomial logistic regression model with lasso penalty for propensity
    propensity.multinom.lasso.bad <- function(x, trt)
    {
        if (!is.factor(trt)) trt <- as.factor(trt)
        gfit <- cv.glmnet(y = trt, x = x, family = "multinomial")

        # predict returns a matrix of probabilities:
        # one column for each treatment level
        propens <- drop(predict(gfit, newx = x, type = "response", s = "lambda.min",
                                nfolds = 5, alpha = 0))

        # return the probability corresponding to the
        # treatment that was observed
        probs <- propens[,match(levels(trt), colnames(propens))][,1:2]

        probs
    }

    subgrp.model <- fit.subgroup(x = x, y = y,
                                 trt = trt,
                                 propensity.func = propensity.multinom.lasso,
                                 loss   = "sq_loss_lasso",
                                 nfolds = 5)              # option for cv.glmnet

    expect_is(subgrp.model, "subgroup_fitted")

    if (Sys.info()[[1]] != "windows")
    {
        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt,
                                     propensity.func = propensity.multinom.lasso,
                                     loss   = "owl_logistic_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))

        expect_warning(fit.subgroup(x = x, y = y,
                                    trt = trt,
                                    propensity.func = propensity.multinom.lasso,
                                    loss   = "owl_logistic_loss_lasso",
                                    method = "a_learning",
                                    nfolds = 5))

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt,
                                  propensity.func = propensity.multinom.lasso,
                                  loss   = "owl_logistic_flip_loss_lasso",
                                  nfolds = 5))


        ## test matching without prop score

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt,
                                     match.id = rep(1:10, each = 10),
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)
        expect_is(subgrp.model, "subgroup_fitted")

        invisible(capture.output(print(subgrp.model, digits = 2)))

        invisible(capture.output(summary(subgrp.model)))

        summ <- summarize.subgroups(subgrp.model)

        expect_is(summ, "data.frame")

        invisible(capture.output(print(summ, digits = 2, p.value = 0.25)))


        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt,
                                     propensity.func = propensity.multinom.lasso,
                                     reference.trt = 3,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")



        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt,
                                  propensity.func = propensity.multinom.lasso.bad,
                                  reference.trt = 3,
                                  loss   = "sq_loss_lasso",
                                  nfolds = 5) )

        # provide incorrect reference trt
        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt,
                                  propensity.func = propensity.multinom.lasso,
                                  reference.trt = "something_wrong",
                                  loss   = "sq_loss_lasso",
                                  nfolds = 5) )

        # no prop func
        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt,
                                     reference.trt = 3,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = as.factor(trt),
                                     propensity.func = propensity.multinom.lasso,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")


        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = as.factor(trt),
                                     reference.trt = "2",
                                     propensity.func = propensity.multinom.lasso,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")


        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = as.factor(trt),
                                     reference.trt = "2",
                                     larger.outcome.better = FALSE,
                                     propensity.func = propensity.multinom.lasso,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")
    }


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
        probs <- propens[cbind(1:nrow(propens), match(levels(trt)[trt], colnames(propens)))]

        probs
    }

    subgrp.model <- fit.subgroup(x = x, y = y,
                                 trt = trt,
                                 propensity.func = propensity.multinom.lasso,
                                 loss   = "sq_loss_lasso",
                                 nfolds = 5)              # option for cv.glmnet

    expect_is(subgrp.model, "subgroup_fitted")

    invisible(capture.output(print(subgrp.model, digits = 2)))

    invisible(capture.output(summary(subgrp.model)))

    summ <- summarize.subgroups(subgrp.model)

    expect_is(summ, "data.frame")

    invisible(capture.output(print(summ, digits = 2, p.value = 0.25)))


    if (Sys.info()[[1]] != "windows")
    {

        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt,
                                     propensity.func = propensity.multinom.lasso,
                                     reference.trt = 3,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")


        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = as.factor(trt),
                                     propensity.func = propensity.multinom.lasso,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")


        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = as.factor(trt),
                                     reference.trt = "2",
                                     propensity.func = propensity.multinom.lasso,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")


        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = as.factor(trt),
                                     reference.trt = "2",
                                     larger.outcome.better = FALSE,
                                     propensity.func = propensity.multinom.lasso,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")




        expect_error(fit.subgroup(x = x, y = y,
                                     trt = trt,
                                     propensity.func = prop.func,
                                     loss   = "sq_loss_lasso_gam",
                                     nfolds = 5))


        # different outcomes

        subgrp.model <- fit.subgroup(x = x, y = Surv(y.time.to.event, status),
                                     trt = as.factor(trt),
                                     propensity.func = propensity.multinom.lasso,
                                     loss   = "cox_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        subgrp.model <- fit.subgroup(x = x, y = y.binary,
                                     trt = as.factor(trt),
                                     propensity.func = propensity.multinom.lasso,
                                     loss   = "logistic_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        expect_error(fit.subgroup(x = x, y = Surv(y.time.to.event, status),
                                  trt = as.factor(trt),
                                  propensity.func = propensity.multinom.lasso,
                                  loss   = "sq_loss_lasso",
                                  nfolds = 5) )

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = as.factor(trt),
                                  propensity.func = propensity.multinom.lasso,
                                  loss   = "cox_loss_lasso",
                                  nfolds = 5) )


        propensity.multinom.lasso.array <- function(x, trt)
        {
            if (!is.factor(trt)) trt <- as.factor(trt)
            gfit <- cv.glmnet(y = trt, x = x, family = "multinomial")

            # predict returns a matrix of probabilities:
            # one column for each treatment level
            propens <- drop(predict(gfit, newx = x, type = "response", s = "lambda.min",
                                    nfolds = 5, alpha = 0))

            # return the probability corresponding to the
            # treatment that was observed
            probs <- array(dim = c(dim(propens), 2, 4))
            probs[is.na(probs)] <- 0.5
            probs
        }

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = as.factor(trt),
                                  propensity.func = propensity.multinom.lasso.array,
                                  loss   = "cox_loss_lasso",
                                  nfolds = 5) )

    }


})

