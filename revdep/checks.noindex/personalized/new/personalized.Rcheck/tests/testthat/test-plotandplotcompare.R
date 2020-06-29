

context("plot.subgroup_fitted and plotCompare and plot.subgroup_validated")

test_that("test plotting for continuous outcomes with various options", {
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
    delta <- 4 * (0.5 + x[,2] - x[,3]  )
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

    subgrp.model <- fit.subgroup(x = x, y = y,
                                 trt = trt01,
                                 propensity.func = prop.func,
                                 loss   = "sq_loss_lasso",
                                 cutpoint = "median",
                                 nfolds = 5)              # option for cv.glmnet

    expect_is(subgrp.model, "subgroup_fitted")

    if (Sys.info()[[1]] != "windows")
    {
        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  loss   = "sq_loss_lasso",
                                  cutpoint = "haha",
                                  nfolds = 5))

        expect_error(fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  loss   = "sq_loss_lasso",
                                  cutpoint = list(x),
                                  nfolds = 5))

        subgrp.model.bin <- fit.subgroup(x = x, y = y.binary,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "logistic_loss_lasso",
                                     cutpoint = "median",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model.bin, "subgroup_fitted")

        subgrp.model.bin2 <- fit.subgroup(x = x, y = y.binary,
                                         trt = trt01,
                                         propensity.func = prop.func,
                                         loss   = "logistic_loss_lasso",
                                         cutpoint = 0,
                                         nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model.bin2, "subgroup_fitted")

        plot(subgrp.model.bin, type = "boxplot")

        pl <- plotCompare(subgrp.model.bin, subgrp.model.bin2, type = "density")

        expect_is(pl, "ggplot")

        pl <- plotCompare(subgrp.model.bin, subgrp.model.bin2, type = "boxplot")

        expect_is(pl, "ggplot")

    }

    subgrp.model2 <- fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  loss   = "owl_logistic_loss_lasso",
                                  cutpoint = "quant50",
                                  nfolds = 5)              # option for cv.glmnet

    expect_is(subgrp.model2, "subgroup_fitted")

    subgrp.val <- validate.subgroup(subgrp.model, B = 3,
                                    method = "training")

    subgrp.val2 <- validate.subgroup(subgrp.model2, B = 3,
                                    method = "training")

    expect_is(subgrp.val, "subgroup_validated")

    # plot fitted
    pl <- plot(subgrp.model, type = "boxplot")

    expect_is(pl, "ggplot")

    pl <- plot(subgrp.model, type = "density")

    expect_is(pl, "ggplot")

    pl <- plot(subgrp.model, type = "density", avg.line = FALSE)

    expect_is(pl, "ggplot")

    pl <- plot(subgrp.model, type = "interaction")

    expect_is(pl, "ggplot")

    pl <- plot(subgrp.model, type = "conditional")

    expect_is(pl, "ggplot")


    # plot validated

    pl <- plot(subgrp.val, type = "boxplot")

    expect_is(pl, "ggplot")

    pl <- plot(subgrp.val, type = "density")

    expect_is(pl, "ggplot")

    pl <- plot(subgrp.val, type = "density", avg.line = FALSE)

    expect_is(pl, "ggplot")

    pl <- plot(subgrp.val, type = "interaction")

    expect_is(pl, "ggplot")

    pl <- plot(subgrp.val, type = "conditional")

    expect_is(pl, "ggplot")

    pl <- plot(subgrp.val, type = "stability")

    expect_is(pl, "plotly")

    # plot compare

    pl <- plotCompare(subgrp.model, subgrp.model2, type = "boxplot")

    expect_is(pl, "ggplot")

    pl <- plotCompare(subgrp.model, subgrp.model2, type = "density")

    expect_is(pl, "ggplot")

    pl <- plotCompare(subgrp.model, subgrp.model2, type = "conditional")

    expect_is(pl, "ggplot")

    expect_error(plotCompare(subgrp.model, subgrp.val, type = "conditional"))

    pl <- plotCompare(subgrp.model, subgrp.model2, type = "density", avg.line = FALSE)

    expect_is(pl, "ggplot")


    pl <- plotCompare(subgrp.model, subgrp.model2, type = "interaction")

    expect_is(pl, "ggplot")




    pl <- plotCompare(subgrp.model, subgrp.val, type = "boxplot")

    expect_is(pl, "ggplot")

    pl <- plotCompare(subgrp.val, subgrp.val2, type = "conditional")

    expect_is(pl, "ggplot")

    pl <- plotCompare(subgrp.model, subgrp.val, type = "density")

    expect_is(pl, "ggplot")

    pl <- plotCompare(subgrp.model, subgrp.val, type = "density", avg.line = FALSE)

    expect_is(pl, "ggplot")


    pl <- plotCompare(subgrp.model, subgrp.val, type = "interaction")

    expect_is(pl, "ggplot")

    subgrp.model2 <- subgrp.model
    pl <- plotCompare(subgrp.model, subgrp.model2, type = "conditional")

    expect_is(pl, "ggplot")

    y2 <- 1 * (y > 0)

    subgrp.model <- fit.subgroup(x = x, y = y2,
                                 trt = trt01,
                                 propensity.func = prop.func,
                                 loss   = "logistic_loss_lasso",
                                 nfolds = 5)              # option for cv.glmnet

    expect_is(subgrp.model, "subgroup_fitted")

    pl <- plot(subgrp.model, type = "conditional")

    expect_is(pl, "ggplot")

    subgrp.model2 <- subgrp.model
    pl <- plotCompare(subgrp.model, subgrp.model2, type = "conditional")

    expect_is(pl, "ggplot")



    subgrp.model <- fit.subgroup(x = x, y = Surv(y.time.to.event, status),
                                 trt = trt01,
                                 propensity.func = prop.func,
                                 loss   = "cox_loss_lasso",
                                 nfolds = 5)              # option for cv.glmnet

    expect_is(subgrp.model, "subgroup_fitted")

    subgrp.val <- validate.subgroup(subgrp.model, B = 3,
                                    method = "train")

    expect_is(subgrp.val, "subgroup_validated")


    pl <- plot(subgrp.model, type = "boxplot")

    expect_is(pl, "ggplot")

    pl <- plot(subgrp.val, type = "boxplot")

    expect_is(pl, "ggplot")


    subgrp.val2 <- validate.subgroup(subgrp.model, B = 3,
                                    method = "boot")

    expect_is(subgrp.val2, "subgroup_validated")

    pl <- plot(subgrp.val2, type = "boxplot")

    expect_is(pl, "ggplot")

    subgrp.model2 <- subgrp.model
    pl <- plotCompare(subgrp.model, subgrp.model2, type = "conditional")

    expect_is(pl, "ggplot")




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

})



