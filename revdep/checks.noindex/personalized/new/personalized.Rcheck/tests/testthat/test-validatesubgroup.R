

context("validate.subgroup")

test_that("test validate.subgroup for continuous outcomes with various options", {
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
    delta <- 5 * (0.5 + x[,2] - x[,3]  )
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

    subgrp.model <- fit.subgroup(x = cbind(x, rbinom(n.obs, 1, 0.25)), y = y,
                                 trt = trt01,
                                 propensity.func = prop.func,
                                 loss   = "sq_loss_lasso",
                                 nfolds = 5)              # option for cv.glmnet

    expect_is(subgrp.model, "subgroup_fitted")

    invisible(capture.output(print(subgrp.model, digits = 2)))

    invisible(capture.output(summary(subgrp.model)))

    summ <- summarize.subgroups(subgrp.model)

    expect_is(summ, "data.frame")

    invisible(capture.output(print(summarize.subgroups(subgrp.model),
                                   digits = 2, p.value = 0.25)))

    subgrp.val <- validate.subgroup(subgrp.model, B = 3,
                                    method = "training")

    expect_is(subgrp.val, "subgroup_validated")

    invisible(capture.output(print(subgrp.val, digits = 2)))

    invisible(capture.output(print(subgrp.val, digits = 2, sample.pct = TRUE)))

    subgrp.val <- validate.subgroup(subgrp.model, B = 3,
                                    method = "boot")

    expect_is(subgrp.val, "subgroup_validated")

    invisible(capture.output(print(subgrp.val, digits = 2)))


    if (Sys.info()[[1]] != "windows")
    {
        expect_error(validate.subgroup(x, B = 3,
                                       method = "training"))

        ## parallel

        subgrp.val <- validate.subgroup(subgrp.model, B = 3,
                                        parallel = TRUE,
                                        method = "training")

        expect_is(subgrp.val, "subgroup_validated")


        subgrp.val <- validate.subgroup(subgrp.model, B = 3,
                                        parallel = TRUE,
                                        method = "boot")

        expect_is(subgrp.val, "subgroup_validated")
    }
})




test_that("test validate.subgroup for binary outcomes and various losses", {
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
    delta <- 5 * (0.5 + x[,2] - x[,3]  )
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
        pmin(pmax(pi.x, 1e-5), 1 - 1e-05)
    }

    subgrp.model <- fit.subgroup(x = x, y = y.binary,
                                 trt = trt01,
                                 propensity.func = prop.func,
                                 loss   = "logistic_loss_lasso",
                                 nfolds = 5)              # option for cv.glmnet

    expect_is(subgrp.model, "subgroup_fitted")

    invisible(capture.output(print(subgrp.model, digits = 2)))

    invisible(capture.output(summary(subgrp.model)))

    summ <- summarize.subgroups(subgrp.model)

    expect_is(summ, "data.frame")

    invisible(capture.output(print(summarize.subgroups(subgrp.model), digits = 2, p.value = 0.25)))

    if (Sys.info()[[1]] != "windows")
    {

        subgrp.model2 <- fit.subgroup(x = x, y = y,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     larger.outcome.better = FALSE,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        subgrp.val <- validate.subgroup(subgrp.model, B = 3,
                                        benefit.score.quantiles = NULL,
                                        method = "training")



        subgrp.val <- validate.subgroup(subgrp.model, B = 3,
                                        benefit.score.quantiles = numeric(0),
                                        method = "training")

        subgrp.val <- validate.subgroup(subgrp.model, B = 3,
                                        method = "training")

        subgrp.val2 <- validate.subgroup(subgrp.model2, B = 3,
                                        method = "training")

        print(subgrp.val)
        print(subgrp.val2)

        expect_error(validate.subgroup(subgrp.val, B = 3, method = "training"))

        expect_error(validate.subgroup(subgrp.model, B = 3, train.fraction = -1,
                                       method = "training"))

        expect_error(validate.subgroup(subgrp.model, B = 3, train.fraction = 2,
                                       method = "training"))

        expect_error(print(subgrp.val, which.quant = 99))

        expect_error(print(subgrp.val, which.quant = 1:10))

        print(subgrp.val, which.quant = c(4, 5))

        print(subgrp.val, which.quant = c(4, 5), sample.pct = TRUE)


        subgrp.model2 <- fit.subgroup(x = x, y = y.binary,
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "logistic_loss_lasso",
                                     retcall = FALSE,
                                     nfolds = 5)

        # retcall must be true
        expect_error(validate.subgroup(subgrp.model2, B = 3,
                                       method = "training"))

        expect_is(subgrp.val, "subgroup_validated")

        invisible(capture.output(print(subgrp.val, digits = 2)))


        subgrp.val <- validate.subgroup(subgrp.model, B = 3,
                                        method = "train")

        expect_is(subgrp.val, "subgroup_validated")

        subgrp.val <- validate.subgroup(subgrp.model, B = 3,
                                        method = "boot")

        expect_is(subgrp.val, "subgroup_validated")

        invisible(capture.output(print(subgrp.val, digits = 2)))




        subgrp.model <- fit.subgroup(x = x, y = Surv(y.time.to.event, status),
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "cox_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        subgrp.val <- validate.subgroup(subgrp.model, B = 3,
                                        method = "train")

        expect_is(subgrp.val, "subgroup_validated")

        print(subgrp.val)


        subgrp.val <- validate.subgroup(subgrp.model, B = 3,
                                        method = "boot")

        expect_is(subgrp.val, "subgroup_validated")

        print(subgrp.val)
    }

    ####### mult trts


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

    if (Sys.info()[[1]] != "windows")
    {
        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt,
                                     propensity.func = propensity.multinom.lasso,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 3)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        print(subgrp.model)

        subgrp.val <- validate.subgroup(subgrp.model, B = 2,
                                        method = "train")

        expect_is(subgrp.val, "subgroup_validated")

        print(subgrp.val)

        print(subgrp.val, which.quant = c(2,4))


        subgrp.model <- fit.subgroup(x = x, y = y,
                                     trt = trt,
                                     propensity.func = propensity.multinom.lasso,
                                     larger.outcome.better = FALSE,
                                     loss   = "sq_loss_lasso",
                                     nfolds = 3)              # option for cv.glmnet

        expect_is(subgrp.model, "subgroup_fitted")

        print(subgrp.model)

        subgrp.val <- validate.subgroup(subgrp.model, B = 2,
                                        method = "train")

        expect_is(subgrp.val, "subgroup_validated")

        print(subgrp.val)

        print(subgrp.val, which.quant = c(2,4))
    }

})
