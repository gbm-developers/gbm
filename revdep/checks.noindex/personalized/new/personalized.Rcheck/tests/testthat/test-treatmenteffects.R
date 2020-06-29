
context("treatment effect calculations")

test_that("test that treatment effect calculations work", {

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

    subgrp.model <- fit.subgroup(x = x, y = y,
                                 trt = trt01,
                                 propensity.func = prop.func,
                                 loss   = "sq_loss_lasso",
                                 nfolds = 5)    # option for cv.glmnet

    trt_eff <- treatment.effects(subgrp.model)

    print(trt_eff)

    expect_true(is.na(trt_eff$gamma))

    subgrp.modela <- fit.subgroup(x = x, y = y,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  loss   = "sq_loss_lasso",
                                  method = "a_learning",
                                  nfolds = 5)    # option for cv.glmnet

    trt_eff <- treatment.effects(subgrp.modela)
    expect_true(is.na(trt_eff$gamma))

    print(trt_eff)


    trt_eff <- treat.effects(subgrp.modela$benefit.scores,
                             subgrp.modela$loss,
                             subgrp.modela$method,
                             subgrp.modela$pi.x)

    expect_error(treat.effects(subgrp.modela$benefit.scores,
                               "poisson_loss_lasso_gam",
                               subgrp.modela$method))

    print(trt_eff)


    expect_warning(treat.effects(subgrp.modela$benefit.scores,
                                 "owl_logistic_flip_loss_gam",
                                 subgrp.modela$method,
                                 subgrp.modela$pi.x))


    library(survival)
    subgrp.model.cox <- fit.subgroup(x = x, y = Surv(y.time.to.event, status),
                                     trt = trt01,
                                     propensity.func = prop.func,
                                     loss   = "cox_loss_lasso",
                                     nfolds = 5)              # option for cv.glmnet

    trt_eff_c <- treatment.effects(subgrp.model.cox)

    expect_true(all(trt_eff_c$gamma >= 0))
    expect_true(is.na(trt_eff_c$delta))

    print(trt_eff_c)


    ## other calculation types

    subgrp.model <- fit.subgroup(x = x, y = y.binary,
                                 trt = trt01,
                                 propensity.func = prop.func,
                                 loss   = "owl_logistic_loss_lasso",
                                 nfolds = 5)    # option for cv.glmnet

    trt_eff <- treatment.effects(subgrp.model)

    expect_true(all(trt_eff$gamma >= 0))
    expect_true(is.na(trt_eff$delta))



    subgrp.model <- fit.subgroup(x = x, y = y.count,
                                 trt = trt01,
                                 propensity.func = prop.func,
                                 loss   = "poisson_loss_lasso",
                                 nfolds = 5)    # option for cv.glmnet

    trt_eff <- treatment.effects(subgrp.model)

    subgrp.modela <- fit.subgroup(x = x, y = y.count,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  loss   = "poisson_loss_lasso",
                                  method = "a_learning",
                                  nfolds = 5)    # option for cv.glmnet

    trt_eff <- treatment.effects(subgrp.modela)



    subgrp.model <- fit.subgroup(x = x, y = y.binary,
                                 trt = trt01,
                                 propensity.func = prop.func,
                                 loss   = "logistic_loss_lasso",
                                 nfolds = 5)    # option for cv.glmnet

    trt_eff <- treatment.effects(subgrp.model)

    subgrp.modela <- fit.subgroup(x = x, y = y.binary,
                                  trt = trt01,
                                  propensity.func = prop.func,
                                  loss   = "logistic_loss_lasso",
                                  method = "a_learning",
                                  nfolds = 5)    # option for cv.glmnet

    trt_eff <- treatment.effects(subgrp.modela)



})
