

context("check.overlap")

test_that("test plot is returned for hist/density/both", {


    set.seed(123)
    n.obs  <- 50
    n.vars <- 5
    x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


    # simulate non-randomized treatment
    xbetat   <- 0.25 + 0.5 * x[,1] - 0.5 * x[,5]
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

    # this one returns matrix with 1 column
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

    pl <- check.overlap(x = x,
                        trt = trt01,
                        propensity.func = prop.func,
                        type = "hist")

    expect_is(pl, "ggplot")


    pl <- check.overlap(x = x,
                        trt = trt01,
                        propensity.func = prop.func2,
                        type = "hist")

    expect_is(pl, "ggplot")

    pl <- check.overlap(x = x,
                        trt = trt01,
                        propensity.func = prop.func3,
                        type = "hist")

    expect_is(pl, "ggplot")

    pl <- check.overlap(x = x,
                        trt = trt01,
                        propensity.func = prop.func,
                        type = "density")

    expect_is(pl, "ggplot")

    pl <- check.overlap(x = x,
                        trt = trt01,
                        propensity.func = prop.func,
                        type = "both")

    expect_is(pl, "ggplot")


    # simulated non-randomized treatment with multiple levels
    xbetat_1   <- 0.15 + 0.5 * x[,1] - 0.25 * x[,5]
    xbetat_2   <- 0.15 - 0.5 * x[,2] + 0.25 * x[,3]
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

    pl <- check.overlap(x = x,
                        trt = trt,
                        type = "histogram",
                        propensity.func = propensity.multinom.lasso)

    expect_is(pl, "ggplot")


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

    pl <- check.overlap(x = x,
                        trt = trt,
                        type = "histogram",
                        propensity.func = propensity.multinom.lasso)

    expect_is(pl, "ggplot")


    # use multinomial logistic regression model with lasso penalty for propensity
    propensity.multinom.lasso.nonames <- function(x, trt)
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

        unname(probs)
    }

    pl <- check.overlap(x = x,
                        trt = trt,
                        type = "histogram",
                        propensity.func = propensity.multinom.lasso.nonames)

    expect_is(pl, "ggplot")

    pl <- check.overlap(x = x,
                        trt = as.factor(trt),
                        type = "histogram",
                        propensity.func = propensity.multinom.lasso.nonames)

    expect_is(pl, "ggplot")

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

    expect_error(check.overlap(x = x,
                               trt = trt,
                               type = "histogram",
                               propensity.func = propensity.multinom.lasso.array))
})
