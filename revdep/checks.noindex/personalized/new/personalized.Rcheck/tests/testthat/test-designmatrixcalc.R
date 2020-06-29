
context("design matrix checking")

test_that("test design matrix for binary trt works - weighting", {

    x <- matrix(1:8, ncol = 2)
    y <- 1:4
    trt <- c(rep(1, 2), rep(0, 2))
    pi.x <- seq(0.2, 0.8, length.out = 4)

    xtilde <- create.design.matrix.binary.trt(x, pi.x, trt, "weighting")

    expect_equivalent(xtilde, (2 * trt - 1)  * cbind(1, x) )
})


test_that("test create.design.matrix for binary trt works - weighting", {

    x <- matrix(1:8, ncol = 2)
    y <- 1:4
    trt <- c(rep(1, 2), rep(0, 2))
    pi.x <- seq(0.2, 0.8, length.out = 4)

    xtilde <- create.design.matrix(x, pi.x, trt, y, "weighting")

    expect_equivalent(xtilde, (2 * trt - 1)  * cbind(1, x) )
})

test_that("test design matrix for multi trt works - weighting", {

    x <- matrix(1:10, ncol = 2)
    y <- 1:5
    trt <- as.factor(c(1, 1, 2, rep(0, 2)))
    pi.x <- seq(0.2, 0.8, length.out = 5)

    xtilde <- create.design.matrix.mult.trt(x, pi.x, trt, "weighting", reference.trt = "2")

    xm <- matrix(0, ncol=4, nrow=5)
    xm[1:2,3:4] <- x[1:2,]
    xm[4:5,1:2] <- x[4:5,]
    xm[3,] <- rep(-x[3,], 2)

    expect_equivalent(xtilde, xm )


    xtilde <- create.design.matrix.mult.trt(x, pi.x, trt, "weighting")

    expect_is(xtilde, "matrix")
})


test_that("test create.design.matrix for multi trt works - weighting", {

    x <- matrix(1:10, ncol = 2)
    y <- 1:5
    trt <- c(1, 1, 2, rep(0, 2))
    pi.x <- seq(0.2, 0.8, length.out = 5)

    xtilde <- create.design.matrix(x, pi.x, trt, y, "weighting", reference.trt = "2")

    x2 <- cbind(1, x)
    xm <- matrix(0, ncol=4 + 2, nrow=5)
    xm[1:2,4:6] <- x2[1:2,]
    xm[4:5,1:3] <- x2[4:5,]
    xm[3,] <- rep(-x2[3,], 2)

    expect_equivalent(xtilde, xm )

    expect_error(xtilde <- create.design.matrix(x, pi.x, trt, y, "a_learning", reference.trt = "2"))
})

test_that("test design matrix for binary trt works - a learning", {

    x <- matrix(1:8, ncol = 2)
    y <- 1:4
    trt <- c(rep(1, 2), rep(0, 2))
    pi.x <- seq(0.2, 0.8, length.out = 4)

    xtilde <- create.design.matrix.binary.trt(x, pi.x, trt, "a_learning")

    expect_equivalent(xtilde, (1 * (trt != 0) - pi.x)   * cbind(1, x) )
})


test_that("test create.design.matrix for binary trt works - a learning", {

    x <- matrix(1:8, ncol = 2)
    y <- 1:4
    trt <- c(rep(1, 2), rep(0, 2))
    pi.x <- seq(0.2, 0.8, length.out = 4)

    xtilde <- create.design.matrix(x, pi.x, trt, y, "a_learning")

    expect_equivalent(xtilde, (1 * (trt != 0) - pi.x)   * cbind(1, x) )
})

test_that("test weights for binary trt works - weighting", {

    trt <- c(rep(1, 2), rep(0, 2))
    pi.x <- seq(0.2, 0.8, length.out = 4)

    wts <- create.weights.binary.trt(pi.x = pi.x, trt = trt, "weighting")

    expect_equivalent(wts, 1 / (pi.x * (trt == 1) + (1 - pi.x) * (trt == 0)) )
})

test_that("test weights for binary trt works - a_learning", {

    trt <- c(rep(1, 2), rep(0, 2))
    pi.x <- seq(0.2, 0.8, length.out = 4)

    wts <- create.weights.binary.trt(pi.x = pi.x, trt = trt, "a_learning")

    expect_equivalent(wts, rep(1, length(trt)) )
})



test_that("test weights for multi trt works - weighting", {

    trt <- c(1:4)
    pi.x <- seq(0.2, 0.8, length.out = 4)

    wts <- create.weights.mult.trt(pi.x = pi.x, trt = trt, "weighting")

    expect_equivalent(wts, 1 / (pi.x) )
})

test_that("test weights for multi trt works - a_learning", {

    trt <- c(1:4)
    pi.x <- seq(0.2, 0.8, length.out = 4)

    wts <- create.weights.mult.trt(pi.x = pi.x, trt = trt, "a_learning")

    expect_equivalent(wts, rep(1, length(trt)) )
})


test_that("test create.weights for multi trt works - weighting", {

    trt <- c(1:4)
    pi.x <- seq(0.2, 0.8, length.out = 4)

    wts <- create.weights(pi.x = pi.x, trt = trt, "weighting")

    expect_equivalent(wts, 1 / (pi.x) )
})

test_that("test create.weights for multi trt works - a_learning", {

    trt <- c(1:4)
    pi.x <- seq(0.2, 0.8, length.out = 4)

    wts <- create.weights(pi.x = pi.x, trt = trt, "a_learning")

    expect_equivalent(wts, rep(1, length(trt)) )
})


test_that("test weights for binary trt works - weighting", {

    trt <- c(rep(1, 2), rep(0, 2))
    pi.x <- seq(0.2, 0.8, length.out = 4)

    wts <- create.weights(pi.x = pi.x, trt = trt, "weighting")

    expect_equivalent(wts, 1 / (pi.x * (trt == 1) + (1 - pi.x) * (trt == 0)) )
})

test_that("test weights for binary trt works - a_learning", {

    trt <- c(rep(1, 2), rep(0, 2))
    pi.x <- seq(0.2, 0.8, length.out = 4)

    wts <- create.weights(pi.x = pi.x, trt = trt, "a_learning")

    expect_equivalent(wts, rep(1, length(trt)) )
})
