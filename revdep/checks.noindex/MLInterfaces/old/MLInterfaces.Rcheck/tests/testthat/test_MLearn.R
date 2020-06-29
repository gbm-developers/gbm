test_that("svmI: y as factor/char", {
    e <- matrix(rnorm(100),
                ncol = 10)
    colnames(e) <- rownames(e) <- LETTERS[1:10]
    e <- ExpressionSet(e)
    e$k <- LETTERS[1:10]
    ## This fails with R-devel (4.0) prior version 1.67.10 with:
    ## Error in svm.default(x, y, scale = scale, ..., na.action
    ## = na.action) : Need numeric dependent variable for regression.
    cl1 <- MLearn(k ~ ., e, svmI, 1:5)
    e$k <- factor(e$k)
    cl2 <- MLearn(k ~ ., e, svmI, 1:5)
    expect_identical(cl1, cl2)
})
