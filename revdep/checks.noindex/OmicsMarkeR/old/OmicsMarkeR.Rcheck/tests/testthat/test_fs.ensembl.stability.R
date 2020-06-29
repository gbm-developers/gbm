library(OmicsMarkeR)
context("fs.ensembl.stability tests")

set.seed(123)
dat.discr <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50,
                             nsamp = 100,
                             st.dev = 1,
                             perturb = 0.2)),
    D = 10
)

vars <- dat.discr$discr.mat
groups <- dat.discr$classes
grp.levs <- levels(groups)

fit <- suppressWarnings(fs.ensembl.stability(vars, 
                                             groups, 
                                             method = c("plsda"), 
                                             f = 10,
                                             k = 3, 
                                             bags = 3,
                                             stability.metric = "canberra",
                                             k.folds = 3, 
                                             verbose = 'none'))

fit2 <- suppressWarnings(fs.ensembl.stability(vars, 
                                             groups, 
                                             method = c("glmnet", "rf"), 
                                             f = NULL,
                                             k = 3, 
                                             bags = 3,
                                             stability.metric = "canberra",
                                             k.folds = 3, 
                                             verbose = 'none'))


fit3 <- suppressWarnings(fs.ensembl.stability(vars, 
                                             groups, 
                                             method = c("svm", "plsda"), 
                                             f = 10,
                                             k = 3, 
                                             bags = 3,
                                             stability.metric = "canberra",
                                             k.folds = 3, 
                                             verbose = 'none'))


test_that("Output is correct", {
    expect_is(fit, "list")
    expect_is(fit2, "list")
    expect_is(fit3, "list")
})

