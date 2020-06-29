library(OmicsMarkeR)
context("fs.stability checks")

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
multi_grps <- as.factor(c(rep(c("A","B","C"), each=33), "C"))

fit <- suppressWarnings(fs.stability(vars,
                                     groups,
                                     method = c("plsda"),
                                     f = 10,
                                     k = 3,
                                     k.folds = 3,
                                     verbose = 'none'))

fit_multi <- suppressWarnings(fs.stability(vars,
                                           groups,
                                           method = c("rf", "svm"),
                                           f = 10,
                                           k = 3,
                                           k.folds = 3,
                                           verbose = 'none'))

fit_multi_grp <- suppressWarnings(fs.stability(vars,
                                               multi_grps,
                                               method = c("svm"),
                                               optimize.resample = TRUE,
                                               f = 10,
                                               k = 3,
                                               k.folds = 3,
                                               verbose = 'none'))

fit_fs <- suppressWarnings(fs.stability(vars,
                                        groups,
                                        method = c("svm", "rf", 
                                                   "glmnet", "plsda"),
                                        model.features = TRUE,
                                        k = 3,
                                        k.folds = 3,
                                        verbose = 'none'))


test_that("Inputs are correct", {
    expect_error(fs.stability(X = c(rnorm(100)),
                              Y = groups,
                              method = c("plsda"),
                              f = 10,
                              k = 3,
                              k.folds = 3,
                              verbose='none'),
                 info = "fs.stability doesn't accept vectors")
    
    
    expect_error(fs.stability(X = vars,
                              Y = groups,
                              method = c("plsda"),
                              model.features = TRUE,
                              k = 3,
                              f = 10,
                              k.folds = 3,
                              stability.metric = "kuncheva",
                              verbose='none'),
                 info = "kuncheva cannot be used with model.features=TRUE")
    
    expect_error(fs.stability(X = matrix(as.character(rnorm(5000)), ncol=50),
                              Y = groups,
                              method = c("plsda"),
                              f = 10,
                              k = 3,
                              k.folds = 3,
                              verbose='none'),
                 info = "fs.stability doesn't accept character matrix")
    
    expect_error(fs.stability(X = as.data.frame(vars),
                              Y = groups,
                              method = c("plsda"),
                              f = 10,
                              k = 3,
                              k.folds = 3,
                              verbose='none'),
                 info = "fs.stability doesn't accept data.frame")
    
    expect_error(fs.stability(X = vars,
                              Y = rnorm(100),
                              method = c("plsda"),
                              f = 10,
                              k = 3,
                              k.folds = 3,
                              verbose='none'),
                 info = "Y cannot be numeric")
    
    expect_error(fs.stability(X = vars,
                              Y = groups,
                              method = c("myfun"),
                              f = 10,
                              k = 3,
                              k.folds = 3,
                              verbose='none'),
                 info = "Custom/Non-implemented methods
                 return an error")
    
    expect_error(fs.stability(X = vars,
                              Y = groups,
                              method = c("plsda"),
                              stability.metric = "kuncheva",
                              k = 3,
                              f = NULL,
                              k.folds = 3,
                              verbose='none'),
                 info = "Kuncheva requires user to specify
                 the number of features expected")
    
    expect_error(fs.stability(X = vars,
                              Y = groups,
                              method = c("plsda"),
                              stability.metric = "kuncheva",
                              k = 3,
                              k.folds = 3,
                              verbose='none'),
                 info = "Kuncheva requires user to specify
                 the number of features expected")
    
    expect_error(fs.stability(X = vars,
                              Y = groups,
                              method = c("plsda"),
                              k = 3,
                              f = 10,
                              k.folds = 3,
                              model.features = TRUE,
                              verbose='none'),
                 info = "You can't specify f when model.features
                 is set to TRUE")
    
    expect_warning(fs.stability(X = vars,
                                Y = groups,
                                method = c("plsda"),
                                stability.metric = "spearman",
                                k = 3,
                                f = 10,
                                k.folds = 3,
                                verbose='none'),
                   info = "Rank correlation sets f to NULL")
    
    expect_error(fs.stability(X = vars,
                              Y = groups,
                              method = c("plsda"),
                              k = 3,
                              f = NULL,
                              k.folds = 3,
                              verbose='none'),
                 info = "Non-rank stability method requires
                 either model.features=TRUE or 'f' to be set")
    
    expect_warning(fs.stability(X = vars,
                              Y = groups,
                              method = c("plsda"),
                              f = 10,
                              k = 1,
                              k.folds = 3,
                              verbose='none'),
                 info = "No replicates were run, stability will be set to NULL.
                    If stability not desired, see 'fit.only.model'")
})

test_that("Outputs are correct", {
    expect_is(fit, "list")
    expect_is(fit_multi_grp, "list")
    expect_is(fit$methods, "character")
    expect_is(fit$performance, "list")
    expect_is(fit$RPT, "numeric")
    expect_is(fit$features, "list")
    expect_is(fit$specs, "list")
    
    # values taken from start of test file
    expect_true(fit$specs$total.samples == 100)
    expect_true(fit$specs$number.features == 50)
    expect_true(all(fit$specs$group.levels == grp.levs))
    expect_true(all(fit$specs$number.observations.group == table(groups)))
    
    expect_true(!is.null(fit_multi_grp$original.best.tunes))
    expect_true(!is.null(fit_multi_grp$final.best.tunes))
})

test_that("feature.table works correctly", {
    ftab <- feature.table(fit, "plsda")
    expect_is(ftab, "data.frame")
    expect_error(feature.table(fit, "rf"))
    
})

test_that("performance.metrics works correctly", {
    perfs <- performance.metrics(fit)
    perfs_multi <- performance.metrics(fit_multi)
    expect_is(perfs, "data.frame")
    expect_is(perfs_multi, "data.frame")
    expect_true(all(perfs >= 0 & perfs <= 1))
    expect_true(all(perfs_multi >= 0 & perfs_multi <= 1))
})

test_that("predictNewClasses works correctly", {
    newdata <- create.discr.matrix(
        create.corr.matrix(
            create.random.matrix(nvar = 50,
                                 nsamp = 100,
                                 st.dev = 1,
                                 perturb = 0.2)),
        D = 10
    )$discr.mat
    
    orig.df <- data.frame(vars, groups)
    
    newClasses <- lapply(c("plsda","glmnet","rf","svm"), function(x) {
        suppressWarnings(
            predictNewClasses(fit_fs, x, orig.df, newdata)
        )
    })
    
    expect_true(all(sapply(newClasses, is.data.frame)))
    expect_true(all(sapply(newClasses, function(x){
        is.factor(x$predictedClass)
    })))
})

test_that("perm.features inputs are correct", {
    expect_error(
        perm.features(fit_fs, vars, groups, "plsda", 
                      sig.level=0.5, nperm=100001,
                      allowParallel=FALSE),
        info = "nperm exceeds currently threshold of 100,000")
    expect_error(
        perm.features(fit_fs, vars, groups, "myfun", 
                      sig.level=0.5, nperm=10,
                      allowParallel=FALSE),
        info = "must provide either a fit_fsted model from
        fs.stability or parameter(s) for desired model")
})

test_that("perm.class inputs are correct", {
    expect_error(
        perm.class(fit_fs, vars, groups, "myfun", k.folds=5,
                   metric="Accuracy", nperm=10, verbose=FALSE),
        info = "must provide either a fit_fsted model from
        fs.stability or one provided in 'modelList()'"
    )
    expect_error(
        perm.class(fs.model = NULL, vars, groups, "plsda", k.folds=5,
                   metric="Accuracy", nperm=10, verbose=FALSE),
        info = "must provide either a fit_fsted model from
        fs.stability or parameter(s) for desired model"
    )
})

test_that("perm.features outputs are correct", {
    featPerms_plsda <- suppressWarnings(perm.features(fit_fs, vars, groups, "plsda",
                                                      sig.level=0.5, nperm=10,
                                                      allowParallel=FALSE, verbose=FALSE))
    
    featPerms_svm <- suppressWarnings(perm.features(fit_fs, vars, groups, "svm",
                                                    sig.level=0.5, nperm=10,
                                                    allowParallel=FALSE, verbose=FALSE))
    
    featPerms_rf <- suppressWarnings(perm.features(fit_fs, vars, groups, "rf",
                                                   sig.level=0.5, nperm=10,
                                                   allowParallel=FALSE, verbose=FALSE))
    
    featPerms_glmnet <- suppressWarnings(perm.features(fit_fs, vars, groups, "glmnet",
                                                       sig.level=0.5, nperm=10,
                                                       allowParallel=FALSE, verbose=FALSE))
    
    expect_is(featPerms_plsda, "list") 
    expect_is(featPerms_rf, "list") 
    expect_is(featPerms_svm, "list") 
    expect_is(featPerms_glmnet, "list") 
})

test_that("perm.class outputs are correct", {    
    classPerms_plsda <- suppressWarnings(
        perm.class(fit_fs, vars, groups, "plsda", k.folds=5,
                   metric="Accuracy", nperm=10, verbose=FALSE))
    
    classPerms_rf <- suppressWarnings(
        perm.class(fit_fs, vars, groups, "rf", k.folds=5,
                   metric="Accuracy", nperm=10, verbose=FALSE))
    
    classPerms_glmnet <- suppressWarnings(
        perm.class(fit_fs, vars, groups, "glmnet", k.folds=5,
                   metric="Accuracy", nperm=10, verbose=FALSE))
    
    classPerms_svm <- suppressWarnings(
        perm.class(fit_fs, vars, groups, "svm", k.folds=5,
                   metric="Accuracy", nperm=10, verbose=FALSE))
    
    expect_is(classPerms_plsda, "data.frame")
    expect_is(classPerms_rf, "data.frame")
    expect_is(classPerms_svm, "data.frame")
    expect_is(classPerms_glmnet, "data.frame")
})

