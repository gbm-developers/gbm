library(OmicsMarkeR)
context("Stability Metrics")

set.seed(123)
x <- LETTERS[1:10]
y <- LETTERS[11:20]
a <- seq(10)
b <- rev(a)

test_that("stability functions' inputs are correct", {
    expect_error(spearman(as.factor(seq(10)), as.factor(seq(10))))
    expect_error(spearman(sample(seq(10)), sample(seq(20))))
    expect_error(canberra(as.factor(seq(10)), as.factor(seq(10))))
    expect_error(jaccard(sample(seq(10)), sample(seq(10))),
                 info = "jaccard requires character inputs")
    expect_error(pof(sample(seq(10)), sample(seq(10))),
                 info = "pof requires character inputs")
    expect_error(sorensen(sample(seq(10)), sample(seq(10))),
                 info = "sorensen requires character inputs")
    expect_error(ochiai(sample(seq(10)), sample(seq(10))),
                 info = "ochiai requires character inputs")
    expect_error(kuncheva(sample(seq(10)), sample(seq(10)), 20),
                 info = "kunchvea requires character inputs")
    expect_error(kuncheva(LETTERS[1:10], LETTERS[5:12], 20),
                 info = "kunchvea requires feature subset to have
                 the same cardinality")
    expect_error(kuncheva(LETTERS[1:10], LETTERS[1:10], "20"),
                 info = "nc must be numeric")
    expect_error(kuncheva(as.character(sample(seq(10))), 
                          as.character(sample(seq(10)))),
                 info = "kuncheva requires number.features to be defined")
})

test_that("stability functions' outputs are correct", {
    jac_max <- jaccard(x,x)
    jac_min <- jaccard(x,y)
    pov_max <- pof(x,x)
    pov_min <- pof(x,y)
    sor_max <- sorensen(x,x)
    sor_min <- sorensen(x,y)
    och_max <- ochiai(x,x)
    och_min <- ochiai(x,y)
    kun_max <- kuncheva(x,x,20)
    kun_min <- kuncheva(x,y,20)
    can_max <- canberra_stability(a,a)
    can_min <- canberra_stability(a,b)
    
    expect_is(spearman(sample(seq(10)), sample(seq(10))), "numeric")
    expect_is(canberra(sample(seq(10)), sample(seq(10))), "numeric")
    expect_is(jac_max, "numeric")
    expect_true(jac_min == 0,
                info = "jaccard minimum value is 0")
    expect_true(jac_max == 1, 
                info = "jaccard maximum value is 1")
    expect_is(pov_max, "numeric")
    expect_true(pov_min == 0,
                info = "pof minimum value is 0")
    expect_true(pov_max == 1, 
                info = "pof maximum vlue is 1")
    expect_is(sor_max, "numeric")
    expect_true(sor_min == 0,
                info = "sorensen minimum value is 0")
    expect_true(sor_max == 1, 
                info = "sorensen maximum value 1")
    expect_is(och_max, "numeric")
    expect_true(och_min == 0,
                info = "ochiai minimum value is 0")
    expect_true(och_max == 1, 
                info = "ochiai maximum value is 1")
    expect_is(kun_max, "numeric")
    expect_true(kun_min == 0,
                info = "kuncheva minimum value is 0")
    expect_true(kun_max == 1, 
                info = "kuncheva maximum value is 1")
    expect_is(can_max, "numeric")
    expect_true(can_min == 0,
                info = "canberra stability minimum value is 0")
    expect_true(can_max == 1, 
                info = "canberra stability maximum value is 1")
    
})

test_that("pairwise.stability requires appropriate inputs", {
    some.numbers <- seq(20)
    
    # matrix of Metabolites identified (e.g. 5 trials)
    features <- 
        replicate(5, paste("Metabolite", sample(some.numbers, 10), sep="_"))
    
    expect_error(pairwise.stability(features, "kuncheva"), 
                 'argument "nc" is missing, with no default')
})

test_that("pairwise.model.stability requires appropriate inputs", {
    some.numbers <- seq(20)
    
    plsda <- 
        replicate(5, paste("Metabolite", sample(some.numbers, 10), sep="_"))
    rf <-
        replicate(5, paste("Metabolite", sample(some.numbers, 10), sep="_"))
    
    features <- list(plsda=plsda, rf=rf)
    
    expect_error(pairwise.model.stability(features, "kuncheva"), 
                 'argument "nc" is missing, with no default')
})

test_that("RPT requires appropriate inputs", {
    expect_error(RPT("A", 0.80), "stability is not of class 'numeric'")
    expect_error(RPT(0.80, "B"), "performance is not of class 'numeric'")
    expect_error(RPT(0.80,0.80,"C"), "beta is not of class 'numeric'")
    expect_error(RPT(1.1,.80), "stability is not in range 0 to 1.")
    expect_error(RPT(.80,1.1), "performance is not in range 0 to 1.")
    expect_error(RPT(.80,.80,-1), "beta is non-positive.")
})
