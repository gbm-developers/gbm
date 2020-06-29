library(OmicsMarkeR)
context("Aggregation Functions")

ranks <- replicate(5, sample(seq(50), 50))
row.names(ranks) <- paste0("V", seq(50))

agg_cla <- aggregation(ranks, "CLA")
agg_es <- aggregation(ranks, "ES")
agg_em <- aggregation(ranks, "EM")
agg_ee <- aggregation(ranks, "EE", f = 50)
agg_cla2 <- aggregation(ranks, "CLA", f = 50)
agg_es2 <- aggregation(ranks, "ES", f = 50)
agg_em2 <- aggregation(ranks, "EM", f = 50)

test_that("output is correct", {
    expect_is(agg_cla, "matrix")
    expect_is(agg_es, "matrix")
    expect_is(agg_em, "matrix")
    expect_is(agg_ee, "matrix")
    expect_is(agg_cla2, "matrix")
    expect_is(agg_es2, "matrix")
    expect_is(agg_em2, "matrix")
})
