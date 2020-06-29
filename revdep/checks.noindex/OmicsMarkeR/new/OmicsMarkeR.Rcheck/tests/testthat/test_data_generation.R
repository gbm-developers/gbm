library(OmicsMarkeR)
context("Data Generation")

mat <- create.random.matrix(nvar = 50,
                            nsamp = 100,
                            st.dev = 1,
                            perturb = 0.2)
corMat <- create.corr.matrix(mat)
discrMat <- create.discr.matrix(corMat, D=10)

test_that("Inputs are correct", {
    mat[,2] <- as.character(rnorm(100))
    expect_error(create.corr.matrix(mat))
})

test_that("return warnings", {
    mat.cor <- matrix(rnorm(5000), nrow=100)
    expect_warning(create.discr.matrix(mat.cor, D=10))
    expect_warning(
        create.discr.matrix(corMat, D=10, num.groups=3))
})

test_that("outputs are correct", {
    expect_is(mat, "matrix")
    expect_is(corMat, "matrix.corr")
    expect_is(discrMat, "list")
})