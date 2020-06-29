context("impute")

test_that("Error for wrong input for option parameter",
          {
            missdata <- SimIm(parkinson, 0.1)
            expect_that(impute(missdata, lmFun = "lassoR", ini = 3)   , throws_error())
          })
