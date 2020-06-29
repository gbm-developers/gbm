skip_on_cran()
skip_on_appveyor()

m <- SDMtune:::bm_maxent
m_cv <- SDMtune:::bm_maxent_cv

test_that("The function returns the expected output", {
  expect_true(tss(m) < 1)
  expect_true(is.numeric(tss(SDMtune:::bm_maxnet)))
  expect_true(tss(m_cv) < 1)
  expect_true(is.numeric(tss(m_cv)))
})

test_that("The function uses the testing dataset", {
  expect_true(tss(m) != tss(m, test = m_cv@models[[1]]@data))
  expect_true(tss(m_cv) > tss(m_cv, test = TRUE))
  expect_true(tss(m_cv) > tss(m_cv, test = SDMtune:::t))
})

test_that("The function raises errors", {
  expect_error(tss(m, SDMtune:::t@data),
               "\"test\" argument invalid, use an SWD object.")
  expect_error(tss(m_cv, SDMtune:::t@data),
               "\"test\" argument invalid, use an SWD object.")
})
