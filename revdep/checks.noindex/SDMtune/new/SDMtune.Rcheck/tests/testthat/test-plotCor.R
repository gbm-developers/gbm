p <- plotCor(SDMtune:::t, cor_th = .8)

test_that("The plot has the correct labels", {
  expect_equal(p$plot_env$label, "Spearman's\ncoefficient")
  expect_equal(p$labels$x, "Var2")
  expect_equal(p$labels$y, "Var1")
})

test_that("The plot as the correct text when using cor_th argument", {
  expect_true(sum(abs(p$layers[[2]]$data$value) >= .8) == nrow(p$layers[[2]]$data))
})

test_that("The plot has the correct text when the argument cor_th is NULL", {
  p <- plotCor(SDMtune:::t)
  expect_identical(length(p$layers[[2]]$data$value), nrow(p$data))
})
