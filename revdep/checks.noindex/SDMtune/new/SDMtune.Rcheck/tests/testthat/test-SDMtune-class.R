context("Plot SDMtune Object")

res <- data.frame(fc = c("l", "q"), reg = c(1, 1), train_AUC = c(0.907, 0.910),
                  test_AUC = c(0.905, 0.904), diff_AUC = c(0.002, 0.006),
                  stringsAsFactors = FALSE)

o <- SDMtune(results = res, models = list(SDMtune:::bm_maxnet))

test_that("Non interactive plot is correct", {
  # Line plot for only one varying hyper
  # AUC
  p <- plot(o)
  expect_equal(p$labels$y, "AUC")
  expect_equal(p$labels$x, "model")
  expect_equal(p$labels$title, "")  # title is not present if not passed
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomPoint")
  expect_equal(class(p$layers[[2]]$geom)[1], "GeomLine")
  expect_equivalent(unique(p$data$type), c("Training", "Validation"))
  # TSS
  colnames(o@results) <- c("fc", "reg", "Train_TSS", "Test_TSS", "diff_TSS")
  p <- plot(o, title = "title")
  expect_equal(p$labels$y, "TSS")
  expect_equal(p$labels$title, "title")  # title present if passed
  # AICc
  o@results <- res[, c(1:5)]
  colnames(o@results) <- c("fc", "reg", "AICc", "delta_AICc")
  p <- plot(o)
  expect_equal(p$labels$y, "AICc")
  # Scatter plot for more than one varying hyper
  o@results$reg <- c(1, 2)
  p <- plot(o)
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomPoint")
  expect_length(p$layers, 1)
})

test_that("Interactive plot is created and can be saved", {
  o@results <- res
  p <- expect_s3_class(plot(o, interactive = TRUE), "SDMtuneChart")
  expect_true(file.exists(p))
  unlink(p, recursive = TRUE)
})
