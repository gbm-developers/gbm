jk_auc <- data.frame(Variable = c("bio1", "bio12", "bio16"),
                     Train_AUC_without = c(0.933, 0.935, 0.934),
                     Train_AUC_withonly = c(0.851, 0.731, 0.712),
                     Test_AUC_without = c(0.930, 0.931, 0.929),
                     Test_AUC_withonly = c(0.928, 0.929, 0.926),
                     stringsAsFactors = FALSE)
jk_tss <- jk_auc
colnames(jk_tss) <- c("Variable", "Train_TSS_without", "Train_TSS_withonly",
                      "Test_TSS_without", "Test_TSS_withonly")
jk_aicc <- jk_auc[, c(1, 2, 4)]
colnames(jk_aicc) <- c("Variable", "AICc_without", "AICc_withonly")


test_that("The plot has the correct labels and elements", {
  # AUC
  p <- plotJk(jk_auc, ref = 0.925)
  expect_equal(p$labels$y, "Train AUC")
  expect_equal(p$layers[[2]]$data$yintercept, 0.925)
  expect_equal(unique(p$data$test), c("With only", "Without"))
  p <- plotJk(jk_auc, type = "test")
  expect_equal(p$labels$y, "Test AUC")
  # Only without
  p <- plotJk(jk_auc[, c(1, 2, 4)])
  expect_equal(unique(p$data$test), "Without")
  # TSS
  p <- plotJk(jk_tss)
  expect_equal(p$labels$y, "Train TSS")
  p <- plotJk(jk_tss, type = "test")
  expect_equal(p$labels$y, "Test TSS")
  # AICc
  p <- plotJk(jk_aicc)
  expect_equal(p$labels$y, "AICc")
  expect_error(plotJk(jk_aicc, type = "test"),
               "Test mode is not available with aicc!")
})

test_that("The function works also when the argument is a list", {
  jk <- list(results = jk_auc, models = list())
  p <- plotJk(jk)
  expect_equal(p$labels$y, "Train AUC")
})
