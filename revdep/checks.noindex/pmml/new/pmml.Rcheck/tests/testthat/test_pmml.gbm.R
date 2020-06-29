library(gbm)
data(audit)


test_that("pmml.gbm final Segment contains modelName attribute", {
  audit_dat <- audit[, -c(1, 4, 6, 9, 10, 11, 12)]
  fit <- gbm(Adjusted ~ ., data = audit_dat, n.trees = 3, interaction.depth = 4, distribution = "multinomial")
  p_fit <- pmml(fit)

  expect_equal(xmlGetAttr(p_fit[[3]][[3]][[7]][[2]], name = "modelName"), "GBM_Model")
  expect_equal(xmlGetAttr(p_fit[[3]][[3]][[7]][[2]], name = "normalizationMethod"), "softmax")
})
