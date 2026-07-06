test_that("multinomial train.error matches fitted multinomial log loss", {
  set.seed(10)

  fit <- gbm(
    Species ~ .,
    data = iris,
    distribution = "multinomial",
    n.trees = 5,
    interaction.depth = 1,
    n.minobsinnode = 5,
    shrinkage = 0.1,
    bag.fraction = 1,
    train.fraction = 1,
    verbose = FALSE
  )

  actual_loss <- vapply(seq_len(fit$n.trees), function(n.trees) {
    prob <- predict(fit, iris, n.trees = n.trees, type = "response")[, , 1]
    -mean(log(prob[cbind(seq_len(nrow(iris)), as.integer(iris$Species))]))
  }, numeric(1))

  expect_equal(fit$train.error, actual_loss)
})

