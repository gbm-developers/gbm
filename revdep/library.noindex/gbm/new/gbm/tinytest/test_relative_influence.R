# Test that relative.influence really does pick out the true predictors
set.seed(1234)
X1 <- matrix(nrow = 1000, ncol = 50)
X1 <- apply(X1, 2, function(x) rnorm(1000))  # random noise
X2 <- matrix(nrow = 1000, ncol = 5)
X2 <- apply(X2, 2, function(x) c(rnorm(500), rnorm(500, 3)))  # real predictors
cls <- rep(c(0, 1), ea = 500) # Class
X <- data.frame(cbind(X1, X2, cls))
mod <- gbm(
  cls ~ .,
  data = X,
  n.trees = 1000,
  cv.folds = 5,
  n.cores = 1,
  shrinkage = .01,
  interaction.depth = 2
)
ri <- rev(sort(relative.influence(mod)))
wh <- names(ri)[1:5]
res <- sum(wh %in% paste("V", 51:55, sep = ""))
expect_identical(
  current = res, 
  target = 5L, 
  info = "Checking if relative influence identifies true predictors."
)
