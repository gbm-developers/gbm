# FOr reproducibility
set.seed(1)

# Create some data
N <- 1000
X1 <- runif(N)
X2 <- runif(N)
X3 <- factor(sample(letters[1:4], N, replace = T))
mu <- c(-1, 0, 1, 2)[as.numeric(X3)]
p <- 1 / (1 + exp(-(sin(3 * X1) - 4 * X2 + mu)))
Y <- rbinom(N, 1, p)
w <- rexp(N)
w <- N * w / sum(w)  # random weights if you want to experiment with them
data <- data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3)

# Fit initial model
gbm1 <- gbm(
  Y ~ X1 + X2 + X3,  # formula
  data = data,  # dataset
  weights = w,
  var.monotone = c(0, 0, 0),  # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
  distribution = "bernoulli",
  n.trees = 3000,  # number of trees
  shrinkage = 0.001,  # shrinkage or learning rate, 0.001 to 0.1 usually work
  interaction.depth = 3,  # 1: additive model, 2: two-way interactions, etc
  bag.fraction = 0.5,  # subsampling fraction, 0.5 is probably best
  train.fraction = 0.5,  # fraction of data for training, first train.fraction*N used for training
  cv.folds = 5,  # do 5-fold cross-validation
  n.cores = 1,
  n.minobsinnode = 10  # minimum total weight needed in each node
)       

# Extract optimal number of trees based on test set performance
best.iter.test <- gbm.perf(gbm1, method = "test", plot.it = FALSE) # returns test set estimate of best number of trees
best.iter <- best.iter.test

# Make some new data
set.seed(2)
N <- 1000
X1 <- runif(N)
X2 <- runif(N)
X3 <- factor(sample(letters[1:4], N, replace = T))
mu <- c(-1, 0, 1, 2)[as.numeric(X3)]
p <- 1 / (1 + exp(-(sin(3 * X1) - 4 * X2 + mu)))
Y <- rbinom(N, 1, p)
data2 <- data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3)

# Predict on the new data using "best" number of trees
# f.predict will be on the canonical scale (logit,log,etc.)
f.1.predict <- predict(gbm1, data2, n.trees = best.iter.test)

# Compute quantity prior to transformation
f.new = sin(3 * X1) - 4 * X2 + mu

# Base the validation tests on observed discrepancies
expect_true(sd(f.new - f.1.predict) < 1.0)
