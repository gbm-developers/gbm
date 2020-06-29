# Load required packages
library(survival)

# Create some data
set.seed(2)
N <- 3000
X1 <- runif(N)
X2 <- runif(N)
X3 <- factor(sample(letters[1:4], N, replace = T))
mu <- c(-1, 0, 1, 2)[as.numeric(X3)]
f <- 0.5 * sin(3 * X1 + 5 * X2 ^ 2 + mu / 10)
tt.surv <- rexp(N, exp(f))
tt.cens <- rexp(N, 0.5)
delta <- as.numeric(tt.surv <= tt.cens)
tt <- apply(cbind(tt.surv, tt.cens), 1, min)

# Throw in some missing values
X1[sample(1:N, size = 100)] <- NA
X3[sample(1:N, size = 300)] <- NA

# Random weights if you want to experiment with them
w <- rep(1, N)

data <- data.frame(
  tt = tt,
  delta = delta,
  X1 = X1,
  X2 = X2,
  X3 = X3
)

# fit initial model
gbm1 <- gbm(
  Surv(tt, delta) ~ X1 + X2 + X3,  # formula
  data = data,  # dataset
  weights = w,
  var.monotone = c(0, 0, 0),  # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
  distribution = "coxph",
  n.trees = 3000, # number of trees
  shrinkage = 0.001, # shrinkage or learning rate, 0.001 to 0.1 usually work
  interaction.depth = 3, # 1: additive model, 2: two-way interactions, etc
  bag.fraction = 0.5, # subsampling fraction, 0.5 is probably best
  train.fraction = 0.5, # fraction of data for training, first train.fraction*N used for training
  cv.folds = 5, # do 5-fold cross-validation
  n.cores = 1, 
  n.minobsinnode = 10, # minimum total weight needed in each node
  keep.data = TRUE
)

# Extract optimal number of trees based on test set performance
best.iter <- gbm.perf(gbm1, method = "test", plot.it = FALSE) # returns test set estimate of best number of trees

# Make some new data
set.seed(2)
N <- 1000
X1 <- runif(N)
X2 <- runif(N)
X3 <- factor(sample(letters[1:4], N, replace = T))
mu <- c(-1, 0, 1, 2)[as.numeric(X3)]

f <- 0.5 * sin(3 * X1 + 5 * X2 ^ 2 + mu / 10)  # -0.5 <= f <= 0.5 via sin fn.
tt.surv <- rexp(N, exp(f))
tt.cens <- rexp(N, 0.5)

data2 <- data.frame(
  tt = apply(cbind(tt.surv, tt.cens), 1, min),
  delta = as.numeric(tt.surv <= tt.cens),
  f = f,
  X1 = X1,
  X2 = X2,
  X3 = X3
)

# predict on the new data using "best" number of trees
# f.predict will be on the canonical scale (logit,log,etc.)
f.predict <- predict(gbm1, newdata = data2, n.trees = best.iter)

#plot(data2$f,f.predict)
# Use observed sd
expect_true(sd(data2$f - f.predict) < 0.4,
            info = "CoxPH: checking if squared error within tolerance.")
