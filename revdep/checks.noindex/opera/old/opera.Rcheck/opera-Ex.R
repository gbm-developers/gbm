pkgname <- "opera"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('opera')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("electric_load")
### * electric_load

flush(stderr()); flush(stdout())

### Name: electric_load
### Title: Electricity forecasting data set
### Aliases: electric_load
### Keywords: datasets

### ** Examples

 data(electric_load)
 # a few graphs to display the data
 attach(electric_load)
 plot(Load, type = 'l')
 plot(Temp, Load, pch = 16, cex = 0.5)
 plot(NumWeek, Load, pch = 16, cex = 0.5)
 plot(Load, Load1, pch = 16, cex = 0.5)
 acf(Load, lag.max = 20)
 detach(electric_load)



cleanEx()
nameEx("mixture")
### * mixture

flush(stderr()); flush(stdout())

### Name: mixture
### Title: Compute an aggregation rule
### Aliases: mixture print.mixture summary.mixture
### Keywords: ~models ~ts

### ** Examples

#' 
library('opera')  # load the package
set.seed(1)

# Example: find the best one week ahead forecasting strategy (weekly data)
# packages
library(mgcv)

# import data
data(electric_load)
idx_data_test <- 620:nrow(electric_load)
data_train <- electric_load[-idx_data_test, ]
data_test <- electric_load[idx_data_test, ]

# Build the expert forecasts 
# ##########################

# 1) A generalized additive model
gam.fit <- gam(Load ~ s(IPI) + s(Temp) + s(Time, k=3) + 
                s(Load1) + as.factor(NumWeek), data = data_train)
gam.forecast <- predict(gam.fit, newdata = data_test)

# 2) An online autoregressive model on the residuals of a medium term model

# Medium term model to remove trend and seasonality (using generalized additive model)
detrend.fit <- gam(Load ~ s(Time,k=3) + s(NumWeek) + s(Temp) + s(IPI), data = data_train)
electric_load$Trend <- c(predict(detrend.fit), predict(detrend.fit,newdata = data_test))
electric_load$Load.detrend <- electric_load$Load - electric_load$Trend

# Residual analysis
ar.forecast <- numeric(length(idx_data_test))
for (i in seq(idx_data_test)) {
 ar.fit <- ar(electric_load$Load.detrend[1:(idx_data_test[i] - 1)])
 ar.forecast[i] <- as.numeric(predict(ar.fit)$pred) + electric_load$Trend[idx_data_test[i]]
}

# Aggregation of experts
###########################

X <- cbind(gam.forecast, ar.forecast)
colnames(X) <- c('gam', 'ar')
Y <- data_test$Load

matplot(cbind(Y, X), type = 'l', col = 1:6, ylab = 'Weekly load', xlab = 'Week')


# How good are the expert? Look at the oracles
oracle.convex <- oracle(Y = Y, experts = X, loss.type = 'square', model = 'convex')
plot(oracle.convex)
oracle.convex

# Is a single expert the best over time ? Are there breaks ?
oracle.shift <- oracle(Y = Y, experts = X, loss.type = 'percentage', model = 'shifting')
plot(oracle.shift)
oracle.shift

# Online aggregation of the experts with BOA
#############################################

# Initialize the aggregation rule
m0.BOA <- mixture(model = 'BOA', loss.type = 'square')

# Perform online prediction using BOA There are 3 equivalent possibilities 1)
# start with an empty model and update the model sequentially
m1.BOA <- m0.BOA
for (i in 1:length(Y)) {
 m1.BOA <- predict(m1.BOA, newexperts = X[i, ], newY = Y[i])
}

# 2) perform online prediction directly from the empty model
m2.BOA <- predict(m0.BOA, newexpert = X, newY = Y, online = TRUE)

# 3) perform the online aggregation directly
m3.BOA <- mixture(Y = Y, experts = X, model = 'BOA', loss.type = 'square')

# These predictions are equivalent:
identical(m1.BOA, m2.BOA)  # TRUE
identical(m1.BOA, m3.BOA)  # TRUE

# Display the results
summary(m3.BOA)
plot(m1.BOA)



cleanEx()
nameEx("opera-package")
### * opera-package

flush(stderr()); flush(stdout())

### Name: opera-package
### Title: Online Prediction by ExpeRt Aggregation
### Aliases: opera-package opera
### Keywords: package

### ** Examples

#' 
library('opera')  # load the package
set.seed(1)

# Example: find the best one week ahead forecasting strategy (weekly data)
# packages
library(mgcv)

# import data
data(electric_load)
idx_data_test <- 620:nrow(electric_load)
data_train <- electric_load[-idx_data_test, ]
data_test <- electric_load[idx_data_test, ]

# Build the expert forecasts 
# ##########################

# 1) A generalized additive model
gam.fit <- gam(Load ~ s(IPI) + s(Temp) + s(Time, k=3) + 
                s(Load1) + as.factor(NumWeek), data = data_train)
gam.forecast <- predict(gam.fit, newdata = data_test)

# 2) An online autoregressive model on the residuals of a medium term model

# Medium term model to remove trend and seasonality (using generalized additive model)
detrend.fit <- gam(Load ~ s(Time,k=3) + s(NumWeek) + s(Temp) + s(IPI), data = data_train)
electric_load$Trend <- c(predict(detrend.fit), predict(detrend.fit,newdata = data_test))
electric_load$Load.detrend <- electric_load$Load - electric_load$Trend

# Residual analysis
ar.forecast <- numeric(length(idx_data_test))
for (i in seq(idx_data_test)) {
 ar.fit <- ar(electric_load$Load.detrend[1:(idx_data_test[i] - 1)])
 ar.forecast[i] <- as.numeric(predict(ar.fit)$pred) + electric_load$Trend[idx_data_test[i]]
}

# Aggregation of experts
###########################

X <- cbind(gam.forecast, ar.forecast)
colnames(X) <- c('gam', 'ar')
Y <- data_test$Load

matplot(cbind(Y, X), type = 'l', col = 1:6, ylab = 'Weekly load', xlab = 'Week')


# How good are the expert? Look at the oracles
oracle.convex <- oracle(Y = Y, experts = X, loss.type = 'square', model = 'convex')
plot(oracle.convex)
oracle.convex

# Is a single expert the best over time ? Are there breaks ?
oracle.shift <- oracle(Y = Y, experts = X, loss.type = 'percentage', model = 'shifting')
plot(oracle.shift)
oracle.shift

# Online aggregation of the experts with BOA
#############################################

# Initialize the aggregation rule
m0.BOA <- mixture(model = 'BOA', loss.type = 'square')

# Perform online prediction using BOA There are 3 equivalent possibilities 1)
# start with an empty model and update the model sequentially
m1.BOA <- m0.BOA
for (i in 1:length(Y)) {
 m1.BOA <- predict(m1.BOA, newexperts = X[i, ], newY = Y[i])
}

# 2) perform online prediction directly from the empty model
m2.BOA <- predict(m0.BOA, newexpert = X, newY = Y, online = TRUE)

# 3) perform the online aggregation directly
m3.BOA <- mixture(Y = Y, experts = X, model = 'BOA', loss.type = 'square')

# These predictions are equivalent:
identical(m1.BOA, m2.BOA)  # TRUE
identical(m1.BOA, m3.BOA)  # TRUE

# Display the results
summary(m3.BOA)
plot(m1.BOA)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
