pkgname <- "horserule"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('horserule')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("HorseRuleFit")
### * HorseRuleFit

flush(stderr()); flush(stdout())

### Name: HorseRuleFit
### Title: Horseshoe RuleFit
### Aliases: HorseRuleFit

### ** Examples

library(MASS)
library(horserule)
data(Boston)
# Split in train and test data
N = nrow(Boston)
train = sample(1:N, 400)
Xtrain = Boston[train,-14]
ytrain = Boston[train, 14]
Xtest = Boston[-train, -14]
ytest = Boston[-train, 14]

# Run the HorseRuleFit with 100 trees
# Increase Number of trees and the number of posterior samples for better modelfit
hrres = HorseRuleFit(X = Xtrain, y=ytrain,
                    thin=1, niter=100, burnin=10,
                    L=5, S=6, ensemble = "both", mix=0.3, ntree=100,
                    intercept=FALSE, linterms=1:13, ytransform = "log",
                    alpha=1, beta=2, linp = 1, restricted = 0)

# Calculate the error
pred = predict(hrres, Xtest, burnin=100, postmean=TRUE)
sqrt(mean((pred-ytest)^2))

# Look at the most important rules/linear effects.
importance_hs(hrres)

# Look at the input variable importance.
Variable_importance(hrres, var_names=colnames(Xtrain))



cleanEx()
nameEx("Variable_importance")
### * Variable_importance

flush(stderr()); flush(stdout())

### Name: Variable_importance
### Title: Variable Importance plot
### Aliases: Variable_importance

### ** Examples

#Fit HorseRuleFit
x = matrix(rnorm(5000), ncol=10)
y = apply(x,1,function(x)sum(x[1:5])+rnorm(1))
hrres = HorseRuleFit(X = x, y=y,
                     thin=1, niter=100, burnin=10,
                     L=5, S=6, ensemble = "both", mix=0.3, ntree=100,
                     intercept=FALSE, linterms=1:10,
                     alpha=1, beta=2, linp = 1, restricted = 0)
Variable_importance(hrres)



cleanEx()
nameEx("convergence_plot")
### * convergence_plot

flush(stderr()); flush(stdout())

### Name: convergence_plot
### Title: convergence_plot
### Aliases: convergence_plot

### ** Examples

library(MASS)
data(Boston)
#Split in train and test data
N = nrow(Boston)
train = sample(1:N, 400)
Xtrain = Boston[train,-14]
ytrain = Boston[train, 14]
Xtest = Boston[-train, -14]
ytest = Boston[-train, 14]

hrres = HorseRuleFit(X = Xtrain, y=ytrain,
                    thin=1, niter=100, burnin=10,
                    L=5, S=6, ensemble = "both", mix=0.3, ntree=100,
                    intercept=FALSE, linterms=1:13, ytransform = "log",
                    alpha=1, beta=2, linp = 1, restricted = 0)

#Check the model convergence out of sample
convergence_plot(hrres, Xtest, ytest, burnin = 10)



cleanEx()
nameEx("hs")
### * hs

flush(stderr()); flush(stdout())

### Name: hs
### Title: Horseshoe regression Gibbs-sampler
### Aliases: hs

### ** Examples

x = matrix(rnorm(1000), ncol=10)
y = apply(x,1,function(x)sum(x[1:5])+rnorm(1))
hsmod = hs(X=x, y=y, niter=100)



cleanEx()
nameEx("importance_hs")
### * importance_hs

flush(stderr()); flush(stdout())

### Name: importance_hs
### Title: Most important Rules/terms
### Aliases: importance_hs

### ** Examples

library(MASS)
library(horserule)
data(Boston)
# Split in train and test data
N = nrow(Boston)
train = sample(1:N, 400)
Xtrain = Boston[train,-14]
ytrain = Boston[train, 14]
Xtest = Boston[-train, -14]
ytest = Boston[-train, 14]
hrres = HorseRuleFit(X = Xtrain, y=ytrain,
                    thin=1, niter=100, burnin=10,
                    L=5, S=6, ensemble = "both", mix=0.3, ntree=100,
                    intercept=FALSE, linterms=1:13, ytransform = "log",
                    alpha=1, beta=2, linp = 1, restricted = 0)

#Create an importance table containing the 10 most important rules and linear terms
importance_hs(hrres, k=10)



cleanEx()
nameEx("predict.HorseRulemodel")
### * predict.HorseRulemodel

flush(stderr()); flush(stdout())

### Name: predict.HorseRulemodel
### Title: predict.hs
### Aliases: predict.HorseRulemodel

### ** Examples

x = matrix(rnorm(1000), ncol=10)
y = apply(x,1,function(x)sum(x[1:5])+rnorm(1))
hrresmod = HorseRuleFit(X=x, y=y, niter=100, burnin=10)
#predict training data to obtain the fitted values
predict(hrresmod, x, burnin=10)



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
