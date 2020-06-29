# test.plotmo.R
# This does a basic sanity test of plotmo.
# For more comprehensive tests, see plotmo/inst/slowtests.
library(plotmo)
library(rpart)
data(kyphosis)
rpart.model <- rpart(Kyphosis~., data=kyphosis)
plotmo(rpart.model, type="vec", trace=1)
