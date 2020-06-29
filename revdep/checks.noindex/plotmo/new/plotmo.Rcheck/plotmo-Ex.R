pkgname <- "plotmo"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('plotmo')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("plot_gbm")
### * plot_gbm

flush(stderr()); flush(stdout())

### Name: plot_gbm
### Title: Plot a gbm model
### Aliases: plot_gbm

### ** Examples

if (require(gbm)) {
    n <- 100                            # toy model for quick demo
    x1 <- 3 * runif(n)
    x2 <- 3 * runif(n)
    x3 <- sample(1:4, n, replace=TRUE)
    y <- x1 + x2 + x3 + rnorm(n, 0, .3)
    data <- data.frame(y=y, x1=x1, x2=x2, x3=x3)
    mod <- gbm(y~., data=data, distribution="gaussian",
               n.trees=300, shrinkage=.1, interaction.depth=3,
               train.fraction=.8, verbose=FALSE)

    plot_gbm(mod)

    # plotres(mod)                      # plot residuals

    # plotmo(mod)                       # plot regression surfaces
}



cleanEx()
nameEx("plot_glmnet")
### * plot_glmnet

flush(stderr()); flush(stdout())

### Name: plot_glmnet
### Title: Plot a glmnet model
### Aliases: plot_glmnet

### ** Examples

if (require(glmnet)) {
    x <- matrix(rnorm(100 * 10), 100, 10)   # n=100 p=10
    y <- x[,1] + x[,2] + 2 * rnorm(100)     # y depends only on x[,1] and x[,2]
    mod <- glmnet(x, y)

    plot_glmnet(mod)

    # plotres(mod)                          # plot the residuals
}



cleanEx()
nameEx("plotmo")
### * plotmo

flush(stderr()); flush(stdout())

### Name: plotmo
### Title: Plot a model's response over a range of predictor values (the
###   model surface)
### Aliases: plotmo
### Keywords: partial dependence regression

### ** Examples

if (require(rpart)) {
    data(kyphosis)
    rpart.model <- rpart(Kyphosis~., data=kyphosis)
    # pass type="prob" to plotmo's internal calls to predict.rpart, and select
    # the column named "present" from the matrix returned by predict.rpart
    plotmo(rpart.model, type="prob", nresponse="present")
}
if (require(earth)) {
    data(ozone1)
    earth.model <- earth(O3 ~ ., data=ozone1, degree=2)
    plotmo(earth.model)
    # plotmo(earth.model, pmethod="partdep") # partial dependence plots
}



cleanEx()
nameEx("plotres")
### * plotres

flush(stderr()); flush(stdout())

### Name: plotres
### Title: Plot the residuals of a regression model
### Aliases: plotres
### Keywords: partial dependence regression

### ** Examples

# we use lm in this example, but plotres is more useful for models
# that don't have a function like plot.lm for plotting residuals

lm.model <- lm(Volume~., data=trees)

plotres(lm.model)



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
