# test.full.R: test earth

print(R.version.string)

library(earth)
library(mda)
source("check.models.equal.R")
set.seed(2016)
data(ozone1)
data(trees)
data(etitanic)
if(!interactive())
    postscript(paper="letter")

PRINT.TIME <- FALSE         # FALSE for no time results (for diff against reference)
PLOT <- TRUE                # TRUE to do plots too, FALSE for speed
options.old <- options()
options(warn=1) # print warnings as they occur
# options(digits=5) # removed because want to check against default

# test that we got an error as expected from a try() call
expect.err <- function(object, expected.msg="")
{
    if(class(object)[1] == "try-error") {
        msg <- attr(object, "condition")$message[1]
        if(length(grep(expected.msg, msg, fixed=TRUE)))
            cat("Got error as expected from ",
                deparse(substitute(object)), "\n", sep="")
        else
            stop(sprintf("Expected: %s\n  Got:      %s",
                         expected.msg, substr(msg, 1, 1000)))
    } else
        stop("did not get expected error for ", expected.msg)
}
dummy.plot <- function()
{
    plot(0, 0, col=0, bty="n", xaxt="n", yaxt="n", xlab="", ylab="", main="")
}
printh <- function(x, expect.warning=FALSE, max.print=0) # like print but with a header
{
    cat("===", deparse(substitute(x)), " ", sep="")
    if(expect.warning)
        cat(" expect warning -->")
    else if (NROW(x) > 1)
        cat("\n")
    if (max.print > 0)
        print(head(x, n=max.print))
    else
        print(x)
}

print(citation("earth"))

#--- test examples from man pages ------------------------------------------------------------

cat("--- earth.Rd -----------------------------\n")
example(earth)

set.seed(2015)

    train.subset <- sample(1:nrow(trees), .8 * nrow(trees))
    test.subset <- (1:nrow(trees))[-train.subset]

    earth.model <- earth(Volume ~ ., data = trees[train.subset,])

    # print R-Squared on the test data
    print(summary(earth.model, newdata=trees[test.subset,]))

    # manually calculate R-Squared on the test data (same as above call to summary)
    yhat <- predict(earth.model, newdata = trees[test.subset,])
    y <- trees$Volume[test.subset]
    printh(1 - sum((y - yhat)^2) / sum((y - mean(y))^2)) # print R-Squared

newrsq <- 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
stopifnot(abs(summary(earth.model, newdata=trees[test.subset,])$newrsq - newrsq) < 1e-10)

cars <- earth(mpg ~ ., data = mtcars, pmethod = "none", trace = 4)

get.used.pred.names <- function(obj) # obj is an earth object
{
  any1 <- function(x) any(x != 0)    # like any but no warning if x is double
  names(which(apply(obj$dirs[obj$selected.terms,,drop=FALSE],2,any1)))
}
printh(get.used.pred.names(cars))

a1 <- earth(survived ~ ., data=etitanic,   # c.f. Harrell "Reg. Mod. Strat." ch. 12
             degree=2, trace=1,
             glm=list(family=binomial))
printh(a1)

a1a <- earth(etitanic[,-2], etitanic[,2],  # equivalent but using earth.default
             degree=2, trace=1,
             glm=list(family=binomial))
printh(a1a)

a2 <- earth(pclass ~ ., data=etitanic, trace=1, glm=list(family=binomial))
printh(a2)

ldose <- rep(0:5, 2) - 2 # V&R 4th ed. p. 191
sex <- factor(rep(c("male", "female"), times=c(6,6)))
numdead <- c(1,4,9,13,18,20,0,2,6,10,12,16)
pair <- cbind(numdead, numalive=20 - numdead)

a3 <- earth(pair ~ sex + ldose, trace=1, pmethod="none",
            glm=list(family=binomial(link=probit), maxit=100))
printh(a3)

numdead2 <- c(2,8,11,12,20,23,0,4,6,16,12,14) # bogus data
doublepair <- cbind(numdead, numalive=20-numdead,
                    numdead2=numdead2, numalive2=30-numdead2)

a4 <- earth(doublepair ~ sex + ldose, trace=1, pmethod="none",
            glm=list(family="binomial"))
printh(a4)

counts <- c(18,17,15,20,10,20,25,13,12) # Dobson 1990 p. 93
outcome <- gl(3,1,9)
treatment <- gl(3,3)

a5 <- earth(counts ~ outcome + treatment, trace=1, pmethod="none",
            glm=list(family=poisson))
printh(a5)

a6 <- earth(numdead ~ sex + ldose, trace=1, pmethod="none",
            glm=list(family=gaussian(link=identity)))
printh(a6$coefficients == a6$glm.coefficients)  # all TRUE
printh(a6)

remove(ldose)
remove(sex)
remove(numdead)
remove(pair)
remove(numdead2)
remove(doublepair)
remove(counts)
remove(outcome)
remove(treatment)

printh(earth(cbind(Volume,lvol=log(Volume)) ~ ., data=trees))
attach(trees)
printh(earth(data.frame(Girth,Height), data.frame(Volume,lvol=log(Volume))))
detach(trees)

lm.fit <- lm(O3 ~ log(temp) + humidity*temp, data=ozone1)
printh(lm.fit)
plotmo(lm.fit, level=.95, trace=-1)
lm.fit2 <- lm(O3 ~ temp+ibh+doy, data=ozone1)
printh(lm.fit2)
plotmo(lm.fit2, all2=TRUE, clip=FALSE, trace=-1)

cat("--- print.default of earth object---------\n")
print.default(cars, digits=3)
cat("--- done print.default of earth object----\n")
if (PLOT)
    plot(cars)
library(mda)
(a <- fda(Species~., data=iris, method=earth, keepxy=TRUE))
if (PLOT)
    plot(a)
printh(summary(a$fit))
expect.err(try(printh(summary(a$fit, none.such1="xxx"))), "unrecognized argument") # summary.earth unrecognized argument "none.such1"
printh(summary(a$fit, style="bf", none.such2="xxx")) # Warning: format.earth ignored unrecognized argument "none.such2"
if (PLOT) {
    plot(a$fit, col.residuals=iris$Species, nresponse=1)
    plotmo(a$fit, nresponse=1, ylim=c(-1.5,1.5), clip=FALSE, trace=-1)
    plotmo(a$fit, nresponse=2, ylim=c(-1.5,1.5), clip=FALSE, trace=-1)
}
a <- update(a, nk=3) # not on man page
printh(a)
printh(summary(a$fit))
head(etitanic) # pclass and sex are unordered factors
earth(pclass ~ ., data=etitanic, trace=2)

cat("--- format.earth.Rd ----------------------\n")
as.func <- function( # convert expression string to func
               object, digits = 8, use.names = TRUE, ...)
  eval(parse(text=paste(
    "function(x)\n",
    "{\n",
    "if(is.vector(x))\n",
    "  x <- matrix(x, nrow = 1, ncol = length(x))\n",
    "with(as.data.frame(x),\n",
    format(object, digits = digits, use.names = use.names, style = "p", ...),
    ")\n",
    "}\n", sep = "")))
a <- earth(Volume ~ ., data = trees)
my.func <- as.func(a, use.names = FALSE)
printh(my.func(c(10,80)))     # yields 17.76888
printh(predict(a, c(10,80)))  # yields 17.76888, but is slower
example(format.earth)
a <- earth(Volume ~ ., data = trees)
cat(format(a)) # basic tests of format.earth
cat(format(a, digits=4))
# cat(format(a, use.names=FALSE))
cat(format(a, style="pmax"))
cat(format(a, style="max"))
cat(format(a, style="bf"))
cat(format(a, use.names=FALSE, style="p"))
cat(format(a, use.names=FALSE, style="m"))
a <- earth(Volume ~ Girth*Height, data = trees, pmethod="none")
cat(format(a))
cat(format(a, colon.char="*"))
a <- lm(Volume ~ ., data = trees)
cat(format(a)) # basic tests of format.lm
cat(format(a, digits=4))
cat(format(a, use.names=FALSE))
cat(format(a, style="p"))
cat(format(a, use.names=FALSE, style="p"))
a <- lm(Volume ~ Girth*Height, data = trees)
cat(format(a))
cat(format(a, colon.char="*"))
cat("--- mars.to.earth.Rd ----------------------\n")
example(mars.to.earth)
library(mda)
a <- mars(trees[,-3], trees[,3])
a <- mars.to.earth(a)
summary(a, digits = 2)
printh(summary(a, digits=2))
printh(summary(a, digits=2, style="bf"))
cat("--- plot.earth.models.Rd ----------------------\n")
if (PLOT)
    example(plot.earth.models)
cat("--- plot.earth.Rd ----------------------\n")
if (PLOT) {
    data(etitanic)
    a <- earth(survived ~ ., data=etitanic, glm=list(family=binomial))
    par(mfrow=c(2,2))
    plot(a$glm.list[[1]], caption="a$glm.list[[1]]")
    example(plot.earth)
}
cat("--- predict.earth.Rd ----------------------\n")
example(predict.earth)
cat("--- update.earth.Rd ----------------------\n")
example(update.earth)

cat("--- evimp.Rd -----------------------------\n")

par(mfrow=c(2,2))
cat('before calling evimp par("mar", "cex"):\n')
print(par("mar", "cex"))

example(evimp)

cat("--- plot.evimp.Rd ------------------------\n")

example(plot.evimp)

rownames(ev)[4] <- "a_long_variable_name"

plot(ev, main="plot.evimp with various options",
    cex.var = .8,
    type.nsubsets = "p",
    col.nsubsets = "red",
    lty.nsubsets = 2, # ignored because type.nsubsets="p"
    type.gcv = "l",
    col.gcv = "green",
    lty.gcv = 3,
    type.rss = "b",
    col.rss = "blue",
    lty.rss = 4,
    cex.legend = .8,
    x.legend = "topright",
    rh.col = "pink")

a <- earth(Volume ~ Girth, data = trees)
plot(evimp(a), main="plot.evimp with single var in model")

cat('after calling evimp par("mar", "cex"):\n')
print(par("mar", "cex"))
par(mfrow=c(1,1))

cat("--- test predict.earth -------------------\n")

a <- earth(Volume ~ ., data = trees)
cat("1a predict(a, c(10,80))\n")
printh(predict(a, c(10,80), trace=1))
cat("1b predict(a, c(10,10,80,80))\n")
printh(predict(a, c(10,10,80,80), trace=1))
cat("1c predict(a, c(10,11,80,81))\n")
printh(predict(a, c(10,11,80,81), trace=1))
cat("2 predict(a)\n")
printh(head(predict(a, trace=1)))
cat("3a predict(a, matrix(c(10,12), nrow=1, ncol=2))\n")
printh(predict(a, matrix(c(10,12), nrow=1, ncol=2), trace=1))
cat("3b predict(a, matrix(c(10,12), nrow=2, ncol=2, byrow=TRUE)\n")
printh(predict(a, matrix(c(10,12), nrow=2, ncol=2, byrow=TRUE), trace=1))
cat("3c predict(a, matrix(c(10,12,80,90), nrow=2, ncol=2))\n")
printh(predict(a, matrix(c(10,12,80,90), nrow=2, ncol=2), trace=1))
xpredict <- matrix(c(10,12,80,90), nrow=2, ncol=2)
colnames(xpredict) <- c("Girth", "Height")
cat("4 predict(a, xpredict with colnames)\n")
printh(predict(a, xpredict, trace=1))
cat("5 predict(a, as.data.frame(xpredict with colnames))\n")
printh(predict(a, as.data.frame(xpredict), trace=1))
# reverse dataframe columns (and their names), predict should deal with it correctly
xpredict <- as.data.frame(cbind(xpredict[,2], xpredict[,1]))
colnames(xpredict) <- c("Height", "Girth")
cat("6a predict(a, xpredict with reversed columns and colnames)\n")
printh(predict(a, xpredict, trace=1))
xpredict2 <- cbind(xpredict[,1], xpredict[,2]) # nameless matrix
cat("6b predict(a, xpredict2)\n")
printh(predict(a, xpredict2, trace=1))

# repeat but with x,y (not formula) call to earth

x1 <- cbind(trees$Girth, trees$Height)
colnames(x1) <- c("Girth", "Height")
a <- earth(x1, trees$Volume)
xpredict <- matrix(c(10,12,80,90), nrow=2, ncol=2)
cat("7a predict(a)\n")
printh(head(predict(a, trace=1)))
cat("7n predict(a, matrix(c(10,12,80,90), nrow=2, ncol=2)\n")
printh(predict(a, matrix(c(10,12,80,90), nrow=2, ncol=2), trace=1))
colnames(xpredict) <- c("Girth", "Height")
cat("8 predict(a, xpredict with colnames)\n")
printh(predict(a, xpredict, trace=1))
cat("9 predict(a, as.data.frame(xpredict with colnames))\n")
printh(predict(a, as.data.frame(xpredict), trace=1))
cat("--Expect warning from predict.earth: the variable names in 'data' do not match those in 'object'\n")
xpredict2 <- cbind(xpredict[,1], xpredict[,2])
colnames(xpredict2) <- c("none.such", "joe")
cat("10a predict(a, xpredict2)\n")
printh(predict(a, xpredict2, trace=1), expect.warning=TRUE)
cat("--Expect warning from predict.earth: the variable names in 'data' do not match those in 'object'\n")
xpredict2 <- cbind(xpredict[,1], xpredict[,2])
colnames(xpredict2) <- c("Height", "Girth") # reversed
cat("10b predict(a, xpredict2)\n")
printh(predict(a, xpredict2, trace=1), expect.warning=TRUE)

cat("--- test predict.earth with multiple response models-------------------\n")

a <- earth(cbind(Volume, Volume + 100) ~ ., data = trees)
cat("1a predict(a, c(10,80))\n")
printh(predict(a, c(10,80), trace=1))
cat("1b predict(a, c(10,10,80,80))\n")
printh(predict(a, c(10,10,80,80), trace=1))
cat("1c predict(a, c(10,11,80,81))\n")
printh(predict(a, c(10,11,80,81), trace=1))
cat("2 predict(a)\n")
printh(head(predict(a, trace=1)))
cat("3a predict(a, matrix(c(10,12), nrow=1, ncol=2))\n")
printh(predict(a, matrix(c(10,12), nrow=1, ncol=2), trace=1))
cat("3b predict(a, matrix(c(10,12), nrow=2, ncol=2, byrow=TRUE)\n")
printh(predict(a, matrix(c(10,12), nrow=2, ncol=2, byrow=TRUE), trace=1))
cat("3c predict(a, matrix(c(10,12,80,90), nrow=2, ncol=2))\n")
printh(predict(a, matrix(c(10,12,80,90), nrow=2, ncol=2), trace=1))
xpredict <- matrix(c(10,12,80,90), nrow=2, ncol=2)
colnames(xpredict) <- c("Girth", "Height")
cat("4 predict(a, xpredict with colnames)\n")
printh(predict(a, xpredict, trace=1))
cat("5 predict(a, as.data.frame(xpredict with colnames))\n")
printh(predict(a, as.data.frame(xpredict), trace=1))
# reverse dataframe columns (and their names), predict should deal with it correctly
xpredict <- as.data.frame(cbind(xpredict[,2], xpredict[,1]))
colnames(xpredict) <- c("Height", "Girth")
cat("6 predict(a, xpredict with reversed columns and colnames)\n")
printh(predict(a, xpredict, trace=1))

# repeat but with x,y (not formula) call to earth

x1 <- cbind(trees$Girth, trees$Height)
colnames(x1) <- c("Girth", "Height")
a <- earth(x1, cbind(trees$Volume, trees$Volume+100))
xpredict <- matrix(c(10,12,80,90), nrow=2, ncol=2)
cat("7a predict(a)\n")
printh(head(predict(a, trace=1)))
cat("7b predict(a, matrix(c(10,12,80,90), nrow=2, ncol=2)\n")
printh(predict(a, matrix(c(10,12,80,90), nrow=2, ncol=2), trace=1))
colnames(xpredict) <- c("Girth", "Height")
cat("8 predict(a, xpredict with colnames)\n")
printh(predict(a, xpredict, trace=1))
cat("9 predict(a, as.data.frame(xpredict with colnames))\n")
printh(predict(a, as.data.frame(xpredict), trace=1))
cat("--Expect warning from predict.earth: the variable names in 'data' do not match those in 'object'\n")
xpredict <- as.data.frame(cbind(xpredict[,2], xpredict[,1]))
colnames(xpredict) <- c("Height", "Girth")
cat("10 predict(a, xpredict)\n")
printh(predict(a, xpredict, trace=1), expect.warning=TRUE)

cat("--- earth.predict with NAs, with formula interface ---\n")

predict.with.message <- function(message, earth.model, newdata) {
    cat("predict.earth  ", message, ":\n", sep="")
    print(predict(earth.model, newdata=newdata, trace=1))
    cat("\n")
}

iris.earth <- earth(Petal.Width ~ Sepal.Length + Sepal.Width + Petal.Length, data=iris)
x <- iris[1,]
predict.with.message("formula interface and vector", iris.earth, newdata=x)
x$Sepal.Width <- as.numeric(NA)
predict.with.message("formula interface and vector with NA", iris.earth, newdata=x)
x <- iris[1,]
x$Petal.Width <- as.numeric(NA) # Petal.Width is unused in the earth model
predict.with.message("formula interface and vector with NA in unused variable", iris.earth, newdata=x)

x <- iris[1:3,]
predict.with.message("formula interface and matrix", iris.earth, newdata=x)
x[2,]$Sepal.Width <- as.numeric(NA)
predict.with.message("formula interface and matrix with NA", iris.earth, newdata=x)
x <- iris[1:3,]
x[2,]$Petal.Width <- as.numeric(NA) # Petal.Width is unused in the earth model
predict.with.message("formula interface and matrix with NA in unused variable", iris.earth, newdata=x)

cat("--- earth.predict with NAs, with matrix interface ---\n")

iris.earth <- earth(iris[,1:3], iris[,4])
x <- iris[1,]
predict.with.message("default interface and vector", iris.earth, newdata=x)
x$Sepal.Width <- as.numeric(NA)
predict.with.message("default interface and vector with NA", iris.earth, newdata=x)
x <- iris[1,]
x$Petal.Width <- as.numeric(NA) # Petal.Width is unused in the earth model
predict.with.message("default interface and vector with NA in unused variable", iris.earth, newdata=x)

x <- iris[1:3,]
predict.with.message("default interface and matrix", iris.earth, newdata=x)
x[2,]$Sepal.Width <- as.numeric(NA)
predict.with.message("default interface and matrix with NA", iris.earth, newdata=x)
x <- iris[1:3,]
x[2,]$Petal.Width <- as.numeric(NA) # Petal.Width is unused in the earth model
predict.with.message("default interface and matrix with NA in unused variable", iris.earth, newdata=x)

cat("--- test reorder.earth ----------------------\n")
a <- earth(O3 ~ ., data = ozone1, degree = 2)
earth:::reorder.earth(a, decomp = "none")
earth:::reorder.earth(a)   # defaults to decomp = "anova"
a$selected.terms[earth:::reorder.earth(a)]

cat("--- tests with ozone data ----------------------\n")

ozone.test <- function(itest, sModel, x, y, degree=2, nk=51,
                    plotit=PLOT, trace=0, smooth.col="red")
{
    fite <- earth(x, y, degree=degree, nk=nk, trace=trace)
    fitm <- mars(x, y, degree=degree, nk=nk)

    cat("itest",
        sprintf("%-3d", itest),
        sprintf("%-32s", sModel),
        "degree", sprintf("%-2d", degree), "nk", sprintf("%-3g", nk),
        "nTerms",  sprintf("%-2d", sum(fite$selected.terms != 0)),
        "of", sprintf("%-3d", nrow(fite$dirs)),
        "GRSq", sprintf("%4.2g", fite$grsq),
        "GRSq ratio", fite$grsq/mars.to.earth(fitm)$grsq,
        "\n")
    caption <- paste("itest ", itest, ": ", sModel, " degree=", degree, " nk=", nk, sep="")
    printh(summary(fite))
    printh(summary(fite, style="bf"))
    if(plotit) {
        fitme <- mars.to.earth(fitm)
        plotmo(fite, caption=paste("EARTH", caption), trace=-1)
        plotmo(fitme, caption=paste("MARS", caption), trace=-1)
        plot(fite, npoints=500, smooth.col=smooth.col, caption=paste("EARTH", caption), info=TRUE)
        plot(fitme, caption=paste("MARS", caption), info=TRUE)
        fitme <- update(fitme)  # generate model selection data
        plot.earth.models(list(fite, fitme), caption=paste(itest, ": Compare earth to mars ", sModel, sep=""))
    }
    fite
}
data(ozone1)
attach(ozone1)

x.global <- cbind(wind, humidity, temp, vis)
y <- doy
itest <- 1; ozone.test(itest, "doy ~ wind+humidity+temp+vis", x.global, y, degree=1, nk=21)

x.global <- cbind(wind, humidity, temp, vis)
y <- doy
itest <- itest+1; a91 <- ozone.test(itest, "doy ~ wind+humidity+temp+vis", x.global, y, degree=2, nk=21)

# this is a basic test of RegressAndFix (because this generates lin dep bx cols)

cat("--Expect warning from mda::mars: NAs introduced by coercion\n") # why do we get a warning?
x.global <- cbind(wind, exp(humidity))
y <- doy
# smooth.col is 0 else get loess errors
# trace==2 so we print the "Fixed rank deficient bx by removing 2 terms, 7 terms remain" message
# TODO why are we getting the rank deficient message?
itest <- itest+1; ozone.test(itest, "doy ~ wind+exp(humidity)", x.global, y, degree=1, nk=21, smooth.col=0, trace=2)

x.global <- cbind(vh,wind,humidity,temp,ibh,dpg,ibt,vis,doy)
y <- O3
itest <- itest+1; ozone.test(itest, "O3~.", x.global, y, degree=2, nk=21)

x.global <- cbind(vh,wind,humidity,temp,ibh,dpg,ibt,vis,doy)
y <- O3
itest <- itest+1; ozone.test(itest, "O3~., nk=51", x.global, y, degree=2, nk=51)

detach(ozone1)

cat("--- fast mars -----------------------------------\n")

printh(earth(O3 ~ ., data=ozone1, degree=2, nk = 31, fast.k = 0, fast.beta = 1))
printh(earth(O3 ~ ., data=ozone1, degree=2, nk = 31, fast.k = 0, fast.beta = 0))
printh(earth(O3 ~ ., data=ozone1, degree=2, nk = 31, fast.k = 5, fast.beta = 1))
printh(earth(O3 ~ ., data=ozone1, degree=2, nk = 31, fast.k = 5, fast.beta = 0))

cat("--- plot.earth and plot.earth.models ------------\n")

a <- earth(O3 ~ ., data=ozone1) # formula interface

if (PLOT)
    plot(a, caption="plot.earth test 1", col.rsq=3, smooth.col=4, qqline.col="pink",
         col.vline=1, col.npreds=0, nresiduals=100, cum.grid="grid",
         grid.col="lightblue", col.sel.grid="lightgreen")

set.seed(1)
if (PLOT) {
    plot(a, caption="plot.earth test 2", which=c(3,4,1), ylim=c(.2,.9),
         id.n=20, legend.pos=c(10,.6), pch=20, lty.vline=1, cex.legend=1,
         grid.col="lightblue")

    plot(a, caption="plot.earth test 3", which=2, main="test main")
}

a1 <- earth(ozone1[,c(2:4,10)], ozone1[,1])     # x,y interface

if (PLOT) {
    plot(a, caption="plot.earth test 4", id.n=1)
    set.seed(1)
    plot.earth.models(a, which=1, ylim=c(.4,.8), jitter=.01)

    plot.earth.models(a1)

    plot.earth.models(list(a, a1), col.cum=c(3,4),  col.grsq=c(1,2), col.rsq=c(3,4),
         col.npreds=1, col.vline=1, lty.vline=3,
         legend.pos=c(5,.4), legend.text=c("a", "b", "c"), cex.legend=1.3)
}

cat("--- plot.earth args -----------------------------\n")

test.plot.earth.args <- function()
{
    caption <- "test earth args"
    printh(caption)

    argtest <- earth(ozone1[,c(2:4,10)], ozone1[,1])

    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    par(mfrow=c(2,3))
    par(cex = 0.8)
    par(mar = c(3, 3, 3, 0.5)) # small margins and text to pack figs in
    par(mgp = c(1.6, 0.6, 0))  # flatten axis elements
    oma <- par("oma") # make space for caption
    oma[3] <- 2.4
    par(oma=oma)
    par(cex.main=1)

    plot(argtest, do.par=FALSE, which=1,
         main="default")

    mtext(caption, outer=TRUE, font=2)

    plot(argtest, do.par=FALSE, which=1,
         col.rsq=3, col.grsq=2,
         col.npreds="blue", grid.col="lightblue",
         main=sprintf("%s\n%s",
            "col.rsq=3, col.grsq=2, ",
            "col.npreds=\"lightblue\", col.sel.grid=\"gray\""))

    plot(argtest, do.par=FALSE, which=1,
         col.vline="pink", legend.pos="topleft",
         lty.grsq=2, lty.npreds=1, lty.vline=1,
         main=sprintf("%s\n%s",
            "col.vline=\"pink\", legend.pos=\"topleft\", ",
            "lty.grsq=2, lty.npreds=1, lty.vline=1"))

    plot(argtest, do.par=FALSE, which=1,
         legend.pos=NA, col.npreds=0,
         main="legend.pos=NA, col.npreds=0")

    plot(argtest, do.par=FALSE, which=1,
         legend.pos=0,
         main="legend.pos=0")
}
test.plot.earth.args()

cat("--- test minspan --------------------------------\n")

a.minspan2 <- earth(O3 ~ ., data=ozone1, minspan=2)
printh(summary(a.minspan2))

a.minspan0 <- earth(O3 ~ ., data=ozone1, minspan=0)
printh(summary(a.minspan0))

a.minspan.minus1 <- earth(O3 ~ ., data=ozone1, minspan=-1)
printh(summary(a.minspan.minus1))

a.minspan.minus3 <- earth(O3 ~ ., data=ozone1, minspan=-3)
printh(summary(a.minspan.minus3))

a.endspan80 <- earth(O3 ~ ., data=ozone1, endspan=80)
printh(summary(a.endspan80))

cat("--- test multiple responses ---------------------\n")

# this uses the global matrix data.global (data.global[,1:2] is the response)

test.two.responses <- function(itest, func1, func2,
    degree=2, nk=51, plotit=PLOT, test.rsq=TRUE, trace=0, minspan=0,
    test.mars.to.earth=FALSE, pmethod="backward")
{
    if(typeof(func1) == "character")
        funcnames <- paste("multiple responses", func1, func2)
    else
        funcnames <- paste("multiple responses", deparse(substitute(func1)), deparse(substitute(func2)))
    cat("itest", sprintf("%-3d", itest), funcnames,
        " degree", sprintf("%-2d", degree), "nk", sprintf("%-3g", nk), "\n\n")
    gc()
    fite <- earth(x=data.global[,c(-1,-2), drop=FALSE], y=data.global[,1:2],
                degree=degree, trace=trace, nk=nk, pmethod=pmethod, minspan=minspan)
    printh(fite)
    caption <- paste("itest ", itest, ": ", funcnames, " degree=", degree, " nk=", nk, sep="")
    if(plotit) {
        if(typeof(func1) == "character") {
            plotmo(fite, caption=caption, nresponse=1, trace=-1)
            plotmo(fite, nresponse=2, trace=-1)
        } else {
            plotmo(fite, func=func1, caption=caption, nresponse=1)
            plotmo(fite, func=func2, nresponse=2)
        }
        plot(fite, caption=caption, nresponse=1)
        plot(fite, nresponse=2)
    }
    cat("\n")
    if(test.mars.to.earth) {
        cat("Testing mars.to.earth with a multiple response model\n")
        fitm <- mars(data.global[,c(-1,-2), drop=FALSE], data.global[,1:2],
                     degree=degree, trace=(trace!=0), nk=nk)
        fitme <- mars.to.earth(fitm)
        printh(fitme)
        printh(summary(fitme))
        if(plotit) {
            plotmo(fitm, func=func1, caption=caption, nresponse=1, clip=FALSE)
            plotmo(fitm, func=func2, nresponse=2, clip=FALSE)
        }
# TODO following code causes error "nk" not found, looking in wrong environment?
#       cat("Expect warnings because of weights in the mars model\n")
#       fitm <- mars(data.global[,c(-1,-2), drop=FALSE], data.global[,1:2],
#                    degree=degree, trace=(trace!=0), nk=nk, wp=c(1,2))
#       fitme <- mars.to.earth(fitm)
#       printh(fitme)
#       printh(summary(fitme))
    }
    fite
}

N <- 100
set.seed(1)
x1 <- runif(N, -1, 1)
x2 <- runif(N, -1, 1)
x3 <- runif(N, -1, 1)
x4 <- runif(N, -1, 1)
x5 <- runif(N, -1, 1)

func1 <- function(x)
{
    sin(3 * x[,1]) + x[,2]
}
func7 <- function(x)    # just one predictor
{
    sin(5 * x[,1])
}
x.global <- cbind(                                     x1, x2)
data.global <- cbind(func1(x.global), func7(x.global), x1, x2)
colnames(data.global) = c("func1", "func7", "x1", "x2")
# expect pmethod="ex" cannot be used with multiple response models
expect.err(try(test.two.responses(itest, func1, func7, nk=51, degree=1, pmethod="ex")), "not allowed with multiple response models")
# expect pmethod="seq" cannot be used with multiple response models
expect.err(try(test.two.responses(itest, func1, func7, nk=51, degree=1, pmethod="seq")), "not allowed with multiple response models")
itest <- itest+1; a <- test.two.responses(itest, func1, func7, nk=51, degree=1)
printh(summary(a))
printh(summary(a, style="bf"))
if (PLOT) {
    plotmo(a, nresponse=1, trace=-1)     # test generation of caption based on response name
    plotmo(a, nresponse=2, trace=-1)
    plot(a, nresponse=1)
    plot(a, nresponse=2)
}
x.global <- cbind(                                     x1, x2)
data.global <- cbind(func1(x.global), func7(x.global), x1, x2)
colnames(data.global) = c("func1",
   "a.very.long.in.fact.extremely.long.response.name",
   "x1.a.very.long.in.fact.extremely.long.predictor.name",
   "x2")
itest <- itest+1; a <- test.two.responses(itest, func1, func7, nk=51, degree=3)
printh(summary(a))
print(evimp(a))
print.default(evimp(a))

eqn56 <- function(x) # Friedman MARS paper equation 56
{
    0.1 * exp(4*x[,1]) +
    4 / (1 + exp(-20*(x[,2]-0.5))) +
    3 * x[,3] +
    2 * x[,4] +
    x[,5]
}
neg.eqn56 <- function(x)
{
    -eqn56(x)
}

eqn56noise <- function(x)
{
    set.seed(ncol(x))
    eqn56(x) + rnorm(nrow(x),0,1)
}

neg.eqn56noise <- function(x)
{
    -eqn56noise(x)
}

robot.arm <- function(x) # Friedman Fast MARS paper
{
    l1     <- x[,1]
    l2     <- x[,2]
    theta1 <- x[,3]
    theta2 <- x[,4]
    phi    <- x[,5]

    x1 <- l1 * cos(theta1) - l2 * cos(theta1 + theta2) * cos(phi)
    y <-  l1 * sin(theta1) - l2 * sin(theta1 + theta2) * cos(phi)
    z <-  l2 *  sin(theta2) * sin(phi)

    sqrt(x1^2 + y^2 + z^2)
}
x.global <- cbind(                                           x1, x2, x3, x4, x5)
data.global <- cbind(eqn56=eqn56(x.global), neg.eqn56noise(x.global), x1, x2, x3, x4, x5)
colnames(data.global) = c("", "neg.eqn56noise", "x1", "x2", "x3", "x4", "x5")
itest <- itest+1; a <- test.two.responses(itest, eqn56, neg.eqn56noise, nk=51, degree=1)
print(evimp(a))
print.default(evimp(a))

x.global <- cbind(                                           x1, x2, x3, x4, x5)
data.global <- cbind(eqn56=eqn56(x.global), neg.eqn56noise(x.global), x1, x2, x3, x4, x5)
colnames(data.global) = NULL
itest <- itest+1; a70 <- test.two.responses(itest, eqn56, neg.eqn56noise, nk=51, degree=2)
printh(summary(a70))
printh(summary(a70, style="bf"))

N1 <- 100
set.seed(1)
x1. <- runif(N1, -1, 1)
x2. <- runif(N1, -1, 1)
x3. <- runif(N1, -1, 1)
x4. <- runif(N1, -1, 1)
x5. <- runif(N1, -1, 1)

x.global <- cbind(                                        (x1.+1)/2, (x2.+2)/2, pi*(x3.+1), pi*(x4.+1), pi*x5./2 )
data.global <- cbind(robot.arm(x.global), eqn56(x.global), (x1.+1)/2, (x2.+2)/2, pi*(x3.+1), pi*(x4.+1), pi*x5./2 )
colnames(x.global)    <- c(                "l1", "l2", "theta1", "theta2", "phi")
colnames(data.global) <- c("arm", "eqn56", "l1", "l2", "theta1", "theta2", "phi")
itest <- itest+1; test.two.responses(itest, robot.arm, eqn56, nk=51, degree=1)
itest <- itest+1; test.two.responses(itest, robot.arm, eqn56, nk=51, degree=2, test.mars.to.earth=TRUE)
itest <- itest+1; test.two.responses(itest, robot.arm, eqn56, nk=201, degree=1)
itest <- itest+1; test.two.responses(itest, robot.arm, eqn56, nk=201, degree=2)
itest <- itest+1; test.two.responses(itest, robot.arm, eqn56, nk=201, degree=10)

attach(ozone1)
x.global <- cbind(                wind, humidity, temp, ibh, dpg, ibt, vis)
data.global <- cbind(O3, doy, vh, wind, humidity, temp, ibh, dpg, ibt, vis)
itest <- itest+1; test.two.responses(itest, "O3", "doy", nk=51, degree=2)
detach(ozone1)

cat("--- formula based multiple response -------------\n")

a2 <- earth(cbind(O3,doy) ~ ., data=ozone1, degree=2)
if (PLOT) {
    plotmo(a2, nresponse=1, trace=-1)                  # TODO1 delete
    plotmo(a2, nresponse=1, trace=-1) # test generation of caption based on response name
    plotmo(a2, nresponse=2, trace=-1)
    plot(a2, nresponse=1) # TODO delete
    plot(a2, nresponse=1)
    plot(a2, nresponse=2)
}

cat("--- test plot.earth.models with multiple responses ---\n")

set.seed(1)
a <- earth(O3 ~ ., data=ozone1, degree=2)
a2 <- earth(cbind(O3,doy) ~ ., data=ozone1, degree=2)
b2 <- earth(cbind(O3,doy) ~ ., data=ozone1, degree=1)
if (PLOT) {
    plot.earth.models(list(a, a2), caption="plot.earth.models with multiple responses, list(a,a2)")
    plot.earth.models(list(a2, a), caption="plot.earth.models with multiple responses, list(a2,a)",
                      col.rsq=c(2,3), col.npreds=c(2,3))
    plot.earth.models(list(a2, b2), caption="plot.earth.models with multiple responses, list(a2,b2)",
                      col.rsq=c(2,3), col.npreds=c(4,5), jitter=.01, legend.pos="topleft")
}

cat("--- subset --------------------------------------\n")

set.seed(9)
train.subset <- sample(1:nrow(ozone1), .8 * nrow(ozone1))
test.subset <- (1:nrow(ozone1))[-train.subset]

# all the following models should be identical
a <- earth(ozone1[,-1], ozone1[,1], subset=train.subset, nprune=7, degree=2)
printh(a)
if (PLOT)
    plotmo(a, caption="test subset: earth(ozone1[,-1], ozone1[,1], subset=train.subset)", trace=-1)

a <- earth(ozone1[train.subset,-1], ozone1[train.subset,1], nprune=7, degree=2)
printh(a)
if (PLOT)
    plotmo(a, caption="test subset: earth(ozone1[train.subset,-1], ozone1[train.subset,1]", trace=-1)

a <- earth(O3 ~ ., data=ozone1, subset=train.subset, nprune=7, degree=2)
printh(a)
if (PLOT)
    plotmo(a, caption="test subset: earth(O3 ~ ., data=ozone1, subset=train.subset", trace=-1)

y <- ozone1[test.subset, 1]
yhat <- predict(a, newdata = ozone1[test.subset, -1])
printh(1 - sum((y - yhat)^2)/sum((y - mean(y))^2)) # print RSquared

cat("--- update -------------------------\n")

a <- earth(O3 ~ ., data=ozone1, degree=2)
printh(update(a, penalty = -1, ponly=TRUE))
printh(update(a, penalty = 10, ponly=TRUE))
a <- earth(O3 ~ ., data=ozone1, nk=31, pmethod="n", degree=2)
a.none <- printh(update(a, nprune=10, pmethod="n"))
printh(update(a.none, pmethod="b"))
printh(update(a.none, nprune=4, pmethod="e"))
a.updated <- update(a.none, nprune=10, pmethod="b")
printh(a.updated)
a.backwards <- update(a, nprune=10, pmethod="b")
printh(a.backwards)
printh(all.equal(a.updated$bx, a.backwards$bx))
a <- earth(O3 ~ ., data=ozone1, nk=31, nprune=10, pmethod="b", degree=2)
printh(a)
printh(all.equal(a$bx, a.backwards$bx))

cat("--- Auto.linpreds -----------------------------\n")

set.seed(2017)
x1 <- runif(10)
x2 <- runif(10)
y <- x1 + x2
data=data.frame(x1=x1, x2=x2, y=y)
OLD.PAR <- par(no.readonly=TRUE)
par(mfrow = c(6, 4), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))

expect.err(try(earth(y~., data=data, Auto.linpr=99)), "Auto.linpreds=99 but it should be FALSE, TRUE, 0, or 1")

a <- earth(y~., data=data, trace=2) # default Auto.linpreds=TRUE
print(summary(a, style="pmax"))
plotmo(a, extend=.3, ylim=c(.2, 1.7),
      do.par=FALSE, pt.col=2, jitter=0,
      main=c("default Auto.linpreds=T", ""))
dummy.plot()
dummy.plot()

a1 <- earth(y~., data=data, trace=2, Auto.linpreds=FALSE)
print(summary(a1, style="pmax"))
plotmo(a1, extend=.3, ylim=c(.2, 1.7),
      do.par=FALSE, pt.col=2, jitter=0,
      main=c("Auto.linpreds=F", ""))
dummy.plot()
dummy.plot()
stopifnot(isTRUE(all.equal(predict(a), predict(a1))))

a2 <- earth(y~., data=data, trace=2, linpreds=TRUE, Auto.linpreds=FALSE)
print(summary(a2, style="pmax"))
plotmo(a2, extend=.3, ylim=c(.2, 1.7),
      do.par=FALSE, pt.col=2, jitter=0,
      main=c("linpreds=T, Auto.linpreds=F", ""))
dummy.plot()
dummy.plot()
stopifnot(isTRUE(all.equal(predict(a), predict(a2))))

a3 <- earth(y~., data=data, linpreds="x1", Auto.linpreds=FALSE)
print(summary(a3, style="pmax"))
plotmo(a3, extend=.3, ylim=c(.2, 1.7),
      do.par=FALSE, pt.col=2, jitter=0,
      main=c("linpreds=x1, Auto.linpreds=F", ""))
dummy.plot()
dummy.plot()
stopifnot(isTRUE(all.equal(predict(a), predict(a3))))

a4 <- earth(y~., data=data, linpreds="x2", Auto.linpreds=FALSE)
print(summary(a4, style="pmax"))
plotmo(a4, extend=.3, ylim=c(.2, 1.7),
      do.par=FALSE, pt.col=2, jitter=0,
      main=c("linpreds=x2, Auto.linpreds=F", ""))
dummy.plot()
dummy.plot()
stopifnot(isTRUE(all.equal(predict(a), predict(a4))))

# x,y interface
a5 <- earth(data[,1:2], data[,3], Auto.linpreds=FALSE)
print(summary(a5, style="pmax"))
plotmo(a5, extend=.3, ylim=c(.2, 1.7),
      do.par=FALSE, pt.col=2, jitter=0,
      main=c("x,y interface", ""))
dummy.plot()
dummy.plot()
stopifnot(isTRUE(all.equal(as.vector(predict(a1)), as.vector(predict(a5)))))
par(OLD.PAR)

# more complicated example (with Auto.linpreds=TRUE, vh enters linearly in a degree2 term)
data(ozone1)
oz <- ozone1[1:50,]
mod.none1 <- earth(O3~., data=oz, degree=2, nk=15, pmethod="none") # default Auto.linpreds=TRUE
print(summary(mod.none1))
mod.none2 <- earth(O3~., data=oz, degree=2, nk=15, pmethod="none", Auto.linpreds=FALSE)
print(summary(mod.none2))
stopifnot(all.equal(predict(mod.none1), predict(mod.none2)))

# example figure in inst/doc
OLD.PAR <- par(no.readonly=TRUE)
par(mfrow=c(2,2), mar=c(4, 3.2, 3, 3), mgp=c(1.6, 0.6, 0), par(cex = 0.7))
set.seed(2017)
offset <- 98
data.autolin <- data.frame(x=offset+(1:10), y=offset+(1:10))
autolinFALSE <- earth(y~x, data=data.autolin, Auto.linpreds=FALSE)
print(summary(autolinFALSE, style="max"))
set.seed(2017) # for same jitter on this and previous graph
plotmo(autolinFALSE, extend=.3, do.par=FALSE, pt.col="red", lwd=2,
       main="Auto.linpreds = FALSE",
       xaxt="n", yaxt="n", jitter=1, cex.main=1,
       xlim=offset+c(-2,13), ylim=offset+c(-3,13))
legend(x="topleft", legend=c("data", "earth model"),
       lty=c(0, 1), lwd=c(0, 2), pch=c(20, NA), col=c("red", 1))
text(x=offset+3.8, y=offset-1.2, cex=.9, "The knot happens to be at the")
text(x=offset+4,   y=offset-2.4, cex=.9, "minimum value of the predictor")

autolinTRUE <- earth(y~x, data=data.autolin) # default Auto.linpreds=TRUE
print(summary(autolinTRUE, style="max"))
set.seed(2017) # for same jitter on this and next graph
plotmo(autolinTRUE, extend=.3, do.par=FALSE, pt.col="red", lwd=2,
       main="Auto.linpreds = TRUE   (default)",
       xaxt="n", yaxt="n", jitter=1, cex.main=1,
       xlim=offset+c(-2,13), ylim=offset+c(-3,13))
legend(x="topleft", legend=c("data", "earth model"),
       lty=c(0, 1), lwd=c(0, 2), pch=c(20, NA), col=c("red", 1))
text(x=offset+4, y=offset-2.4, cex=.9, "Same data as previous graph")
stopifnot(isTRUE(all.equal(predict(autolinTRUE), predict(autolinFALSE))))
par(OLD.PAR)

# test Auto.linpreds with data sent in by a user
ndata <- matrix(data=c(
-0.0781, -0.6109, -0.216, -1.5172, 0.8184, -1.1242,
-0.0781, -0.5885, -0.216, -1.3501, 0.8184, -0.8703,
-0.0781, -0.5885, -0.216, -1.3501, 0.8184, -0.9549,
-0.0781, -0.5885, -0.216, -1.3501, 1.4136, -0.8703,
-2.5759, -0.5885, 1.1665, -1.3501, 2.0089, -0.9549,
-2.5759, -0.5885, 1.1665, -1.3501, 2.0089, -0.8703,
-0.0781, -0.4937, -0.216, -0.9949, -0.372, -1.0396,
-0.0781, -0.4463, -0.216, -0.8278, -0.372, -0.447,
-0.0781, -0.4463, -0.216, -0.8278, -0.372, -0.701,
-0.0781, -0.4463, -0.216, -0.8278, -0.372, -0.6163,
-0.0781, -0.4463, -0.216, -0.8278, 0.8184, -0.447,
-0.0781, -0.4463, -0.216, -0.8278, 0.8184, -0.6163,
-0.0781, -0.4463, 1.1665, -0.8278, 0.8184, -0.447,
-0.0781, -0.4379, 1.1665, 0.2585, -0.372, -0.1085,
-0.0781, -0.2147, 1.1665, 0.0496, -0.372, -0.1085,
-0.0781, -0.2147, -0.216, 0.2585, -0.372, -0.0238,
-0.0781, -0.1589, -0.216, 0.2585, -0.372, -0.1931,
-0.0781, -0.1589, -0.216, 0.2585, -0.372, -0.1085,
-0.0781, -0.1589, 1.1665, 0.2585, -0.372, -0.1931,
-0.0781, -0.1589, -0.216, 0.2585, 0.8184, -0.1085,
-0.0781, -0.1589, -0.216, 0.2585, 0.8184, 0.0608,
-0.0781, -0.1589, -0.216, 1.0942, 0.8184, -0.0238,
-0.0781, 0.0643, 1.1665, 1.0942, -0.372, 0.2301,
-0.0781, 0.0643, -0.216, 1.0942, -1.5624, 0.3148,
-0.0781, 0.0643, -0.216, 1.0942, -0.9672, 0.1455,
-0.0781, 0.0643, 1.1665, 1.4284, 0.2232, 0.4841,
-0.0781, 0.1563, -0.216, 1.0942, -0.372, 0.5687,
2.4197, 0.3432, -0.216, 1.0942, -1.5624, 1.0766,
-0.0781, 0.3432, -0.216, 1.0942, -1.5624, 1.1613,
-0.0781, 0.3432, 1.1665, 1.0942, 0.2232, 0.738,
2.4197, 2.7145, -2.9811, 1.0942, -1.5624, 2.5156,
2.4197, 4.3884, -2.9811, 1.0942, -1.5624, 3.5314),
ncol=6)
colnames(ndata) <- c("x1", "x2", "x3", "x4", "x5", "y")
ndata <- as.data.frame(ndata)

cat("Auto.linpreds=TRUE pmethod=\"none\":\n")
auto.linpreds.true.pmethod.none <- earth(y~., data=ndata, degree=2, nk=21, trace=2, pmethod="none")
print(summary(auto.linpreds.true.pmethod.none, decomp="none"))
cat("\nAuto.linpreds=FALSE pmethod=\"none\":\n")
auto.linpreds.false.pmethod.none <- earth(y~., data=ndata, degree=2, nk=21, trace=2, Auto.linpreds=FALSE, pmethod="none")
print(summary(auto.linpreds.false.pmethod.none, decomp="none"))
stopifnot(isTRUE(all.equal(predict(auto.linpreds.true.pmethod.none), predict(auto.linpreds.false.pmethod.none))))

cat("\nAuto.linpreds=TRUE:\n")
auto.linpreds.true <- earth(y~., data=ndata, degree=2, nk=21, trace=2)
print(summary(auto.linpreds.true, decomp="none"))
cat("\nAuto.linpreds=FALSE:\n")
auto.linpreds.false <- earth(y~., data=ndata, degree=2, nk=21, trace=2, Auto.linpreds=FALSE)
print(summary(auto.linpreds.false, decomp="none"))
# following fails because of different pruning because of different term count
# stopifnot(isTRUE(all.equal(predict(auto.linpreds.true), predict(auto.linpreds.false))))

cat("--- Force.xtx.prune -----------------------------\n")

expect.err(try(earth(Volume ~ ., data = trees, Force.xtx.prune=TRUE, pmethod="ex")), "not allowed with") # pmethod="ex" cannot be used with Force.xtx.prune

m1 <- earth(Volume ~ ., data = trees)
m2 <- earth(Volume ~ ., data = trees, Force.xtx.prune=TRUE)
check.models.equal(m1, m2, "Force.xtx.prune test 1", check.subsets=FALSE)

m1 <- earth(O3 ~ wind+temp, data = ozone1, nk=51)
m2 <- earth(O3 ~ wind+temp, data = ozone1, nk=51, Force.xtx.prune=TRUE)
check.models.equal(m1, m2, "Force.xtx.prune test 2", check.subsets=FALSE)

# TODO The following exposes a bug in leaps(?).  It is described in
# check.one.term.per.step in the earth R code.  The test is commented out
# because this bug causes a discrepancy with Force.xtx.prune (although
# usually the bug does not cause any problems).
#
# m1 <- earth(O3 ~ ., data = ozone1, nk=51, degree=2, trace=5)
# m2 <- earth(O3 ~ ., data = ozone1, nk=51, degree=2, Force.xtx.prune=TRUE)
# check.models.equal(m1, m2, "Force.xtx.prune test 3", check.subsets=FALSE)

cat("--- extractAIC.earth ----------------------------\n")

a <-earth(O3 ~ ., data=ozone1, degree=2)
cat("Ignore 10 warnings: extractAIC.earth: using GCV instead of AIC\n")
printh(drop1(a), expect.warning=TRUE)
printh(drop1(a, warn=FALSE)) # repeat but with warnings suppressed

cat("--- fda and mda with earth -----------------------------------\n")

am <- fda(Species ~ ., data=iris, method=mars, degree=1, keepxy=TRUE)
printh(am)
a <- fda(Species ~ ., data=iris, method=earth, degree=1, keepxy=TRUE)
printh(a)
printh(confusion(a))
if (PLOT) {
    par(mar=c(3, 3, 2, .5))  # small margins and text to pack figs in
    par(mgp=c(1.6, 0.6, 0))  # flatten axis elements
    par(oma=c(0,0,4,0))      # make space for caption
    layout(rbind(c(1,1,0,0), c(2,3,4,5), c(6,7,8,9)), heights=c(2,1,1))
    plot(a)
    plotmo(a$fit, nresponse=1, ylim=c(-1.5,1.5), clip=FALSE, do.par=FALSE, trace=-1)
    plotmo(a$fit, nresponse=2, ylim=c(-1.5,1.5), clip=FALSE, do.par=FALSE, trace=-1)
    mtext("fda test", outer=TRUE, font=2, line=1.5, cex=1)
}

data(glass)
set.seed(123)
samp <- sample(c(1:214), size=100, replace=FALSE)
glass.train <- glass[samp,]
glass.test <- glass[-samp,]
am <- mda(Type ~ ., data=glass.train, method=mars,  keepxy=TRUE, degree=2)
a <-  mda(Type ~ ., data=glass.train, method=earth, keepxy=TRUE, degree=2, keep.fitted=TRUE)
printh(am)
printh(a)
cat("mda with mars  ", attr(confusion(am), "error"), "\n")
cat("mda with earth ", attr(confusion(a),  "error"), "\n")
if (PLOT) {
    plot(a$fit, caption="mda on glass data", nresponse=1)
    plotmo(a$fit, nresponse=9, clip=FALSE, ylim=NA, caption="mda on glass data", trace=-1)
}

cat("\n---- update and keepxy, formula interface --------------------------\n")

new.trees <- trees + c(1,2,3,4)
new.trees <- new.trees[, -c(20:23)]
a.formula <- earth(Volume ~ ., subset=rep(TRUE, nrow(trees)), data = trees)
cat("\nupdate(a, trace=1)\n")
a.formula.1update <- update(a.formula, trace=1)
a.formula.1  <- earth(Volume ~ ., subset=rep(TRUE, nrow(trees)), data = trees)
check.models.equal(a.formula.1update, a.formula.1, msg="a1update a1")

cat("\nupdate(a.formula, data=new.trees, trace=1)\n")
a.formula.2update <- update(a.formula, data=new.trees, trace=1)
a.formula.2  <- earth(Volume ~ ., subset=rep(TRUE, nrow(trees)), data = new.trees)
check.models.equal(a.formula.2update, a.formula.2, msg="a2update a2")

cat("\nupdate(a.formula, wp=2, trace=1)\n")
a.formula.3update <- update(a.formula, wp=2, trace=1)
a.formula.3  <- earth(Volume ~ ., subset=rep(TRUE, nrow(trees)), data = trees, wp=2)
check.models.equal(a.formula.3update, a.formula.3, msg="a3update a3")

cat("\nupdate(a.formula, subset=subset.new, trace=1)\n")
subset.new <- rep(TRUE, nrow(trees))
subset.new[1:4] = FALSE
a.formula.4update <- update(a.formula, subset=subset.new, trace=1)
a.formula.4  <- earth(Volume ~ ., data = trees, subset=subset.new)
check.models.equal(a.formula.4update, a.formula.4, msg="a4update a4")

# now use keepxy=TRUE

a.formula <- earth(Volume ~ ., wp=1, data = trees, keepxy=TRUE)

cat("\nupdate(a.formula, trace=1)\n")
a.formula.5update <- update(a.formula, trace=1)
a.formula.5  <- earth(Volume ~ ., wp=1, data = trees, keepxy=TRUE)
check.models.equal(a.formula.5update, a.formula.5, msg="a5update a5")

cat("\nupdate(a.formula, data=new.trees, trace=1)\n")
a.formula.6update <- update(a.formula, data=new.trees, trace=1)
a.formula.6  <- earth(Volume ~ ., wp=1, data = new.trees, keepxy=TRUE)
check.models.equal(a.formula.6update, a.formula.6, msg="a6update a6")

cat("\nupdate(a.formula, wp=2, trace=1)\n")
a.formula.7update <- update(a.formula, wp=2, trace=1)
a.formula.7  <- earth(Volume ~ ., wp=2, data = trees, keepxy=TRUE)
check.models.equal(a.formula.7update, a.formula.7, msg="a7update a7")

cat("\n----- update and keepxy, matrix interface--------------------------\n")

Volume <- trees$Volume
x <- cbind(trees$Height, trees$Volume)
colnames(x) <- c("Height", "Volume")

new.x <- cbind(new.trees$Height, new.trees$Volume)
colnames(new.x) <- c("Height", "Volume")

a <- earth(x, Volume, subset=rep(TRUE, nrow(trees)))
cat("\nupdate(a, trace=1)\n")
a1update <- update(a, trace=1)
a1  <- earth(x, Volume, subset=rep(TRUE, nrow(trees)))
check.models.equal(a1update, a1, msg="a1update a1")

cat("\nupdate(a, x=new.x, trace=1)\n")
a2update <- update(a, x=new.x, trace=1)
a2  <- earth(new.x, Volume, subset=rep(TRUE, nrow(trees)))
check.models.equal(a2update, a2, msg="a2update a2")

cat("\nupdate(a, wp=2, trace=0)\n")
a3update <- update(a, wp=2, trace=0)
a3  <- earth(x, Volume, subset=rep(TRUE, nrow(trees)), wp=2)
check.models.equal(a3update, a3, msg="a3update a3")

cat("\nupdate(a, subset=subset.new, trace=4)\n")
subset.new <- rep(TRUE, nrow(trees))
subset.new[1:4] = FALSE
a4update <- update(a, subset=subset.new, trace=4)
a4  <- earth(x, Volume, subset=subset.new)
check.models.equal(a4update, a4, msg="a4update a4")

# now use keepxy=TRUE

a <- earth(x, Volume, wp=1, keepxy=TRUE)

cat("\nupdate(a, trace=4)\n")
a5update <- update(a, trace=4)
a5  <- earth(x, Volume, wp=1, keepxy=TRUE)
check.models.equal(a5update, a5, msg="a5update a5")

cat("\nupdate(a, x=new.x, trace=4)\n")
a6update <- update(a, x=new.x, trace=4)
a6  <- earth(new.x, Volume, wp=1, keepxy=TRUE)
check.models.equal(a6update, a6, msg="\"allowed\" function a6update a6")

cat("\nupdate(a, wp=2)\n")
a7update <- update(a, wp=2)
a7  <- earth(x, Volume, wp=2, keepxy=TRUE)
check.models.equal(a7update, a7, msg="\"allowed\" function a7update a7")

cat("--- \"allowed\" argument -----------------\n")

example1  <- function(degree, pred, parents)
{
    pred != 2  # disallow predictor 2, which is "Height"
}
a1 <- earth(Volume ~ ., data = trees, allowed = example1)
printh(summary(a1))

example1a  <- function(degree, pred, parents, namesx)
{
    namesx[pred] != "Height"  # disallow "Height"
}
a1a <- earth(Volume ~ ., data = trees, allowed = example1a)
check.models.equal(a1, a1a, msg="\"allowed\" function a1 a1a")

iheight <- 0
example1b  <- function(degree, pred, parents, namesx, first)
{
    if (first) {
        iheight <<- which(namesx == "Height") # note use of <<- not <-
        if (length(iheight) != 1)
            stop("could not find Height in ", paste(namesx, collapse=" "))
    }
    pred != iheight
}
a1b <- earth(Volume ~ ., data = trees, allowed = example1b)
check.models.equal(a1, a1a, msg="\"allowed\" function a1 a1b")

example2 <- function(degree, pred, parents)
{
    # disallow humidity in terms of degree > 1
    # 3 is the "humidity" column in the input matrix
    if (degree > 1 && (pred == 3 || parents[3]))
        return(FALSE)
    TRUE
}
a2 <- earth(O3 ~ ., data = ozone1, degree = 2, allowed = example2)
printh(summary(a2))
example3 <- function(degree, pred, parents)
{
    # allow only humidity and temp in terms of degree > 1
    # 3 and 4 are the "humidity" and "temp" columns
    allowed.set = c(3,4)
    if (degree > 1 && (all(allowed.set != pred) || any(parents[-allowed.set])))
        return(FALSE)
    TRUE
}
a3 <- earth(O3 ~ ., data = ozone1, degree = 2, allowed = example3)
printh(summary(a3))

ihumidity <- 0
itemp <- 0
example3a  <- function(degree, pred, parents, namesx, first)
{
    if (first) {
        ihumidity <<- which(namesx == "humidity") # note use of <<- not <-
        if (length(ihumidity) != 1)
            stop("could not find humidity in ", paste(namesx, collapse=" "))
        itemp <<- which(namesx == "temp") # note use of <<- not <-
        if (length(itemp) != 1)
            stop("could not find temp in ", paste(namesx, collapse=" "))
    }
    # allow only humidity and temp in terms of degree > 1
    allowed.set = c(ihumidity, itemp)
    if (degree > 1 &&
           (all(allowed.set != pred) || any(parents[-allowed.set])))
        return(FALSE)
    TRUE
}
a3a <- earth(O3 ~ ., data = ozone1, degree = 2, allowed = example3)
check.models.equal(a3, a3a, msg="\"allowed\" function a3 a3a")

#--- no predictor in PREDICTORS is allowed to interact with any predictor in PARENTS
#--- but all other interactions are allowed

PREDICTORS <- c("age")
PARENTS <- c("survived", "parch")

example4 <- function(degree, pred, parents, namesx)
{
    if (degree > 1) {
        predictor <- namesx[pred]
        parents   <- namesx[parents != 0]
        if((any(predictor %in% PREDICTORS) && any(parents %in% PARENTS)) ||
           (any(predictor %in% PARENTS)    && any(parents %in% PREDICTORS))) {
            return(FALSE)
        }
    }
    TRUE
}
a4.allowed <- earth(sex~., data=etitanic, degree=2, allowed=example4)
printh(summary(a4.allowed))
plotmo(a4.allowed, caption="a4.allowed")

#--- predictors in PREDICTORS are allowed to interact with predictors in PARENTS
#--- but no other interactions are allowed

PREDICTORS <- c("age")
PARENTS    <- c("survived", "parch")

example5 <- function(degree, pred, parents, namesx)
{
    if (degree <= 1)
        return(TRUE)
    predictor <- namesx[pred]
    parents   <- namesx[parents != 0]
    if((any(predictor %in% PREDICTORS) && any(parents %in% PARENTS)) ||
       (any(predictor %in% PARENTS)    && any(parents %in% PREDICTORS))) {
        return(TRUE)
    }
    FALSE
}
a5.allowed <- earth(sex~., data=etitanic, degree=2, allowed=example5)
printh(summary(a5.allowed))
plotmo(a5.allowed, caption="a5.allowed")

# "allowed" function checks, these check error handling by forcing an error

expect.err(try(earth(Volume ~ ., data = trees, allowed = 99)), "argument is not a function")

example7  <- function(degree, pred) pred!=2
expect.err(try(earth(Volume ~ ., data = trees, allowed = example7)), "function does not have the correct number of arguments")

example8  <- function(degree, pred, parents99) pred!=2
expect.err(try(earth(Volume ~ ., data = trees, allowed = example8)), "function needs the following arguments")

example9  <- function(degree, pred, parents, namesx99) pred!=2
expect.err(try(earth(Volume ~ ., data = trees, allowed = example9)), "function needs the following arguments")

cat("--- beta cache -------------------------\n")

a1 <- earth(O3 ~ ., data = ozone1, degree = 3)
a2 <- earth(O3 ~ ., data = ozone1, degree = 3, Use.beta.cache=FALSE)
a1$call <- NULL
a2$call <- NULL
stopifnot(identical(a1, a2))

cat("--- test \"call\" printing in earth.default and summary.earth ---\n")
# we want to make sure that long x or y aren't printed but short ones are

x = c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
      0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
      0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
      0,1,2,3,4,5,6,7,8,9,0)

y = c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
      0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
      0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
      0,1,2,3,4,5,6,7,8,9,0)

a <- earth(x = x, y=y, trace=5)

a.longx  <- earth(x = c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
                        0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
                        0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
                        0,1,2,3,4,5,6,7,8,9,0),
                  y=y,
                  trace=4)

a.longy  <- earth(x = x,
                  y = c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
                        0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
                        0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
                        0,1,2,3,4,5,6,7,8,9,0),
                  trace=4)

a.longxy <- earth(x = c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
                        0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
                        0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
                        0,1,2,3,4,5,6,7,8,9,0),
                  y = c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
                        0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
                        0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,
                        0,1,2,3,4,5,6,7,8,9,0),
                  trace=4)
printh(summary(a))
printh(summary(a.longx))
printh(summary(a.longy))
printh(summary(a.longxy))
printh(summary(a.longxy, style="bf"))

cat("--- factors with x,y interface -------------------------\n")
# this also tests for integer variables in the input matrix
data(etitanic)
attach(etitanic)
a1 <- earth(pclass, sex, degree=2, trace=2)        # x=unordered y=unordered
printh(summary(a1))
printh(summary(a1, style="bf"))
if (PLOT)
    plot(a1)
a2 <- earth(sex, pclass, degree=2, trace=2)        # x=unordered y=unordered
printh(summary(a2))
if (PLOT)
    plot(a2, nresponse=1)
a3 <- earth(pclass, age, degree=2, trace=2)        # x=unordered y=numeric
printh(summary(a3))
if (PLOT)
    plot(a3, nresponse=1)
a4 <- earth(age, pclass, degree=2, trace=2)        # x=numeric y=unordered
printh(summary(a4))
if (PLOT)
    plot(a4, nresponse=1)
a5 <- earth(etitanic[,c(2:4)], pclass, degree=2, trace=2)  # x=mixed  y=unordered
printh(summary(a5))
if (PLOT)
    plot(a5, nresponse=1)
a6 <- earth(etitanic[,c(1,3,4,5,6)], survived, degree=2, trace=2)  # x=mixed y=unordered
printh(summary(a6))
if (PLOT)
    plot(a6)
a7 <- earth(etitanic[,c(2,3,5,6)], etitanic[,c(1,4)], degree=2, trace=2)  # x=mixed y=mixed
printh(summary(a7))
if (PLOT)
    plot(a7, nresponse=1)

cat("--- factors with formula interface -------------------------\n")
# these correspond to the models above (except a7 which is a multiple response model)
a1f <- earth(sex ~ pclass, degree=2, trace=2)        # x=unordered y=unordered
printh(summary(a1f))
printh(summary(a1f, style="bf"))
if (PLOT)
    plot(a1f)
a2f <- earth(pclass ~ sex, degree=2, trace=2)        # x=unordered y=unordered
printh(summary(a2f))
if (PLOT)
    plot(a2f, nresponse=1)
a3f <- earth(age ~ pclass, degree=2, trace=2)        # x=unordered y=numeric
printh(summary(a3f))
if (PLOT)
    plot(a3f)
a4f <- earth(pclass ~ age, degree=2, trace=2)        # x=numeric y=unordered
printh(summary(a4f))
if (PLOT)
    plot(a4f, nresponse=1)
a5f <- earth(pclass ~ survived + sex + age, data=etitanic, degree=2, trace=2)  # x=mixed y=unordered
printh(summary(a5f))
if (PLOT)
    plot(a5f, nresponse=1)
a6f <- earth(survived ~ ., data=etitanic, degree=2, trace=2)  # x=mixed y=unordered
printh(summary(a6f))
if (PLOT)
    plot(a6f)
detach(etitanic)

# basic test with ordered factors
ff <- factor(substring("statistics", 1:10, 1:10), levels=letters, ordered=TRUE)
ff <- c(ff, ff, ff)
vowels = (ff == 1 | ff == 9) * 3
printh(head(ff))
printh(head(vowels))
a8 <- earth(ff, vowels, degree=1, trace=2)        # x=ordered y=numeric
printh(summary(a8))
if (PLOT)
    plot(a8, nresponse=1)
a9 <- earth(vowels, ff, degree=1, trace=2)        # x=numeric y=ordered
if (PLOT)
    plot(a9, nresponse=1)
printh(summary(a9))

cat("--- wp argument---------------------------------\n")
set.seed(79)
NWP <- 100
x1 <- runif(NWP)
x2 <- runif(NWP)
y1 <- (x1 > .5) + .3 * runif(1)
y2 <- sin(3 * x2) + .3 * runif(1)
myw <- 10
m <- mars(cbind(x1,x2), cbind(y1, y2))
me1 <- mars.to.earth(m)
printh(me1)
e1 <- earth(cbind(x1,x2), cbind(y1, y2))
printh(e1)
e2 <- earth(cbind(x1,x2), cbind(y1, y2),  wp=c(1,1))
printh(e2)
e1$call <- NULL
e2$call <- NULL
e1$wp <- NULL
e2$wp <- NULL
stopifnot(identical(e1, e2))
e3 <- earth(cbind(x1,x2), cbind(y1, y2),  wp=c(.001,1))
printh(e3)
wp <- c(1, 2)
e3 <- earth(cbind(x1,x2), cbind(y1, y2),  wp=wp)
printh(e3)
m3 <- mars(cbind(x1,x2), cbind(y1, y2),  wp=wp)
cat("response weights: wp", wp, "earth gcv", e3$gcv,
    "mars gcv", m3$gcv, "mars gcv*length(wp)",
    m3$gcv * length(wp), "\n")

e4 <- earth(cbind(O3, O3) ~ ., data=ozone1, wp=c(1, .01))
printh(e4) # both sub models should be the same
printh(summary(e4))

# wp with formula interface
e5 <- earth(cbind(O3, wind) ~ ., data=ozone1, wp=c(1, 1))
printh(e5)
printh(summary(e5))
e5 <- earth(cbind(O3, wind) ~ ., data=ozone1, wp=c(.3, 1))
printh(e5)
printh(summary(e5))
# wp with factors
e6 <- earth(pclass ~ ., data=etitanic, degree=2, wp=c(.001,.001,1))
printh(e6)
printh(summary(e6))
e7 <- earth(pclass ~ ., data=etitanic, degree=2, wp=c(1,.001,.001))
printh(e7)
printh(summary(e7))
if (PLOT)
    plot(e7, pt.col=as.numeric(etitanic$pclass)+1, nresponse=1)

cat("--- earth.regress ---------------------------------\n")

msg = "earth.regress with trees data, single response, no weights"
cat("Test:", msg, "\n")

data(trees)
y <- trees$Volume
x <- cbind(trees$Girth, trees$Height)
colnames(x) <- c("girth", "height")

a.lm <- lm(y ~ x)
a.lm.rss <- sum((a.lm$fitted.values - y)^2)
if (is.null(dim(a.lm$coefficients)))
    dim(a.lm$coefficients) <- c(length(a.lm$coefficients), 1)
a <- earth:::earth.regress(x, y)
rownames(a.lm$coefficients) <- rownames(a$coefficients)
check.fuzzy.equal(a.lm$coefficients, a$coefficients, msg=paste("coefficients [", msg, "]", sep=""))
check.fuzzy.equal(a.lm.rss, a$rss, msg=paste("rss [", msg, "]"))
check.fuzzy.equal(a.lm$residuals, a$residuals, msg=paste("residuals [", msg, "]"))

msg = "earth.regress with ozone1 data, multiple responses, no weights"
cat("Test:", msg, "\n")

data(ozone1)
y <- cbind(ozone1$O3, ozone1$O3 ^ 2)
colnames(y) <- c("O3", "O32")
x <- cbind(ozone1$wind, ozone1$humidity, ozone1$temp)
colnames(x) <- c("wind", "humidity", "temp")

a.lm <- lm(y ~ x)
a.lm.rss <- sum((a.lm$fitted.values - y)^2)
a <- earth:::earth.regress(x, y)
rownames(a.lm$coefficients) <- rownames(a$coefficients)
check.fuzzy.equal(a.lm$coefficients, a$coefficients, msg=paste("coefficients [", msg, "]"))
check.fuzzy.equal(a.lm.rss, a$rss, msg=paste("rss [", msg, "]", sep=""))
check.fuzzy.equal(a.lm$residuals, a$residuals, msg=paste("residuals [", msg, "]", sep=""))

# msg = "earth.regress with ozone1 data, multiple responses with case weights"
# cat("Test:", msg, "\n")
#
# # options(digits=10)
# weights. <- rep(.5, nrow(x))
# weights.[1] <- 1
# weights.[2] <- 2
# weights.[3] <- 3
# weights.[4] <- 4
# weights.[5] <- 5
# a.lm <- lm(y ~ x, weights=weights.)
# # a.lm.rss <- sum((a.lm$fitted.values - y)^2) # line below is equivalent
# a.lm.rss <- sum(a.lm$residuals^2)
# a <- earth:::earth.regress(x, y, weights=weights.)
# rownames(a.lm$coefficients) <- rownames(a$coefficients)
# check.fuzzy.equal(a.lm$coefficients, a$coefficients, msg=paste("coefficients [", msg, "]", sep=""))
# check.fuzzy.equal(a.lm.rss, a$rss, msg=paste("rss [", msg, "]", sep=""))
# check.fuzzy.equal(a.lm$residuals, a$residuals, msg=paste("residuals [", msg, "]", sep=""))

# msg = "earth.regress case weights with zero weights 1"
# cat("Test:", msg, "\n")
#
# weights. <- rep(1, nrow(x))
# weights.[2] <- 0
# weights.[4] <- 0
# a.lm <- lm(y ~ x, weights=weights.)
# # a.lm.rss <- sum((a.lm$fitted.values - y)^2) # line below is equivalent
# a.lm.rss <- sum(a.lm$residuals^2)
# a <- earth:::earth.regress(x, y, weights=weights.)
# rownames(a.lm$coefficients) <- rownames(a$coefficients)
# # options(digits=10)
# check.fuzzy.equal(a.lm$coefficients, a$coefficients, msg=paste("coefficients [", msg, "]", sep=""))
# check.fuzzy.equal(a.lm.rss, a$rss, msg=paste("rss [", msg, "]", sep=""))
# check.fuzzy.equal(a.lm$residuals, a$residuals, fuzz=1e-6, msg=paste("residuals [", msg, "]", sep=""))
#
# msg = "earth.regress case weights with zero weights 2"
# cat("Test:", msg, "\n")
# weights. <- rep(1, nrow(x))
# weights.[5] <- 0
# weights.[6] <- 0
# weights.[7] <- 0
# weights.[21] <- 0
# weights.[22] <- 0
# weights.[23] <- 0
# weights.[24] <- 0
# weights.[25] <- 0
# weights.[26] <- 0
# weights.[27] <- 0
# a.lm <- lm(y ~ x, weights=weights.)
# # a.lm.rss <- sum((a.lm$fitted.values - y)^2) # line below is equivalent
# a.lm.rss <- sum(a.lm$residuals^2)
# a <- earth:::earth.regress(x, y, weights=weights.)
# rownames(a.lm$coefficients) <- rownames(a$coefficients)
# check.fuzzy.equal(a.lm$coefficients, a$coefficients, msg=paste("coefficients [", msg, "]", sep=""))
# check.fuzzy.equal(a.lm.rss, a$rss, msg=paste("rss [", msg, "]", sep=""))
# check.fuzzy.equal(a.lm$residuals, a$residuals, fuzz=1e-6, msg=paste("residuals [", msg, "]", sep=""))
#
# msg = "earth.regress case weights with zero weights and missing columns 1"
# cat("Test:", msg, "\n")
# x <- cbind(ozone1$wind, ozone1$humidity, ozone1$temp, ozone1$wind^2, ozone1$humidity^2, ozone1$temp^2)
# weights. <- rep(1, nrow(x))
# weights.[5] <- 0
# weights.[6] <- 0
# weights.[7] <- 0
# weights.[21] <- 0
# weights.[22] <- 0
# weights.[23] <- 0
# weights.[24] <- 0
# weights.[25] <- 0
# weights.[26] <- 0
# weights.[27] <- 0
# colnames(x) <- c("wind", "humidity", "temp", "wind2", "humidity2", "temp2")
# used.cols = as.logical(c(1,0,1,0,1,1))
# x.missing <- x[,used.cols]
# a.lm <- lm(y ~ x.missing, weights=weights.)
# a.lm.rss <- sum((a.lm$fitted.values - y)^2) # line below is equivalent
# a.lm.rss <- sum(a.lm$residuals^2)
# a <- earth:::earth.regress(x, y, weights=weights., used.cols=used.cols)
# rownames(a.lm$coefficients) <- rownames(a$coefficients)
# check.fuzzy.equal(a.lm$coefficients, a$coefficients, msg=paste("coefficients [", msg, "]", sep=""))
# check.fuzzy.equal(a.lm.rss, a$rss, msg=paste("rss [", msg, "]", sep=""))
# check.fuzzy.equal(a.lm$residuals, a$residuals, fuzz=1e-6, msg=paste("residuals [", msg, "]", sep=""))
#
# msg = "earth.regress case weights with zero weights and missing columns 2"
# cat("Test:", msg, "\n")
# x <- cbind(ozone1$wind, ozone1$humidity, ozone1$temp, ozone1$wind^2, ozone1$humidity^2, ozone1$temp^2)
# weights. <- rep(1, nrow(x))
# weights.[5] <- .1
# weights.[6] <- .2
# weights.[7] <- 1.9
# weights.[21] <- .59
# colnames(x) <- c("wind", "humidity", "temp", "wind2", "humidity2", "temp2")
# used.cols = as.logical(c(0,1,0,0,1,0))
# x.missing <- x[,used.cols]
# a.lm <- lm(y ~ x.missing, weights=weights.)
# a.lm.rss <- sum((a.lm$fitted.values - y)^2) # line below is equivalent
# a.lm.rss <- sum(a.lm$residuals^2)
# a <- earth:::earth.regress(x, y, weights=weights., used.cols=used.cols)
# rownames(a.lm$coefficients) <- rownames(a$coefficients)
# check.fuzzy.equal(a.lm$coefficients, a$coefficients, msg=paste("coefficients [", msg, "]", sep=""))
# check.fuzzy.equal(a.lm.rss, a$rss, msg=paste("rss [", msg, "]", sep=""))
# check.fuzzy.equal(a.lm$residuals, a$residuals, fuzz=1e-6, msg=paste("residuals [", msg, "]", sep=""))

cat("---standard method functions ------------------------\n")

short.etitanic <- etitanic[seq(from=1, to=1000, by=20),]
a1 <- earth(pclass ~ ., data=short.etitanic, glm=list(family=binomial), trace=0)
printh(variable.names(a1))
printh(case.names(a1))
printh(case.names(a1, use.names=FALSE))

named.short.etitanic <- short.etitanic
rownames(named.short.etitanic) <- paste("xx", 1:nrow(named.short.etitanic))
a2 <- earth(pclass ~ ., data=named.short.etitanic, glm=list(family=binomial), trace=0)
printh(variable.names(a2))
printh(case.names(a2))
printh(case.names(a2, use.names=FALSE))

printh(deviance(a1), expect.warning=TRUE)
printh(deviance(a1, warn=FALSE))
printh(effects(a1), expect.warning=TRUE)
printh(effects(a1, warn=FALSE))
printh(family(a1))
printh(anova(a1), expect.warning=TRUE)
printh(anova(a1, warn=FALSE))
printh(family(a1))

# TODO removed because causes different results on different machines
# cat("--- thresh=0 -----------------------------------------\n")
#
# a.no.thresh <- earth(O3 ~ ., data = ozone1, thresh=0, nk=1000, degree=2, trace=4)
# printh(a.no.thresh)
# printh(summary(a.no.thresh))
# plotmo(a.no.thresh, degree1=1, degree2=c(4,9,16), clip=0, , caption="test with thresh=0", trace=-1)

# test the way plotmo gets the data with earth with a formula interface
# use strange data name se to make sure eval gets correct environment (don't pick up se in plotmo)
se <- ozone1
a <- earth(O3 ~ ., data=se, degree=2, keepxy=0)
printh(summary(a))
plotmo(a, trace=2, caption="getdata earth test1")
a <- earth(O3 ~ ., data=se, degree=2, keepxy=1)
printh(summary(a))
plotmo(a, trace=1, caption="getdata earth test2")
a <- earth(O3 ~ ., data=se, degree=2, keepxy=1)
se <- NULL
printh(summary(a))
plotmo(a, trace=2, caption="getdata earth test3")
se <- ozone1
a <- earth(O3 ~ ., data=se, degree=2, keepxy=0)
se <- NULL
printh(summary(a))
expect.err(try(plotmo(a, trace=0, caption="getdata earth test4")), "cannot get the original model predictors")

# test the way plotmo gets the data with earth with the default interface
se <- ozone1
a <- earth(se[,2:10], se[,1], degree=2, keepxy=0)
printh(summary(a))
plotmo(a, trace=0, caption="getdata earth test5")
a <- earth(se[,2:10], se[,1], degree=2, keepxy=1)
printh(summary(a))
plotmo(a, trace=0, caption="getdata earth test6")
a <- earth(se[,2:10], se[,1], degree=2, keepxy=1)
se <- NULL
printh(summary(a))
plotmo(a, trace=0, caption="getdata earth test7")
se <- ozone1
a <- earth(se[,2:10], se[,1], degree=2, keepxy=0)
se <- NULL
expect.err(try(plotmo(a, trace=0, caption="getdata earth test8")), "cannot get the original model predictors")
se <- ozone1
a <- earth(se[,2:10], se[,1], degree=2, keepxy=0)
# TODO error message could be improved here
se$vh <- NULL # vh is unused (but plotmo still needs it --- why?)
expect.err(try(plotmo(a, trace=0, caption="getdata earth test9")), "cannot get the original model predictors") # plotmo.x.default cannot get the x matrix
se <- ozone1
a <- earth(se[,2:10], se[,1], degree=2, keepxy=TRUE)
se$vh <- NULL # vh is unused (but plotmo still needs it --- why?)
printh(summary(a))
plotmo(a, trace=0, caption="getdata earth test9")

# test the way plotmo gets the data with lm
se <- ozone1
a <- lm(O3 ~ ., data=se)
printh(summary(a))
plotmo(a, trace=0, caption="getdata lm test1")
a <- lm(O3 ~ ., data=se, x=1)
printh(summary(a))
plotmo(a, trace=0, caption="getdata lm test2")
a <- lm(O3 ~ ., data=se, y=1)
printh(summary(a))
plotmo(a, trace=0, caption="getdata lm test3")
a <- lm(O3 ~ ., data=se, x=1, y=1)
printh(summary(a))
plotmo(a, trace=0, caption="getdata lm test3")
a <- lm(O3 ~ ., data=se, x=0, y=1, model=F)
se <- 99
expect.err(try(plotmo(a, trace=0, caption="getdata lm test4")), "cannot get the original model predictors")
se <- ozone1
a <- lm(O3 ~ ., data=se, x=1, y=1)
se <- 77
printh(summary(a))
plotmo(a, trace=0, caption="getdata lm test5")
se <- ozone1
a <- lm(O3 ~ ., data=se, model=F)
se$wind <- NULL
expect.err(try(plotmo(a, trace=0, caption="getdata lm test6")), "cannot get the original model predictors")

cat("test fixed.point warning in print.summary.earth\n")
options(digits=3)
et <- etitanic
et$age <- 1000 * et$age
a <- earth(survived~., data=et)
print(summary(a))
print(summary(a, fixed.point=FALSE))
options(digits=7) # back to default

cat("--- summary earth with new data ----------------------\n")
a.trees <- earth(Volume~., data=trees)
cat("summary(a.trees, newdata=trees)\n")
print(summary(a.trees, newdata=trees))
cat("summary(a.trees, newdata=trees[1:5,])\n")
a.trees.summary <- print(summary(a.trees, newdata=trees[1:5,]))

a.xy.trees <- earth(trees[,1:2], trees[,3])
cat("summary(a.xy.trees, newdata=trees[1:5,])\n")
a.xy.trees.summary <- print(summary(a.xy.trees, newdata=trees[1:5,]))
stopifnot(a.xy.trees.summary$newrsq == a.trees.summary$newrsq)

a.xy1.trees <- earth(trees[,1:2], trees$Volume)
cat("summary(a.xy1.trees, newdata=trees[1:5,])\n")
a.xy1.trees.summary <- print(summary(a.xy1.trees, newdata=trees[1:5,]))
stopifnot(a.xy1.trees.summary$newrsq == a.trees.summary$newrsq)

cat("--- ../../tests/test.earth.R -------------------------\n")

options(options.old)
source("../../tests/test.earth.R")

cat("--- check that spurious warn gone: non-integer #successes in a binomial glm ---\n")

library(segmented) # for down data
data(down)
fit.e <- earth(cases/births~age, data=down, weights=down$births, glm=list(family="binomial"))
summary(fit.e)

# test nk=1, 2, and 3
cat("nk=1\n")
par(mfrow = c(2, 2), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
a.nk1 <- earth(Volume~., data=trees, nk=1)
plot(a.nk1, which=1, main="nk=1")
print(a.nk1)
cat("nk=2\n")
a.nk2 <- earth(Volume~., data=trees, nk=2)
print(summary(a.nk2))
plot(a.nk2, which=1, main="nk=2")
cat("nk=3\n")
a.nk3 <- earth(Volume~., data=trees, nk=3)
plot(a.nk3, which=1, main="nk=3")

cat("\ntest model.matrix.earth\n")

check.model.matrix <- function(msg, xnew, bx1, bx2)
{
    cat("check.model.matrix", msg, ":\n")
    print(xnew)
    if(!identical(bx1, bx2)) {
       cat("\nnot identical\n")
       cat(deparse(substitute(bx1)), ":\n", sep="")
       print(bx1)
       cat(deparse(substitute(bx2)), ":\n", sep="")
       print(bx2)
       stop("check.model.matrix ", msg, " failed")
    }
}

data(trees)
earth.trees.formula <- earth(Volume ~ ., data=trees, subset=1:20)
bx <- model.matrix(earth.trees.formula)
check.model.matrix("earth.trees.formula formula 1", NULL, bx, earth.trees.formula$bx)

# nprune so only Girth is used, not Height
earth.girth.formula <- earth(Volume ~ ., data=trees, nprune=3)

# model.matrix where xnew is a data.frame

xnew <- trees[,1:2]
bx <- model.matrix(earth.girth.formula, xnew)
lm.mod <- lm(trees$Volume ~ bx[,-1])  # -1 to drop intercept
stopifnot(coef(earth.girth.formula) == coef(lm.mod))

colnames(xnew) <- NULL
bx <- model.matrix(earth.girth.formula, xnew)
lm.mod2 <- lm(trees$Volume ~ bx[,-1])
stopifnot(coef(earth.girth.formula) == coef(lm.mod2))

xnew <- data.frame(Girth=c(8.3, 8.6), Height=c(70, 65))
bx <- model.matrix(earth.girth.formula, xnew)
check.model.matrix("earth.girth.formula formula 2", xnew, bx, earth.girth.formula$bx[1:2,])

# test what happens when variables are missing
predict.girth.height <- predict(earth.girth.formula, xnew)
predict.girth  <- predict(earth.girth.formula, newdata=data.frame(Girth=c(8.3, 8.6)))
stopifnot(all.equal(predict.girth.height, predict.girth))
predict.height <- predict(earth.girth.formula, newdata=data.frame(Height=c(70, 65)))
stopifnot(all(is.na(predict.height)))

xnew <- trees[1:2,]
bx <- model.matrix(earth.girth.formula, xnew)
check.model.matrix("earth.girth.formula formula 3", xnew, bx, earth.girth.formula$bx[1:2,])

xnew <- trees[1:2,1:2] # exclude Volume column
colnames(xnew) <- NULL
bx <- model.matrix(earth.girth.formula, xnew)
check.model.matrix("earth.girth.formula formula 4", xnew, bx, earth.girth.formula$bx[1:2,])

xnew <- trees[1:2,3:1] # change order of columns
bx <- model.matrix(earth.girth.formula, xnew)
check.model.matrix("earth.girth.formula formula 5", xnew, bx, earth.girth.formula$bx[1:2,])

xnew <- trees[1:2,1,drop=FALSE]  # include only Girth
bx <- model.matrix(earth.girth.formula, xnew)
check.model.matrix("earth.girth.formula formula 6", xnew, bx, earth.girth.formula$bx[1:2,])

xnew <- trees[1,2:1]
bx <- model.matrix(earth.girth.formula, xnew)
check.model.matrix("earth.girth.formula formula 7", xnew, bx, earth.girth.formula$bx[1,,drop=FALSE])

xnew <- trees[1,1:2]
names(xnew) <- NULL
bx <- model.matrix(earth.girth.formula, xnew)
check.model.matrix("earth.girth.formula formula 8", xnew, bx, earth.girth.formula$bx[1,,drop=FALSE])

# model.matrix where xnew is a matrix (same as above code but with as.matrix)

xnew <- as.matrix(trees[,1:2])
bx <- model.matrix(earth.girth.formula, xnew)
lm.mod <- lm(trees$Volume ~ bx[,-1])  # -1 to drop intercept
stopifnot(coef(earth.girth.formula) == coef(lm.mod))

colnames(xnew) <- NULL
bx <- model.matrix(earth.girth.formula, xnew)
lm.mod2 <- lm(trees$Volume ~ bx[,-1])
stopifnot(coef(earth.girth.formula) == coef(lm.mod2))

xnew <- as.matrix(data.frame(Girth=c(8.3, 8.6), Height=c(70, 65)))
bx <- model.matrix(earth.girth.formula, xnew)
check.model.matrix("earth.girth.formula formula 9", xnew, bx, earth.girth.formula$bx[1:2,])

xnew <- as.matrix(trees[1:2,])
bx <- model.matrix(earth.girth.formula, xnew)
check.model.matrix("earth.girth.formula formula 10", xnew, bx, earth.girth.formula$bx[1:2,])

xnew <- as.matrix(trees[1:2,1:2]) # exclude Volume column
colnames(xnew) <- NULL
bx <- model.matrix(earth.girth.formula, xnew)
check.model.matrix("earth.girth.formula formula 11", xnew, bx, earth.girth.formula$bx[1:2,])

xnew <- as.matrix(trees[1:2,3:1]) # change order of columns
bx <- model.matrix(earth.girth.formula, xnew)
check.model.matrix("earth.girth.formula formula 12", xnew, bx, earth.girth.formula$bx[1:2,])

xnew <- as.matrix(trees[1:2,1,drop=FALSE]) # include only Girth
bx <- model.matrix(earth.girth.formula, xnew)
check.model.matrix("earth.girth.formula formula 13", xnew, bx, earth.girth.formula$bx[1:2,])

xnew <- as.matrix(trees[1,2:1])
bx <- model.matrix(earth.girth.formula, xnew, trace=2)
check.model.matrix("earth.girth.formula formula 14", xnew, bx, earth.girth.formula$bx[1,,drop=FALSE])

xnew <- as.matrix(trees[3,1:2])
names(xnew) <- NULL
bx <- model.matrix(earth.girth.formula, xnew)
check.model.matrix("earth.girth.formula formula 15", xnew, bx, earth.girth.formula$bx[3,,drop=FALSE])

#--- model.matrix with an x,y model

data(trees)
earth.trees.xy.subset <- earth(trees[,1:2], trees[,3], subset=1:20)
bx <- model.matrix(earth.trees.xy.subset)
check.model.matrix("earth.trees.xy.subset x,y 1", NULL, bx, earth.trees.xy.subset$bx)

# nprune so only Girth is used, not Height
earth.girth.xy <- earth(trees[,1:2], trees[,3], nprune=3)

# model.matrix where xnew is a data.frame

xnew <- trees[,1:2]
bx <- model.matrix(earth.girth.xy, xnew)
lm.mod <- lm(trees$Volume ~ bx[,-1])  # -1 to drop intercept
stopifnot(coef(earth.girth.xy) == coef(lm.mod))

colnames(xnew) <- NULL
bx <- model.matrix(earth.girth.xy, xnew)
lm.mod2 <- lm(trees$Volume ~ bx[,-1])
stopifnot(coef(earth.girth.xy) == coef(lm.mod2))

xnew <- data.frame(Girth=c(8.3, 8.6), Height=c(70, 65))
bx <- model.matrix(earth.girth.xy, xnew)
check.model.matrix("earth.girth.xy x,y 2", xnew, bx, earth.girth.xy$bx[1:2,])

# test what happens when variables are missing
predict.girth.height <- predict(earth.girth.xy, xnew)
predict.girth  <- predict(earth.girth.xy, newdata=data.frame(Girth=c(8.3, 8.6)))
stopifnot(all.equal(predict.girth.height, predict.girth))
predict.height <- predict(earth.girth.xy, newdata=data.frame(Height=c(70, 65)))
stopifnot(all(is.na(predict.height)))

xnew <- trees[1:2,]
bx <- model.matrix(earth.girth.xy, xnew)
check.model.matrix("earth.girth.xy x,y 3", xnew, bx, earth.girth.xy$bx[1:2,])

xnew <- trees[1:2,1:2] # exclude Volume column
colnames(xnew) <- NULL
bx <- model.matrix(earth.girth.xy, xnew)
check.model.matrix("earth.girth.xy x,y 4", xnew, bx, earth.girth.xy$bx[1:2,])

# # TODO fails
# xnew <- trees[1:2,3:1] # change order of columns
# bx <- model.matrix(earth.girth.xy, xnew)
# check.model.matrix("earth.girth.xy x,y 5", xnew, bx, earth.girth.xy$bx[1:2,])

xnew <- trees[1:2,1,drop=FALSE]  # include only Girth
bx <- model.matrix(earth.girth.xy, xnew)
check.model.matrix("earth.girth.xy x,y 6", xnew, bx, earth.girth.xy$bx[1:2,])

xnew <- trees[1,2:1]
bx <- model.matrix(earth.girth.xy, xnew)
check.model.matrix("earth.girth.xy x,y 7", xnew, bx, earth.girth.xy$bx[1,,drop=FALSE])

xnew <- trees[1,1:2]
names(xnew) <- NULL
bx <- model.matrix(earth.girth.xy, xnew)
check.model.matrix("earth.girth.xy x,y 8", xnew, bx, earth.girth.xy$bx[1,,drop=FALSE])

# model.matrix where xnew is a matrix (same as above code but with as.matrix)

xnew <- as.matrix(trees[,1:2])
bx <- model.matrix(earth.girth.xy, xnew)
lm.mod <- lm(trees$Volume ~ bx[,-1])  # -1 to drop intercept
stopifnot(coef(earth.girth.xy) == coef(lm.mod))

colnames(xnew) <- NULL
bx <- model.matrix(earth.girth.xy, xnew)
lm.mod2 <- lm(trees$Volume ~ bx[,-1])
stopifnot(coef(earth.girth.xy) == coef(lm.mod2))

xnew <- as.matrix(data.frame(Girth=c(8.3, 8.6), Height=c(70, 65)))
bx <- model.matrix(earth.girth.xy, xnew)
check.model.matrix("earth.girth.xy x,y 9", xnew, bx, earth.girth.xy$bx[1:2,])

xnew <- as.matrix(trees[1:2,])
bx <- model.matrix(earth.girth.xy, xnew)
check.model.matrix("earth.girth.xy x,y 10", xnew, bx, earth.girth.xy$bx[1:2,])

xnew <- as.matrix(trees[1:2,1:2]) # exclude Volume column
colnames(xnew) <- NULL
bx <- model.matrix(earth.girth.xy, xnew)
check.model.matrix("earth.girth.xy x,y 11", xnew, bx, earth.girth.xy$bx[1:2,])

# # TODO fails
# xnew <- as.matrix(trees[1:2,3:1]) # change order of columns
# bx <- model.matrix(earth.girth.xy, xnew)
# check.model.matrix("earth.girth.xy x,y 12", xnew, bx, earth.girth.xy$bx[1:2,])

xnew <- as.matrix(trees[1:2,1,drop=FALSE]) # include only Girth
bx <- model.matrix(earth.girth.xy, xnew)
check.model.matrix("earth.girth.xy x,y 13", xnew, bx, earth.girth.xy$bx[1:2,])

xnew <- as.matrix(trees[1,2:1])
bx <- model.matrix(earth.girth.xy, xnew, trace=2)
check.model.matrix("earth.girth.xy x,y 14", xnew, bx, earth.girth.xy$bx[1,,drop=FALSE])

xnew <- as.matrix(trees[3,1:2])
names(xnew) <- NULL
bx <- model.matrix(earth.girth.xy, xnew)
check.model.matrix("earth.girth.xy x,y 15", xnew, bx, earth.girth.xy$bx[3,,drop=FALSE])

cat("--- compare backward, none, exhaustive, forward, seqrep ---------------------\n")
data(ozone1)
oz <- ozone1[1:50,]
cat("--mod.none\n")
mod.none <- earth(O3~., data=oz, degree=2, trace=4, pmethod="none")
print(summary(mod.none))
cat("--mod.backward\n")
mod.backward <- earth(O3~., data=oz, degree=2, trace=4, pmethod="backward")
print(summary(mod.backward))
cat("--mod.forward\n")
mod.forward <- earth(O3~., data=oz, degree=2, trace=4, pmethod="forward")
print(summary(mod.forward))
cat("--mod.exhaustive\n")
mod.exhaustive <- earth(O3~., data=oz, degree=2, trace=4, pmethod="exhaustive")
print(summary(mod.exhaustive))
cat("--mod.seqrep\n")
mod.seqrep <- earth(O3~., data=oz, degree=2, trace=4, pmethod="seqrep")
print(summary(mod.seqrep))
tab <- data.frame(pmethod=c("none", "backward", "forward", "exhaustive", "seqrep"),
                  grsq=c(mod.none$grsq,
                         mod.backward$grsq,
                         mod.forward$grsq,
                         mod.exhaustive$grsq,
                         mod.seqrep$grsq),
                  nterms=c(length(mod.none$selected.terms),
                           length(mod.backward$selected.terms),
                           length(mod.forward$selected.terms),
                           length(mod.exhaustive$selected.terms),
                           length(mod.seqrep$selected.terms)))
cat("\n")
print(tab)

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}
