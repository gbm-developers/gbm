### R code from vignette source 'getstart.Rnw'

###################################################
### code chunk number 1: getstart.Rnw:5-6
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 2: getstart.Rnw:25-32
###################################################
library(spatstat)
spatstat.options(image.colfun=function(n) { grey(seq(0,1,length=n)) })
sdate <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Date")
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 3: getstart.Rnw:56-58
###################################################
getOption("SweaveHooks")[["fig"]]()
data(redwood)
plot(redwood, pch=16, main="")


###################################################
### code chunk number 4: getstart.Rnw:80-82
###################################################
getOption("SweaveHooks")[["fig"]]()
data(longleaf)
plot(longleaf, main="")


###################################################
### code chunk number 5: getstart.Rnw:139-142
###################################################
data(finpines)
mypattern <- unmark(finpines)
mydata <- round(as.data.frame(finpines), 2)


###################################################
### code chunk number 6: getstart.Rnw:157-158 (eval = FALSE)
###################################################
## mydata <- read.csv("myfile.csv")


###################################################
### code chunk number 7: getstart.Rnw:169-170
###################################################
head(mydata)


###################################################
### code chunk number 8: getstart.Rnw:185-186 (eval = FALSE)
###################################################
##   mypattern <- ppp(mydata[,3], mydata[,7], c(100,200), c(10,90))


###################################################
### code chunk number 9: getstart.Rnw:189-190 (eval = FALSE)
###################################################
## ppp(x.coordinates, y.coordinates, x.range, y.range)


###################################################
### code chunk number 10: getstart.Rnw:199-200
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(mypattern)


###################################################
### code chunk number 11: getstart.Rnw:207-208 (eval = FALSE)
###################################################
## summary(mypattern)


###################################################
### code chunk number 12: getstart.Rnw:212-213
###################################################
options(SweaveHooks=list(fig=function() par(mar=rep(4,4)+0.1)))


###################################################
### code chunk number 13: getstart.Rnw:215-216
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(Kest(mypattern))


###################################################
### code chunk number 14: getstart.Rnw:222-223 (eval = FALSE)
###################################################
## plot(envelope(mypattern,Kest))


###################################################
### code chunk number 15: getstart.Rnw:225-226
###################################################
env <- envelope(mypattern,Kest, nsim=39)


###################################################
### code chunk number 16: getstart.Rnw:228-229
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(env, main="envelope(mypattern, Kest)")


###################################################
### code chunk number 17: getstart.Rnw:231-232
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 18: getstart.Rnw:238-239
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(density(mypattern))


###################################################
### code chunk number 19: getstart.Rnw:249-250 (eval = FALSE)
###################################################
## marks(mypattern) <- mydata[, c(5,9)]


###################################################
### code chunk number 20: getstart.Rnw:252-253
###################################################
mypattern <-finpines


###################################################
### code chunk number 21: getstart.Rnw:256-257 (eval = FALSE)
###################################################
## plot(Smooth(mypattern))


###################################################
### code chunk number 22: getstart.Rnw:260-261
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(Smooth(mypattern, sigma=1.2), main="Smooth(mypattern)")


