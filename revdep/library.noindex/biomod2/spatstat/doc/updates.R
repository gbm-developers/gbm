### R code from vignette source 'updates.Rnw'

###################################################
### code chunk number 1: updates.Rnw:20-24
###################################################
library(spatstat)
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 2: updates.Rnw:38-74
###################################################
readSizeTable <- function(fname) {
  if(is.null(fname) || !file.exists(fname)) return(NULL)
  a <- read.table(fname, header=TRUE)
  a$date <- as.Date(a$date)
  return(a)
}
getSizeTable <- function(packagename="spatstat", tablename="packagesizes.txt") {
  fname <- system.file("doc", tablename, package=packagename)
  readSizeTable(fname)
}
counts <- c("nhelpfiles", "nobjects", "ndatasets", "Rlines", "srclines")
mergeSizeTables <- function(a, b) {
  if(is.null(b)) return(a)
  for(i in seq_len(nrow(a))) {
    j <- which(b$date <= a$date[i])
    if(length(j) > 0) 
      a[i,counts] <- a[i,counts] + b[max(j), counts]
  }
  return(a)
}
z <- getSizeTable()
zutils <- getSizeTable("spatstat.utils")
zdata <- getSizeTable("spatstat.data")
zlocal <- getSizeTable("spatstat", "spatstatlocalsize.txt")
z <- mergeSizeTables(z, zutils)
z <- mergeSizeTables(z, zdata)
z <- mergeSizeTables(z, zlocal)
#
currentcount <- z[nrow(z), counts]
bookcount    <- z[z$version == "1.42-0", counts]
changes <- currentcount - bookcount
newobj <- changes[["nobjects"]]
newdat <- changes[["ndatasets"]] + 1  # counting rule doesn't detect redwood3
newcode  <- changes[["Rlines"]] + changes[["srclines"]]
bookcode <- bookcount[["Rlines"]] + bookcount[["srclines"]]
growth <- signif((100 * newcode)/bookcode, digits=2)


###################################################
### code chunk number 3: updates.Rnw:91-96
###################################################
options(SweaveHooks=list(fig=function() par(mar=0.2+c(2,4,2,0))))
Plot <- function(fmla, ..., dat=z) {
  yvals <- eval(as.expression(fmla[[2]]), envir=dat)
  plot(fmla, ..., data=dat, type="l", xlab="", lwd=2, ylim=c(0, max(yvals)))
}


###################################################
### code chunk number 4: updates.Rnw:102-107
###################################################
getOption("SweaveHooks")[["fig"]]()
Plot((Rlines + srclines)/1000 ~ date, ylab="Lines of code (x 1000)", 
     main="Spatstat growth")
lines(srclines/1000 ~ date, data=z)
text(as.Date("2015-01-01"), 9.5, "C code")
text(as.Date("2015-01-01"), 60, "R code")


