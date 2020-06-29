pkgname <- "simulatorZ"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('simulatorZ')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("cvSubsets")
### * cvSubsets

flush(stderr()); flush(stdout())

### Name: cvSubsets
### Title: cvSubsets
### Aliases: cvSubsets

### ** Examples



library(curatedOvarianData)


data(E.MTAB.386_eset)





id <- cvSubsets(E.MTAB.386_eset, 3)


subsets <- lapply(1:3, function(i){E.MTAB.386_eset[1:10, id[[i]]]})


sapply(subsets, dim)


rm(subsets)





## Number of observations in the set does not need to be a multiple of


## the fold parameter


id2 <- cvSubsets(E.MTAB.386_eset, 5)


subsets <- lapply(1:5, function(j){E.MTAB.386_eset[1:10, id2[[j]]]})


sapply(subsets, dim)


rm(subsets)





cleanEx()
nameEx("funCV")
### * funCV

flush(stderr()); flush(stdout())

### Name: funCV
### Title: funCV
### Aliases: funCV

### ** Examples



library(curatedOvarianData)


library(GenomicRanges)


set.seed(8)


data( E.MTAB.386_eset )


eset <- E.MTAB.386_eset[1:100, 1:30]


rm(E.MTAB.386_eset)





time <- eset$days_to_death


cens.chr <- eset$vital_status


cens <- rep(0, length(cens.chr))


cens[cens.chr=="living"] <- 1


y <- Surv(time, cens)  


y1 <- cbind(time, cens)





nrows <- 200; ncols <- 6


counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)


rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),


                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),


                     strand=sample(c("+", "-"), 200, TRUE))


colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),


                     row.names=LETTERS[1:6])


sset <- SummarizedExperiment(assays=SimpleList(counts=counts),


                             rowRanges=rowRanges, colData=colData)


time <- c(1588,1929,1813,1542,1830,1775)  


cens <- c(1,0,1,1,1,1)


y.vars <- Surv(time, cens)





funCV(eset, 3, y)


funCV(exprs(eset), 3, y1)


funCV(sset, 3, y.vars)


## any training function will do as long as it takes the gene expression matrix X


## and response variable y(matrix, data.frame or Surv object) as parameters, and


## return the coefficients as its value





cleanEx()
nameEx("geneFilter")
### * geneFilter

flush(stderr()); flush(stdout())

### Name: geneFilter
### Title: geneFilter
### Aliases: geneFilter

### ** Examples



set.seed(8)


library(curatedOvarianData)


library(GenomicRanges)


data(GSE17260_eset)


data(E.MTAB.386_eset)


data(GSE14764_eset)


## to save time, we take a small subset from each dataset


esets.list <- list(GSE17260=GSE17260_eset[1:50, 1:10], 


                   E.MTAB.386=E.MTAB.386_eset[1:50, 1:10], 


                   GSE14764=GSE14764_eset[1:50, 1:10])


rm(E.MTAB.386_eset, GSE14764_eset, GSE17260_eset)





result.set <- geneFilter(esets.list, 0.1)


dim(result.set[[1]])





## as we cannot calculate correlation with one set, this function just 


## delivers the same set if esets has length 1


result.oneset <- geneFilter(esets.list[1])


dim(result.oneset[[1]])








## Support matrices


X.list <- lapply(esets.list, function(eset){


  return(exprs(eset)) ## Columns represent samples!


})


result.set <- geneFilter(X.list, 0.1)


dim(result.set[[1]])





## Support RangedSummarizedExperiment


nrows <- 200; ncols <- 6


counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)


rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),


                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),


                     strand=sample(c("+", "-"), 200, TRUE))


colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),


                     row.names=LETTERS[1:6])


sset <- SummarizedExperiment(assays=SimpleList(counts=counts),


                             rowRanges=rowRanges, colData=colData)


s.list <- list(sset, sset)


result.set <- geneFilter(s.list, 0.9) 


## the same set should resemble each other, no genes filtered


dim(assay(result.set[[1]]))





cleanEx()
nameEx("getTrueModel")
### * getTrueModel

flush(stderr()); flush(stdout())

### Name: getTrueModel
### Title: getTrueModel
### Aliases: getTrueModel

### ** Examples



library(curatedOvarianData)


data(GSE14764_eset)


data(E.MTAB.386_eset)


esets.list <- list(GSE14764=GSE14764_eset[1:500, 1:20], 


                   E.MTAB.386=E.MTAB.386_eset[1:500, 1:20])


rm(E.MTAB.386_eset, GSE14764_eset)





## simulate on multiple ExpressionSets


set.seed(8) 





y.list <- lapply(esets.list, function(eset){


  time <- eset$days_to_death


  cens.chr <- eset$vital_status


  cens <- rep(0, length(cens.chr))


  cens[cens.chr=="living"] <- 1


  return(Surv(time, cens))


})


   


res1 <- getTrueModel(esets.list, y.list, 100)


## Get true model from one set


res2 <- getTrueModel(esets.list[1], y.list[1], 100)


names(res2)


res2$lp


## note that y.list[1] cannot be replaced by y.list[[1]]





cleanEx()
nameEx("masomenos")
### * masomenos

flush(stderr()); flush(stdout())

### Name: masomenos
### Title: masomenos
### Aliases: masomenos

### ** Examples



set.seed(8)


library(curatedOvarianData)


data( E.MTAB.386_eset )


eset <- E.MTAB.386_eset[1:100, 1:30]


rm(E.MTAB.386_eset)





X <- t(exprs(eset))





time <- eset$days_to_death


cens <- sample(0:1, 30, replace=TRUE)


y <- Surv(time, cens)





beta <- masomenos(X=X, y=y)


beta





cleanEx()
nameEx("plusMinus")
### * plusMinus

flush(stderr()); flush(stdout())

### Name: plusMinus
### Title: plusMinus
### Aliases: plusMinus

### ** Examples



set.seed(8)


library(curatedOvarianData)


data( E.MTAB.386_eset )


eset <- E.MTAB.386_eset[1:100, 1:30]


rm(E.MTAB.386_eset)





X <- t(exprs(eset))





time <- eset$days_to_death


cens <- sample(0:1, 30, replace=TRUE)


y <- Surv(time, cens)





beta <- plusMinus(X, y)


beta





cleanEx()
nameEx("rowCoxTests")
### * rowCoxTests

flush(stderr()); flush(stdout())

### Name: rowCoxTests
### Title: rowCoxTests
### Aliases: rowCoxTests

### ** Examples



#test


##regressor-matrix (gene expressions)


X<-matrix(rnorm(1e6),nrow=10000)


#seed


set.seed(123)


#times


time<-rnorm(n=ncol(X),mean=100)


#censoring(1->death)


status<-rbinom(n=ncol(X),size=1, prob=0.8)





##survival object


y<-Surv(time,status)





## Do 10,000 Cox regressions:


system.time(output <- rowCoxTests(X=X,y=y, option="fast"))





cleanEx()
nameEx("simBootstrap")
### * simBootstrap

flush(stderr()); flush(stdout())

### Name: simBootstrap
### Title: simBootstrap
### Aliases: simBootstrap

### ** Examples



library(curatedOvarianData)


library(GenomicRanges)


data(E.MTAB.386_eset)


data(GSE14764_eset)


esets.list <- list(E.MTAB.386=E.MTAB.386_eset[1:200, 1:20], GSE14764=GSE14764_eset[1:200, 1:20])


rm(E.MTAB.386_eset, GSE14764_eset)





## simulate on multiple ExpressionSets


set.seed(8) 





y.list <- lapply(esets.list, function(eset){


  time <- eset$days_to_death


  cens.chr <- eset$vital_status


  cens <- rep(0, length(cens.chr))


  cens[cens.chr=="living"] <- 1


  return(Surv(time, cens))


})





simmodels <- simBootstrap(obj=esets.list, y.vars=y.list, 10, 100)


simmodels$obj.list[[1]]





# balance covariates


simmodels <- simBootstrap(obj=esets.list, y.vars=y.list, 10, 100,


                          balance.variables="tumorstage")


rm(esets.list, simmodels)





## Support RangedSummarizedExperiment


nrows <- 200; ncols <- 10


counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)


rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),


                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),


                     strand=sample(c("+", "-"), 200, TRUE))


colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 5),


                     row.names=LETTERS[1:10])


sset <- SummarizedExperiment(assays=SimpleList(counts=counts),


                             rowRanges=rowRanges, colData=colData)





s.list <- list(sset[,1:5], sset[,6:10])


time <- c(540, 527, 668, 587, 620, 540, 527, 668, 587, 620)


cens <- c(1, 0, 0, 1, 0, 1, 0, 0, 1, 0)


y.vars <- Surv(time, cens)


y.vars <- list(y.vars[1:5,],y.vars[1:5,])


simmodels <- simBootstrap(obj=s.list, y.vars=y.vars, 20, 100) 





cleanEx()
nameEx("simData")
### * simData

flush(stderr()); flush(stdout())

### Name: simData
### Title: simData
### Aliases: simData

### ** Examples



library(curatedOvarianData)


library(GenomicRanges)





data(E.MTAB.386_eset)


data(GSE14764_eset)


esets.list <- list(E.MTAB.386=E.MTAB.386_eset[1:100, 1:10], GSE14764=GSE14764_eset[1:100, 1:10])


rm(E.MTAB.386_eset, GSE14764_eset)





## simulate on multiple ExpressionSets


set.seed(8)


# one-step bootstrap: skip resampling set labels


simmodels <- simData(esets.list, 20, type="one-step")  


# two-step-non-parametric bootstrap


simmodels <- simData(esets.list, 10, type="two-steps")





## simulate one set


simmodels <- simData(list(esets.list[[1]]), 10, type="two-steps")





## balancing covariates


# single covariate


simmodels <- simData(list(esets.list[[1]]), 5, balance.variables="tumorstage")





# multiple covariates


simmodels <- simData(list(esets.list[[1]]), 5, 


                     balance.variables=c("tumorstage", "age_at_initial_pathologic_diagnosis"))  





## Support matrices


X.list <- lapply(esets.list, function(eset){


  return(exprs(eset))


})


simmodels <- simData(X.list, 20, type="two-steps")





## Support RangedSummarizedExperiment


nrows <- 200; ncols <- 6


counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)


rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),


                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),


                     strand=sample(c("+", "-"), 200, TRUE))


colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),


                     row.names=LETTERS[1:6])


sset <- SummarizedExperiment(assays=SimpleList(counts=counts),


                             rowRanges=rowRanges, colData=colData)





s.list <- list(sset[,1:3], sset[,4:6])


simmodels <- simData(s.list, 20, type="two-steps")





cleanEx()
nameEx("simTime")
### * simTime

flush(stderr()); flush(stdout())

### Name: simTime
### Title: simTime
### Aliases: simTime

### ** Examples



library(curatedOvarianData)


data(E.MTAB.386_eset)


data(GSE14764_eset)


esets.list <- list(E.MTAB.386=E.MTAB.386_eset[1:100, 1:20], GSE14764=GSE14764_eset[1:100, 1:20])


rm(E.MTAB.386_eset, GSE14764_eset)





## simulate on multiple ExpressionSets


set.seed(8) 





y.list <- lapply(esets.list, function(eset){


  time <- eset$days_to_death


  cens.chr <- eset$vital_status


  cens <- rep(0, length(cens.chr))


  cens[cens.chr=="living"] <- 1


  return(Surv(time, cens))


})





# To perform both parametric and non-parametric bootstrap, you can call simBootstrap()


# or, you can divide the steps into:


res <- getTrueModel(esets.list, y.list, 100)


simmodels <- simData(obj=esets.list, y.vars=y.list, n.samples=10)





# Then, use this function


simmodels <- simTime(simmodels=simmodels, original.yvars=y.list, result=res) 





# it also supports performing only the parametrc bootstrap step on a list of expressionsets


# but you need to construct the parameter by scratch


res <- getTrueModel(esets.list, y.list, 100)


setsID <- seq_along(esets.list)


indices <- list()


for(i in setsID){


  indices[[i]] <- seq_along(sampleNames(esets.list[[i]])) 


}


simmodels <- list(obj=esets.list, y.vars=y.list, indices=indices, setsID=setsID)





new.simmodels <- simTime(simmodels=simmodels, original.yvars=y.list, result=res)  





cleanEx()
nameEx("zmatrix")
### * zmatrix

flush(stderr()); flush(stdout())

### Name: zmatrix
### Title: zmatrix
### Aliases: zmatrix

### ** Examples



library(curatedOvarianData)


library(GenomicRanges)


data(E.MTAB.386_eset)


data(GSE14764_eset)


esets.list <- list(E.MTAB.386=E.MTAB.386_eset[1:100, 1:30], GSE14764=GSE14764_eset[1:100, 1:30])


rm(E.MTAB.386_eset, GSE14764_eset)





## simulate on multiple ExpressionSets


set.seed(8) 





y.list <- lapply(esets.list, function(eset){


  time <- eset$days_to_death


  cens.chr <- eset$vital_status


  cens <- rep(0, length(cens.chr))


  cens[cens.chr=="living"] <- 1


  return(Surv(time, cens))


})





# generate on original ExpressionSets


z <- zmatrix(esets.list, y.list, 3)





# generate on simulated ExpressionSets


simmodels <- simBootstrap(esets.list, y.list, 100, 100)


z <- zmatrix(simmodels$obj.list, simmodels$y.vars.list, 3)





# support matrix


X.list <- lapply(esets.list, function(eset){


  return(exprs(eset)) ### columns represent samples !!


}) 


z <- zmatrix(X.list, y.list, 3)





# support RangedSummarizedExperiment


nrows <- 200; ncols <- 6


counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)


rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),


                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),


                     strand=sample(c("+", "-"), 200, TRUE))


colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),


                     row.names=LETTERS[1:6])


sset <- SummarizedExperiment(assays=SimpleList(counts=counts),


                             rowRanges=rowRanges, colData=colData)





time <- sample(4500:4700, 6, replace=TRUE)


cens <- sample(0:1, 6, replace=TRUE)


y.vars <- Surv(time, cens)





z <- zmatrix(list(sset[,1:3], sset[,4:6]), list(y.vars[1:3,],y.vars[4:6,]), 3)





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
