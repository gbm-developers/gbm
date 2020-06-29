pkgname <- "crimelinkage"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('crimelinkage')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("clusterPath")
### * clusterPath

flush(stderr()); flush(stdout())

### Name: clusterPath
### Title: Follows path of one crime up a dendrogram
### Aliases: clusterPath

### ** Examples

# See vignette: "Crime Series Identification and Clustering" for usage.



cleanEx()
nameEx("compareCrimes")
### * compareCrimes

flush(stderr()); flush(stdout())

### Name: compareCrimes
### Title: Creates evidence variables by calculating 'distance' between
###   crime pairs
### Aliases: compareCrimes

### ** Examples

data(crimes)
 pairs = t(combn(crimes$crimeID[1:4],m=2))   # make some crime pairs

 varlist = list(
   spatial = c("X", "Y"),
   temporal = c("DT.FROM","DT.TO"),
   categorical = c("MO1",  "MO2", "MO3"))    # crime variables list

 compareCrimes(pairs,crimes,varlist,binary=TRUE)



cleanEx()
nameEx("crimeClust_bayes")
### * crimeClust_bayes

flush(stderr()); flush(stdout())

### Name: crimeClust_bayes
### Title: Bayesian model-based partially-supervised clustering for crime
###   series identification
### Aliases: crimeClust_bayes

### ** Examples

# Toy dataset with 12 crimes and three criminals.

 # Make IDs: Criminal 1 committed crimes 1-4, etc.
 id <- c(1,1,1,1,
         2,2,2,2,
                 3,3,3,3)

 # spatial locations of the crimes:
 s <- c(0.8,0.9,1.1,1.2,
        1.8,1.9,2.1,2.2,
        2.8,2.9,3.1,3.2)
 s <- cbind(0,s)

 # Categorical crime features, say mode of entry (1=door, 2=other) and
 # type of residence (1=apartment, 2=other)
 Mode <- c(1,1,1,1,  #Different distribution by criminal
           1,2,1,2,
           2,2,2,2)
 Type <- c(1,2,1,2,  #Same distribution for all criminals
           1,2,1,2,
           1,2,1,2)
 Xcat <- cbind(Mode,Type)

 # Times of the crimes
 t <- c(1,2,3,4,
        2,3,4,5,
        3,4,5,6)

 # Now let's pretend we don't know the criminal for crimes 1, 4, 6, 8, and 12.
 id <- c(NA,1,1,NA,2,NA,2,NA,3,3,3,NA)

 # Fit the model (nb: use much larger iters and burn on real problem)
 fit <- crimeClust_bayes(crimeID=id, spatial=s, t1=t,t2=t, Xcat=Xcat,
                   maxcriminals=12,iters=500,burn=100,update=100)

 # Plot the posterior probability matrix that each pair of crimes was
 # committed by the same criminal:
 if(require(fields,quietly=TRUE)){
 fields::image.plot(1:12,1:12,fit$p.equal,
            xlab="Crime",ylab="Crime",
            main="Probability crimes are from the same criminal")
 }

 # Extract the crimes with the largest posterior probability
 bayesPairs(fit$p.equal)
 bayesProb(fit$p.equal[1,])



cleanEx()
nameEx("crimeClust_hier")
### * crimeClust_hier

flush(stderr()); flush(stdout())

### Name: crimeClust_hier
### Title: Agglomerative Hierarchical Crime Series Clustering
### Aliases: crimeClust_hier

### ** Examples

data(crimes)
 #- cluster the first 10 crime incidents
 crimedata = crimes[1:10,]
 varlist = list(spatial = c("X", "Y"), temporal = c("DT.FROM","DT.TO"),
     categorical = c("MO1",  "MO2", "MO3"))
 estimateBF <- function(X) rnorm(NROW(X))   # random estimation of log Bayes Factor
 HC = crimeClust_hier(crimedata,varlist,estimateBF)
 plot_hcc(HC,yticks=-2:2)

 # See vignette: "Crime Series Identification and Clustering" for more examples.



cleanEx()
nameEx("crimes")
### * crimes

flush(stderr()); flush(stdout())

### Name: crimes
### Title: Ficticious dataset of crime events
### Aliases: crimes
### Keywords: datasets

### ** Examples

head(crimes)



cleanEx()
nameEx("getBF")
### * getBF

flush(stderr()); flush(stdout())

### Name: getBF
### Title: Estimates the bayes factor for continous and categorical
###   predictors.
### Aliases: getBF

### ** Examples

# See vignette: "Statistical Methods for Crime Series Linkage" for usage.



cleanEx()
nameEx("getCrimeSeries")
### * getCrimeSeries

flush(stderr()); flush(stdout())

### Name: getCrimeSeries
### Title: Generate a list of offenders and their associated crime series.
### Aliases: getCrimeSeries

### ** Examples

data(offenders)

 getCrimeSeries("O:40",offenders)
 getCrimeSeries(c("O:40","O:3"),offenders)  # list of crime series from multiple offenders



cleanEx()
nameEx("getCrimes")
### * getCrimes

flush(stderr()); flush(stdout())

### Name: getCrimes
### Title: Generate a list of crimes for a specific offender
### Aliases: getCrimes

### ** Examples

data(crimes)
 data(offenders)

 getCrimes("O:40",crimes,offenders)



cleanEx()
nameEx("getCriminals")
### * getCriminals

flush(stderr()); flush(stdout())

### Name: getCriminals
### Title: Lookup the offenders responsible for a set of solved crimes
### Aliases: getCriminals

### ** Examples

data(offenders)

 getCriminals("C:1",offenders)

 getCriminals("C:78",offenders)                      # shows co-offenders

 getCriminals(c("C:26","C:78","85","110"),offenders) # all offenders from a crime series



cleanEx()
nameEx("getROC")
### * getROC

flush(stderr()); flush(stdout())

### Name: getROC
### Title: Cacluate ROC like metrics.
### Aliases: getROC

### ** Examples

f = 1:10
y = rep(0:1,length=10)
getROC(f,y)



cleanEx()
nameEx("linkage")
### * linkage

flush(stderr()); flush(stdout())

### Name: linkage
### Title: Hierarchical Based Linkage
### Aliases: linkage

### ** Examples

# See vignette: "Crime Series Identification and Clustering" for usage.



cleanEx()
nameEx("makeGroups")
### * makeGroups

flush(stderr()); flush(stdout())

### Name: makeGroups
### Title: Generates crime groups from crime series data
### Aliases: makeGroups

### ** Examples

data(crimes)
 data(offenders)
 seriesData = makeSeriesData(crimedata=crimes,offenderTable=offenders)
 groups = makeGroups(seriesData,method=1)
 head(groups,10)



cleanEx()
nameEx("makePairs")
### * makePairs

flush(stderr()); flush(stdout())

### Name: makePairs
### Title: Generates indices of linked and unlinked crime pairs (with
###   weights)
### Aliases: makeLinked makePairs makeUnlinked

### ** Examples

data(crimes)
 data(offenders)
 seriesData = makeSeriesData(crimedata=crimes,offenderTable=offenders)
 allPairs = makePairs(seriesData,thres=365,m=40)



cleanEx()
nameEx("makeSeriesData")
### * makeSeriesData

flush(stderr()); flush(stdout())

### Name: makeSeriesData
### Title: Make crime series data
### Aliases: makeSeriesData

### ** Examples

data(crimes)
 data(offenders)

 seriesData = makeSeriesData(crimedata=crimes,offenderTable=offenders)
 head(seriesData)

 nCrimes = table(seriesData$offenderID)  # length of each crime series
 table(nCrimes)                  # distribution of crime series length
 mean(nCrimes>1)                 # proportion of offenders with multiple crimes

 nCO = table(seriesData$crimeID) # number of co-offenders per crime
 table(nCO)                      # distribution of number of co-offenders
 mean(nCO>1)                     # proportion of crimes with multiple co-offenders



cleanEx()
nameEx("naiveBayes")
### * naiveBayes

flush(stderr()); flush(stdout())

### Name: naiveBayes
### Title: Naive bayes classifier using histograms and shrinkage
### Aliases: naiveBayes naiveBayes.fit

### ** Examples

# See vignette: "Statistical Methods for Crime Series Linkage" for usage.



cleanEx()
nameEx("offenders")
### * offenders

flush(stderr()); flush(stdout())

### Name: offenders
### Title: Ficticious offender data
### Aliases: offenders
### Keywords: datasets

### ** Examples

head(offenders)



cleanEx()
nameEx("plot.naiveBayes")
### * plot.naiveBayes

flush(stderr()); flush(stdout())

### Name: plot.naiveBayes
### Title: Plots for Naive Bayes Model
### Aliases: plot.naiveBayes

### ** Examples

# See vignette: "Statistical Methods for Crime Series Linkage" for usage.



cleanEx()
nameEx("plotBF")
### * plotBF

flush(stderr()); flush(stdout())

### Name: plotBF
### Title: plots 1D bayes factor
### Aliases: plotBF

### ** Examples

# See vignette: "Statistical Methods for Crime Series Linkage" for usage.



cleanEx()
nameEx("plot_hcc")
### * plot_hcc

flush(stderr()); flush(stdout())

### Name: plot_hcc
### Title: Plot a hierarchical crime clustering object
### Aliases: plot_hcc

### ** Examples

# See vignette: "Crime Series Identification and Clustering" for usage.



cleanEx()
nameEx("predict.naiveBayes")
### * predict.naiveBayes

flush(stderr()); flush(stdout())

### Name: predict.naiveBayes
### Title: Generate prediction (sum of log bayes factors) from a
###   'naiveBayes' object
### Aliases: predict.naiveBayes

### ** Examples

# See vignette: "Statistical Methods for Crime Series Linkage" for usage.



cleanEx()
nameEx("predictBF")
### * predictBF

flush(stderr()); flush(stdout())

### Name: predictBF
### Title: Generate prediction of a component bayes factor
### Aliases: predictBF

### ** Examples

# See vignette: "Statistical Methods for Crime Series Linkage" for usage.



cleanEx()
nameEx("seriesID")
### * seriesID

flush(stderr()); flush(stdout())

### Name: seriesID
### Title: Crime series identification
### Aliases: seriesID

### ** Examples

# See vignette: "Crime Series Identification and Clustering" for usage.



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
