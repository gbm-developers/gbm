pkgname <- "featurefinder"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('featurefinder')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("addFeatures")
### * addFeatures

flush(stderr()); flush(stdout())

### Name: addFeatures
### Title: addFeatures
### Aliases: addFeatures
### Keywords: addFeatures

### ** Examples


require(featurefinder)
data(futuresdata)
data=futuresdata
data$SMIfactor=paste("smi",as.matrix(data$SMIfactor),sep="")
n=length(data$DAX)
nn=floor(length(data$DAX)/2)

# Can we predict the relative movement of DAX and SMI?
data$y=data$DAX*0 # initialise the target to 0
data$y[1:(n-1)]=((data$DAX[2:n])-(data$DAX[1:(n-1)]))/
  (data$DAX[1:(n-1)])-(data$SMI[2:n]-(data$SMI[1:(n-1)]))/(data$SMI[1:(n-1)])

# Fit a simple model
thismodel=lm(formula=y ~ .,data=data)
expected=predict(thismodel,data)
actual=data$y
residual=actual-expected
data=cbind(data,expected, actual, residual)

CSVPath=tempdir()
fcsv=paste(CSVPath,"/futuresdata.csv",sep="")
write.csv(data[(nn+1):(length(data$y)),],file=fcsv,row.names=FALSE)
exclusionVars="\"residual\",\"expected\", \"actual\",\"y\""
factorToNumericList=c()

# Now the dataset is prepared, try to find new features
tempDir=findFeatures(outputPath="NoPath", fcsv, exclusionVars,
factorToNumericList,                     
treeGenerationMinBucket=50,
treeSummaryMinBucket=20,
useSubDir=FALSE)  
         
newfeat1=((data$SMIfactor==0) & (data$CAC < 2253) & (data$CAC< 1998) & (data$CAC>=1882)) * 1.0
newfeat2=((data$SMIfactor==1) & (data$SMI < 7837) & (data$SMI >= 7499)) * 1.0
newfeatures=cbind(newfeat1, newfeat2) # create columns for the newly found features
datanew=cbind(data,newfeatures)
thismodel=lm(formula=y ~ .,data=datanew)
expectednew=predict(thismodel,datanew)

requireNamespace("Metrics")
OriginalRMSE = Metrics::rmse(data$y,expected)
NewRMSE = Metrics::rmse(data$y,expectednew)

print(paste("OriginalRMSE = ",OriginalRMSE))
print(paste("NewRMSE = ",NewRMSE))

# Append new features to a dataframe automatically
dataWithNewFeatures = addFeatures(df=data, path=tempDir, prefix="auto_")
head(df)



cleanEx()
nameEx("dat")
### * dat

flush(stderr()); flush(stdout())

### Name: dat
### Title: dat
### Aliases: dat
### Keywords: dat

### ** Examples

data(dat)
head(dat)



cleanEx()
nameEx("dat0")
### * dat0

flush(stderr()); flush(stdout())

### Name: dat0
### Title: dat0
### Aliases: dat0
### Keywords: dat0

### ** Examples

data(dat0)
head(dat0)



cleanEx()
nameEx("data")
### * data

flush(stderr()); flush(stdout())

### Name: data
### Title: data
### Aliases: data
### Keywords: data

### ** Examples

data(data)
head(data)



cleanEx()
nameEx("doAllFactors")
### * doAllFactors

flush(stderr()); flush(stdout())

### Name: doAllFactors
### Title: doAllFactors
### Aliases: doAllFactors
### Keywords: doAllFactors

### ** Examples

data(doAllFactors)
head(doAllFactors)



cleanEx()
nameEx("expr")
### * expr

flush(stderr()); flush(stdout())

### Name: expr
### Title: expr
### Aliases: expr
### Keywords: expr

### ** Examples

data(expr)
head(expr)



cleanEx()
nameEx("fileConn")
### * fileConn

flush(stderr()); flush(stdout())

### Name: fileConn
### Title: fileConn
### Aliases: fileConn
### Keywords: fileConn

### ** Examples

data(fileConn)
head(fileConn)



cleanEx()
nameEx("filename")
### * filename

flush(stderr()); flush(stdout())

### Name: filename
### Title: filename
### Aliases: filename
### Keywords: filename

### ** Examples

data(filename)
head(filename)



cleanEx()
nameEx("findFeatures")
### * findFeatures

flush(stderr()); flush(stdout())

### Name: findFeatures
### Title: findFeatures
### Aliases: findFeatures
### Keywords: findFeatures

### ** Examples


require(featurefinder)
data(futuresdata)
data=futuresdata
data$SMIfactor=paste("smi",as.matrix(data$SMIfactor),sep="")
n=length(data$DAX)
nn=floor(length(data$DAX)/2)

# Can we predict the relative movement of DAX and SMI?
data$y=data$DAX*0 # initialise the target to 0
data$y[1:(n-1)]=((data$DAX[2:n])-(data$DAX[1:(n-1)]))/
  (data$DAX[1:(n-1)])-(data$SMI[2:n]-(data$SMI[1:(n-1)]))/(data$SMI[1:(n-1)])

# Fit a simple model
thismodel=lm(formula=y ~ .,data=data)
expected=predict(thismodel,data)
actual=data$y
residual=actual-expected
data=cbind(data,expected, actual, residual)

CSVPath=tempdir()
fcsv=paste(CSVPath,"/futuresdata.csv",sep="")
write.csv(data[(nn+1):(length(data$y)),],file=fcsv,row.names=FALSE)
exclusionVars="\"residual\",\"expected\", \"actual\",\"y\""
factorToNumericList=c()

# Now the dataset is prepared, try to find new features
findFeatures(outputPath="NoPath", fcsv, exclusionVars,factorToNumericList,                     
         treeGenerationMinBucket=50,
         treeSummaryMinBucket=20,
         useSubDir=FALSE)  
         
newfeat1=((data$SMIfactor==0) & (data$CAC < 2253) & (data$CAC< 1998) & (data$CAC>=1882)) * 1.0
newfeat2=((data$SMIfactor==1) & (data$SMI < 7837) & (data$SMI >= 7499)) * 1.0
newfeatures=cbind(newfeat1, newfeat2) # create columns for the newly found features
datanew=cbind(data,newfeatures)
thismodel=lm(formula=y ~ .,data=datanew)
expectednew=predict(thismodel,datanew)

requireNamespace("Metrics")
OriginalRMSE = Metrics::rmse(data$y,expected)
NewRMSE = Metrics::rmse(data$y,expectednew)

print(paste("OriginalRMSE = ",OriginalRMSE))
print(paste("NewRMSE = ",NewRMSE))



cleanEx()
nameEx("futuresdata")
### * futuresdata

flush(stderr()); flush(stdout())

### Name: futuresdata
### Title: futuresdata
### Aliases: futuresdata
### Keywords: futuresdata

### ** Examples

data(futuresdata)
head(futuresdata)



cleanEx()
nameEx("generateResidualCutoffCode")
### * generateResidualCutoffCode

flush(stderr()); flush(stdout())

### Name: generateResidualCutoffCode
### Title: generateResidualCutoffCode
### Aliases: generateResidualCutoffCode
### Keywords: saveTree

### ** Examples


require(featurefinder)
data(examples)
generateResidualCutoffCode(data=dat0,"treesAll.txt",treesAll,mainfaclevels, runname,
  treeGenerationMinBucket=treeGenerationMinBucket,
  treeSummaryMinBucket=treeSummaryMinBucket,
  treeSummaryResidualThreshold=treeSummaryResidualThreshold,
  treeSummaryResidualMagnitudeThreshold=treeSummaryResidualMagnitudeThreshold,
  doAllFactors=doAllFactors,
  maxFactorLevels=maxFactorLevels)



cleanEx()
nameEx("generateTrees")
### * generateTrees

flush(stderr()); flush(stdout())

### Name: generateTrees
### Title: generateTrees
### Aliases: generateTrees
### Keywords: generateTrees

### ** Examples


require(featurefinder)
data(examples)
treesThisvar=generateTrees(data=dat0,vars,expr,outputPath=tempdir(),runname,
  treeGenerationMinBucket=treeGenerationMinBucket,
  treeSummaryMinBucket=treeSummaryMinBucket,
  treeSummaryResidualThreshold=treeSummaryResidualThreshold,
  treeSummaryResidualMagnitudeThreshold=treeSummaryResidualMagnitudeThreshold,
  doAllFactors=doAllFactors,
  maxFactorLevels=maxFactorLevels)



cleanEx()
nameEx("getVarAv")
### * getVarAv

flush(stderr()); flush(stdout())

### Name: getVarAv
### Title: getVarAv
### Aliases: getVarAv
### Keywords: saveTree

### ** Examples


require(featurefinder)
data(examples)
av=getVarAv(dat,"expected",pathterms)



cleanEx()
nameEx("i")
### * i

flush(stderr()); flush(stdout())

### Name: i
### Title: i
### Aliases: i
### Keywords: i

### ** Examples

data(i)
head(i)



cleanEx()
nameEx("mainfaclevels")
### * mainfaclevels

flush(stderr()); flush(stdout())

### Name: mainfaclevels
### Title: mainfaclevels
### Aliases: mainfaclevels
### Keywords: mainfaclevels

### ** Examples

data(mainfaclevels)
head(mainfaclevels)



cleanEx()
nameEx("maxFactorLevels")
### * maxFactorLevels

flush(stderr()); flush(stdout())

### Name: maxFactorLevels
### Title: maxFactorLevels
### Aliases: maxFactorLevels
### Keywords: maxFactorLevels

### ** Examples

data(maxFactorLevels)
head(maxFactorLevels)



cleanEx()
nameEx("mpgdata")
### * mpgdata

flush(stderr()); flush(stdout())

### Name: mpgdata
### Title: mpgdata
### Aliases: mpgdata
### Keywords: mpgdata

### ** Examples

data(mpgdata)
head(mpgdata)



cleanEx()
nameEx("names")
### * names

flush(stderr()); flush(stdout())

### Name: names
### Title: names
### Aliases: names
### Keywords: names

### ** Examples

data(names)
head(names)



cleanEx()
nameEx("parseSplits")
### * parseSplits

flush(stderr()); flush(stdout())

### Name: parseSplits
### Title: parseSplits
### Aliases: parseSplits
### Keywords: saveTree

### ** Examples


require(featurefinder)
data(examples)
parseSplits(treesAll[[1]][[2]])



cleanEx()
nameEx("pathterms")
### * pathterms

flush(stderr()); flush(stdout())

### Name: pathterms
### Title: pathterms
### Aliases: pathterms
### Keywords: pathterms

### ** Examples

data(pathterms)
head(pathterms)



cleanEx()
nameEx("printResiduals")
### * printResiduals

flush(stderr()); flush(stdout())

### Name: printResiduals
### Title: printResiduals
### Aliases: printResiduals
### Keywords: saveTree

### ** Examples


require(featurefinder)
data(examples)
printResiduals(fileConn,splitlist[t][[1]],dat, runname, names[t],
  treeSummaryResidualThreshold,treeSummaryMinBucket,
  treeSummaryResidualMagnitudeThreshold)



cleanEx()
nameEx("runname")
### * runname

flush(stderr()); flush(stdout())

### Name: runname
### Title: runname
### Aliases: runname
### Keywords: runname

### ** Examples

data(runname)
head(runname)



cleanEx()
nameEx("saveTree")
### * saveTree

flush(stderr()); flush(stdout())

### Name: saveTree
### Title: saveTree
### Aliases: saveTree
### Keywords: saveTree

### ** Examples


require(featurefinder)
data(examples)
fit1=saveTree(data,vars,expr,i,outputPath=tempdir(),runname,mainfaclevels[1],
     treeGenerationMinBucket)



cleanEx()
nameEx("splitlist")
### * splitlist

flush(stderr()); flush(stdout())

### Name: splitlist
### Title: splitlist
### Aliases: splitlist
### Keywords: splitlist

### ** Examples

data(splitlist)
head(splitlist)



cleanEx()
nameEx("t")
### * t

flush(stderr()); flush(stdout())

### Name: t
### Title: t
### Aliases: t
### Keywords: t

### ** Examples

data(t)
head(t)



cleanEx()
nameEx("tree")
### * tree

flush(stderr()); flush(stdout())

### Name: tree
### Title: tree
### Aliases: tree
### Keywords: tree

### ** Examples

data(tree)
head(tree)



cleanEx()
nameEx("treeGenerationMinBucket")
### * treeGenerationMinBucket

flush(stderr()); flush(stdout())

### Name: treeGenerationMinBucket
### Title: treeGenerationMinBucket
### Aliases: treeGenerationMinBucket
### Keywords: treeGenerationMinBucket

### ** Examples

data(treeGenerationMinBucket)
head(treeGenerationMinBucket)



cleanEx()
nameEx("treeSummaryMinBucket")
### * treeSummaryMinBucket

flush(stderr()); flush(stdout())

### Name: treeSummaryMinBucket
### Title: treeSummaryMinBucket
### Aliases: treeSummaryMinBucket
### Keywords: treeSummaryMinBucket

### ** Examples

data(treeSummaryMinBucket)
head(treeSummaryMinBucket)



cleanEx()
nameEx("treeSummaryResidualMagnitudeThreshold")
### * treeSummaryResidualMagnitudeThreshold

flush(stderr()); flush(stdout())

### Name: treeSummaryResidualMagnitudeThreshold
### Title: treeSummaryResidualMagnitudeThreshold
### Aliases: treeSummaryResidualMagnitudeThreshold
### Keywords: treeSummaryResidualMagnitudeThreshold

### ** Examples

data(treeSummaryResidualMagnitudeThreshold)
head(treeSummaryResidualMagnitudeThreshold)



cleanEx()
nameEx("treeSummaryResidualThreshold")
### * treeSummaryResidualThreshold

flush(stderr()); flush(stdout())

### Name: treeSummaryResidualThreshold
### Title: treeSummaryResidualThreshold
### Aliases: treeSummaryResidualThreshold
### Keywords: treeSummaryResidualThreshold

### ** Examples

data(treeSummaryResidualThreshold)
head(treeSummaryResidualThreshold)



cleanEx()
nameEx("trees")
### * trees

flush(stderr()); flush(stdout())

### Name: trees
### Title: trees
### Aliases: trees
### Keywords: trees

### ** Examples

data(trees)
head(trees)



cleanEx()
nameEx("treesAll")
### * treesAll

flush(stderr()); flush(stdout())

### Name: treesAll
### Title: treesAll
### Aliases: treesAll
### Keywords: treesAll

### ** Examples

data(treesAll)
head(treesAll)



cleanEx()
nameEx("vars")
### * vars

flush(stderr()); flush(stdout())

### Name: vars
### Title: vars
### Aliases: vars
### Keywords: vars

### ** Examples

data(vars)
head(vars)



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
