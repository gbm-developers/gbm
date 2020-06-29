pkgname <- "fscaret"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('fscaret')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("dataPreprocess")
### * dataPreprocess

flush(stderr()); flush(stdout())

### Name: dataPreprocess
### Title: dataPreprocess
### Aliases: dataPreprocess
### Keywords: univar robust

### ** Examples



library(fscaret)

# Create data sets and labels data frame
trainMatrix <- matrix(rnorm(150*120,mean=10,sd=1), 150, 120)

# Adding some near-zero variance attributes

temp1 <- matrix(runif(150,0.0001,0.0005), 150, 12)

# Adding some highly correlated attributes

sampleColIndex <- sample(ncol(trainMatrix), size=10)

temp2 <- matrix(trainMatrix[,sampleColIndex]*2, 150, 10)

# Output variable

output <- matrix(rnorm(150,mean=10,sd=1), 150, 1)

trainMatrix <- cbind(trainMatrix,temp1,temp2, output)

colnames(trainMatrix) <- paste("X",c(1:ncol(trainMatrix)),sep="")

# Subset test data set

testMatrix <- trainMatrix[sample(round(0.1*nrow(trainMatrix))),]

labelsDF <- data.frame("Labels"=paste("X",c(1:(ncol(trainMatrix)-1)),sep=""))

lk_col <- ncol(trainMatrix)
lk_row <- nrow(trainMatrix)

with.labels = TRUE

testRes <- dataPreprocess(trainMatrix, testMatrix,
			  labelsDF, lk_col, lk_row, with.labels)
			  
summary(testRes)

# Selected attributes after data set preprocessing
testRes$labelsDF

# Training and testing data sets after preprocessing
testRes$trainMatryca
testRes$testMatryca




cleanEx()
nameEx("fscaret")
### * fscaret

flush(stderr()); flush(stdout())

### Name: fscaret
### Title: feature selection caret
### Aliases: fscaret
### Keywords: methods iteration optimize array

### ** Examples


if((Sys.info()['sysname'])!="SunOS"){

library(fscaret)

# Load data sets
data(dataset.train)
data(dataset.test)

requiredPackages <- c("R.utils", "gsubfn", "ipred", "caret", "parallel", "MASS")

if(.Platform$OS.type=="windows"){

myFirstRES <- fscaret(dataset.train, dataset.test, installReqPckg=FALSE,
                  preprocessData=FALSE, with.labels=TRUE, classPred=FALSE,
                  regPred=TRUE, skel_outfile=NULL,
                  impCalcMet="RMSE&MSE", myTimeLimit=4,
                  Used.funcRegPred=c("lm"), Used.funcClassPred=NULL,
                  no.cores=1, method="boot", returnResamp="all",
                  supress.output=TRUE,saveModel=FALSE)

} else {

myCores <- 2

myFirstRES <- fscaret(dataset.train, dataset.test, installReqPckg=FALSE,
                  preprocessData=FALSE, with.labels=TRUE, classPred=FALSE,
                  regPred=TRUE, skel_outfile=NULL,
                  impCalcMet="RMSE&MSE", myTimeLimit=4,
                  Used.funcRegPred=c("lm","ppr"), Used.funcClassPred=NULL,
                  no.cores=myCores, method="boot", returnResamp="all",
                  supress.output=TRUE,saveModel=FALSE)

}



# Results
myFirstRES

}




cleanEx()
nameEx("impCalc")
### * impCalc

flush(stderr()); flush(stdout())

### Name: impCalc
### Title: impCalc
### Aliases: impCalc
### Keywords: design models

### ** Examples


## Not run: 
##D # 
##D # Hashed to comply with new CRAN check
##D # 
##D library(fscaret)
##D 
##D # Load dataset
##D data(dataset.train)
##D data(dataset.test)
##D 
##D # Make objects
##D trainDF <- dataset.train
##D testDF <- dataset.test
##D model <- c("lm","Cubist")
##D fitControl <- trainControl(method = "boot", returnResamp = "all") 
##D myTimeLimit <- 5
##D no.cores <- 2
##D supress.output <- TRUE
##D skel_outfile <- paste("_default_",sep="")
##D mySystem <- .Platform$OS.type
##D with.labels <- TRUE
##D redPred <- TRUE
##D classPred <- FALSE
##D saveModel <- FALSE
##D lvlScale <- FALSE
##D 
##D if(mySystem=="windows"){
##D no.cores <- 1
##D }
##D 
##D # Scan dimensions of trainDF [lk_row x lk_col]
##D lk_col = ncol(trainDF)
##D lk_row = nrow(trainDF)
##D 
##D # Read labels of trainDF
##D labelsFrame <- as.data.frame(colnames(trainDF))
##D labelsFrame <-cbind(c(1:ncol(trainDF)),labelsFrame)
##D # Create a train data set matrix
##D trainMatryca_nr <- matrix(data=NA,nrow=lk_row,ncol=lk_col)
##D 
##D row=0
##D col=0
##D 
##D for(col in 1:(lk_col)) {
##D    for(row in 1:(lk_row)) {
##D      trainMatryca_nr[row,col] <- (as.numeric(trainDF[row,col]))
##D     }
##D }
##D 
##D # Pointing standard data set train
##D xTrain <- data.frame(trainMatryca_nr[,-lk_col])
##D yTrain <- as.vector(trainMatryca_nr[,lk_col])
##D 
##D 
##D #--------Scan dimensions of trainDataFrame1 [lk_row x lk_col]
##D lk_col_test = ncol(testDF)
##D lk_row_test = nrow(testDF)
##D 
##D testMatryca_nr <- matrix(data=NA,nrow=lk_row_test,ncol=lk_col_test)
##D 
##D row=0
##D col=0
##D 
##D for(col in 1:(lk_col_test)) {
##D    for(row in 1:(lk_row_test)) {
##D      testMatryca_nr[row,col] <- (as.numeric(testDF[row,col]))
##D     }
##D }
##D 
##D # Pointing standard data set test
##D xTest <- data.frame(testMatryca_nr[,-lk_col])
##D yTest <- as.vector(testMatryca_nr[,lk_col])
##D 
##D 
##D # Calling low-level function to create models to calculate on
##D myVarImp <- regVarImp(model, xTrain, yTrain, xTest,
##D 	    fitControl, myTimeLimit, no.cores, lk_col,
##D 	    supress.output, mySystem)
##D 
##D 
##D myImpCalc <- impCalc(skel_outfile, xTest, yTest,
##D               lk_col,labelsFrame,with.labels,redPred,classPred,saveModel,lvlScale)
##D 
## End(Not run)




cleanEx()
nameEx("imputeMean")
### * imputeMean

flush(stderr()); flush(stdout())

### Name: imputeMean
### Title: imputeMean
### Aliases: impute.mean
### Keywords: math logic

### ** Examples


library(fscaret)

# Make sample matrix
testData <- matrix(data=rep(1:5),ncol=10,nrow=15)

# Replace random values with NA's
n <- 15
replace <- TRUE
set.seed(1)

rand.sample <- sample(length(testData), n, replace=replace)
testData[rand.sample] <- NA 

# Print out input matrix
testData

# Record cols with missing values
missing.colsTestMatrix <- which(colSums(is.na(testData))>0)

for(i in 1:length(missing.colsTestMatrix)){

rowToReplace <- missing.colsTestMatrix[i]
testData[,rowToReplace] <- impute.mean(testData[,rowToReplace])

}

# Print out matrix with replaced NA's by column mean 
testData




cleanEx()
nameEx("requiredPackages")
### * requiredPackages

flush(stderr()); flush(stdout())

### Name: requiredPackages
### Title: requiredPackages
### Aliases: requiredPackages
### Keywords: datasets

### ** Examples

data(requiredPackages)



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
