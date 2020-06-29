pkgname <- "MiDA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('MiDA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("MiBiClassGBODT")
### * MiBiClassGBODT

flush(stderr()); flush(stdout())

### Name: MiBiClassGBODT
### Title: Binary classification using gradient boosting over desicion
###   trees
### Aliases: MiBiClassGBODT

### ** Examples


#get gene expression and specimen data
data("IMexpression");data("IMspecimen")
#sample expression matrix and specimen data for binary classification,
#only "NORM" and "EBV" specimens are left
SampleMatrix<-MiDataSample(IMexpression, IMspecimen$diagnosis,"norm", "ebv")
SampleSpecimen<-MiSpecimenSample(IMspecimen$diagnosis, "norm", "ebv")
#Fitting, low tuning for faster running
BoostRes<-MiBiClassGBODT(SampleMatrix, SampleSpecimen, n.crossval = 3,
                       ntrees = 10, shrinkage = 1, intdepth = 2)
BoostRes[[1]] # QC values for n.crossval = 3 models and its summary
length(BoostRes[[2]]) # n.crossval = 3 data frames of probes feature importance for classification
head(BoostRes[[2]][[1]])




cleanEx()
nameEx("MiDataSample")
### * MiDataSample

flush(stderr()); flush(stdout())

### Name: MiDataSample
### Title: Select matrix columns based on values of attendant vector
### Aliases: MiDataSample

### ** Examples


#get gene expression and specimen data
data("IMexpression");data("IMspecimen")
dim(IMexpression) # 100 columns (genes/transcripts) - 89 specimens
colnames(IMexpression)[1:10] # look at first 10 columns of matrix - specimens IDs
IMspecimen[1:10,] # specimens IDs and group factor - diagnoses in attendant vector
# note that specimens in matrix columns are in the same order as specimens in description data
# select specimens with only EBV and NORM diagnoses (and sample the description data as well)
SampleMatrix<-MiDataSample(IMexpression, IMspecimen$diagnosis, "ebv", "norm")
SampleSamples<-MiSpecimenSample(IMspecimen$diagnosis, "ebv", "norm")
dim(SampleMatrix)# only 68 specimens with EBV and NORM diagnoses left
colnames(SampleMatrix)[1:10]
SampleSamples[1:10] # corresponding diagnoses




cleanEx()
nameEx("MiFracData")
### * MiFracData

flush(stderr()); flush(stdout())

### Name: MiFracData
### Title: Subset an expression matrix based on probe's feature importance
### Aliases: MiFracData

### ** Examples

# get gene expression and specimen data
data("IMexpression");data("IMspecimen")
#sample expression matrix and specimen data for binary classification,
#only "NORM" and "EBV" specimens are left
SampleMatrix<-MiDataSample(IMexpression, IMspecimen$diagnosis,"norm", "ebv")
dim(SampleMatrix) # 100 probes
SampleSpecimen<-MiSpecimenSample(IMspecimen$diagnosis, "norm", "ebv")
#Fitting, low tuning for faster running
ClassRes<-MiBiClassGBODT(SampleMatrix, SampleSpecimen, n.crossval = 3,
                         ntrees = 10, shrinkage = 1, intdepth = 2)
# List of influence data frames for all 3 models build using cross-validation
# is the 2nd element of BiClassGBODT results
# take 10 most important probes from each model
Sample2Matrix<-MiFracData(SampleMatrix, importance.list = ClassRes[[2]], 10)
dim(Sample2Matrix) # less than 100 probes left




cleanEx()
nameEx("MiInflCount")
### * MiInflCount

flush(stderr()); flush(stdout())

### Name: MiInflCount
### Title: Mean microarray probes' feature importance from binary
###   classification
### Aliases: MiInflCount

### ** Examples


# get gene expression and specimen data
data("IMexpression");data("IMspecimen")
# sample expression matrix and specimen data for binary classification,
# only "NORM" and "EBV" specimens are left
SampleMatrix<-MiDataSample(IMexpression, IMspecimen$diagnosis,"norm", "ebv")
SampleSpecimen<-MiSpecimenSample(IMspecimen$diagnosis, "norm", "ebv")
#Fitting, low tuning for faster running
ClassRes<-MiBiClassGBODT(SampleMatrix, SampleSpecimen, n.crossval = 3,
                        ntrees = 10, shrinkage = 1, intdepth = 2)
# List of influence data frames for all 3 models build using cross-validation
# is the 2nd element of BiClassGBODT results
Importances<-MiInflCount(ClassRes[[2]])
Importances[[1]][1:10,] # mean and sd. 0s are for low feature importance
Importances[[2]][1:10,] # original values for n.crossval = 3 models




cleanEx()
nameEx("MiIntDepthAjust")
### * MiIntDepthAjust

flush(stderr()); flush(stdout())

### Name: MiIntDepthAjust
### Title: Ajust maximum depth parameter for fitting generalized boosted
###   regression models
### Aliases: MiIntDepthAjust

### ** Examples

#get gene expression and specimen data
data("IMexpression");data("IMspecimen")
#sample expression matrix and specimen data for binary classification,
#only "NORM" and "EBV" specimens are left
SampleMatrix<-MiDataSample(IMexpression, IMspecimen$diagnosis,"norm", "ebv")
SampleSpecimen<-MiSpecimenSample(IMspecimen$diagnosis, "norm", "ebv")
#Fitting, low tuning for faster running. Test intdepth
set.seed(1)
ClassRes<-MiIntDepthAjust(SampleMatrix, SampleSpecimen, test.frac = 5, times=3,
                          ntrees = 10, shrinkage = 1, intdepth =  c(1,2))
ClassRes[[1]] # train accuracy
ClassRes[[2]] # test accuracy




cleanEx()
nameEx("MiNTreesAjust")
### * MiNTreesAjust

flush(stderr()); flush(stdout())

### Name: MiNTreesAjust
### Title: Ajust number of trees parameter for fitting generalized boosted
###   regression models
### Aliases: MiNTreesAjust

### ** Examples

#get gene expression and specimen data
data("IMexpression");data("IMspecimen")
#sample expression matrix and specimen data for binary classification,
#only "NORM" and "EBV" specimens are left
SampleMatrix<-MiDataSample(IMexpression, IMspecimen$diagnosis,"norm", "ebv")
SampleSpecimen<-MiSpecimenSample(IMspecimen$diagnosis, "norm", "ebv")
#Fitting, low tuning for faster running. Test ntrees
set.seed(1)
ClassRes<-MiNTreesAjust(SampleMatrix, SampleSpecimen, test.frac = 5, times = 3,
                       ntrees = c(10, 20), shrinkage = 1, intdepth = 2)
ClassRes[[1]] # train accuracy
ClassRes[[2]] # test accuracy





cleanEx()
nameEx("MiNorm")
### * MiNorm

flush(stderr()); flush(stdout())

### Name: MiNorm
### Title: Microarray data normalization
### Aliases: MiNorm

### ** Examples

data("IMexpression")
# Loess normalization
LoMatrix<-MiNorm(IMexpression, method="Loess")
par(mfrow=c(1,2))
boxplot(log2(IMexpression),main="Before normalization")
boxplot(log2(LoMatrix),main="Loess normalization")
par(mfrow=c(1,1))





graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("MiSelectSignif")
### * MiSelectSignif

flush(stderr()); flush(stdout())

### Name: MiSelectSignif
### Title: Select biological markers with high fold change and
###   classification importance
### Aliases: MiSelectSignif

### ** Examples

probes<-paste("probe", 1:50, sep="") #probes
mean1<-rnorm(50, mean=0, sd=1) #means
mean2<-rnorm(50, mean=5, sd=1)
infl<-c(1:50) # influence
stat.val<-rep(c(0.05, 0.04), c(20, 30))
Result<-MiSelectSignif(probes, mean1, mean2, FC.method="absolute", infl, stat.val,
                      tresh.FC=0.75, tresh.infl=0.75, tresh.stat=0.05)
Result[1:5,]




cleanEx()
nameEx("MiShrinkAjust")
### * MiShrinkAjust

flush(stderr()); flush(stdout())

### Name: MiShrinkAjust
### Title: Ajust learning rate parameter for fitting generalized boosted
###   regression modelsfor fitting generalized boosted regression models
### Aliases: MiShrinkAjust

### ** Examples

data("IMexpression");data("IMspecimen")
#sample expression matrix and specimen data for binary classification,
#only "NORM" and "EBV" specimens are left
SampleMatrix<-MiDataSample(IMexpression, IMspecimen$diagnosis,"norm", "ebv")
SampleSpecimen<-MiSpecimenSample(IMspecimen$diagnosis, "norm", "ebv")
#Fitting, low tuning for faster running. Test shrinkage
set.seed(1)
ClassRes<-MiShrinkAjust(SampleMatrix, SampleSpecimen, test.frac = 5, times = 3,
                        ntrees = 10, shrinkage = c(0.1, 1), intdepth = 2)
ClassRes[[1]] # train accuracy
ClassRes[[2]] # test accuracy




cleanEx()
nameEx("MiSpecimenSample")
### * MiSpecimenSample

flush(stderr()); flush(stdout())

### Name: MiSpecimenSample
### Title: Select values from factor vector
### Aliases: MiSpecimenSample

### ** Examples

#get gene expression and specimen data
data("IMexpression");data("IMspecimen")
dim(IMexpression) # 100 columns (genes/transcripts) - 89 specimens
colnames(IMexpression)[1:10] # look at first 10 columns of matrix - specimens IDs
IMspecimen[1:10,] # specimens IDs and group factor - diagnoses in attendant vector
# note that specimens in matrix columns are in the same order as specimens in description data
# select specimens with only EBV and NORM diagnoses (and sample the description data as well)
SampleMatrix<-MiDataSample(IMexpression, IMspecimen$diagnosis, "ebv", "norm")
SampleSamples<-MiSpecimenSample(IMspecimen$diagnosis, "ebv", "norm")
dim(SampleMatrix)# only 68 specimens with EBV and NORM diagnoses left
colnames(SampleMatrix)[1:10]
SampleSamples[1:10] # corresponding diagnoses




cleanEx()
nameEx("MiStatCount")
### * MiStatCount

flush(stderr()); flush(stdout())

### Name: MiStatCount
### Title: FDR for microarray gene expression data
### Aliases: MiStatCount

### ** Examples

data("IMexpression"); data("IMspecimen") # load data and specimen information
#sampling data and specimen information
ExpData<-MiDataSample(IMexpression, IMspecimen$diagnosis,"ebv", "norm")
Specimens<-MiSpecimenSample(IMspecimen$diagnosis, "ebv", "norm")
#Counting statistics
StatRes<-MiStatCount(ExpData, Specimens)
head(StatRes)




cleanEx()
nameEx("MiSummarize")
### * MiSummarize

flush(stderr()); flush(stdout())

### Name: MiSummarize
### Title: Microarray data summarization
### Aliases: MiSummarize

### ** Examples

data("IMexpression") # load data
# See 5 zonds to AGTR2.NM_000686
IMexpression [1:10, 1:5]
SumMatrix<-MiSummarize(IMexpression, sep=".")
# now there is median expression for AGTR2.NM_000686
SumMatrix[ 1:10, 1:5]




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
