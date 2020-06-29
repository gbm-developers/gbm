pkgname <- "CMA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('CMA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ElasticNetCMA")
### * ElasticNetCMA

flush(stderr()); flush(stdout())

### Name: ElasticNetCMA
### Title: Classfication and variable selection by the ElasticNet
### Aliases: ElasticNetCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression
golubX <- as.matrix(golub[,-1])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run ElasticNet - penalized logistic regression (no tuning)
result <- ElasticNetCMA(X=golubX, y=golubY, learnind=learnind, norm.fraction = 0.2, alpha=0.5)
show(result)
ftable(result)
plot(result)



cleanEx()
nameEx("GeneSelection")
### * GeneSelection

flush(stderr()); flush(stdout())

### Name: GeneSelection
### Title: General method for variable selection with various methods
### Aliases: GeneSelection
### Keywords: multivariate

### ** Examples

# load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression from first 10 genes
golubX <- as.matrix(golub[,-1])
### Generate five different learningsets
set.seed(111)
five <- GenerateLearningsets(y=golubY, method = "CV", fold = 5, strat = TRUE)
### simple t-test:
selttest <- GeneSelection(golubX, golubY, learningsets = five, method = "t.test")
### show result:
show(selttest)
toplist(selttest, k = 10, iter = 1)
plot(selttest, iter = 1)



cleanEx()
nameEx("GenerateLearningsets")
### * GenerateLearningsets

flush(stderr()); flush(stdout())

### Name: GenerateLearningsets
### Title: Repeated Divisions into learn- and tets sets
### Aliases: GenerateLearningsets
### Keywords: multivariate

### ** Examples

# LOOCV
loo <- GenerateLearningsets(n=40, method="LOOCV")
show(loo)
# five-fold-CV
CV5 <- GenerateLearningsets(n=40, method="CV", fold=5)
show(loo)
# MCCV
mccv <- GenerateLearningsets(n=40, method = "MCCV", niter=3, ntrain=30)
show(mccv)
# Bootstrap
boot <- GenerateLearningsets(n=40, method="bootstrap", niter=3)
# stratified five-fold-CV
set.seed(113)
classlabels <- sample(1:3, size = 50, replace = TRUE, prob = c(0.3, 0.5, 0.2))
CV5strat <- GenerateLearningsets(y = classlabels, method="CV", fold=5, strat = TRUE)
show(CV5strat)



cleanEx()
nameEx("LassoCMA")
### * LassoCMA

flush(stderr()); flush(stdout())

### Name: LassoCMA
### Title: L1 penalized logistic regression
### Aliases: LassoCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression
golubX <- as.matrix(golub[,-1])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run L1 penalized logistic regression (no tuning)
lassoresult <- LassoCMA(X=golubX, y=golubY, learnind=learnind, norm.fraction = 0.2)
show(lassoresult)
ftable(lassoresult)
plot(lassoresult)



cleanEx()
nameEx("Planarplot")
### * Planarplot

flush(stderr()); flush(stdout())

### Name: Planarplot
### Title: Visualize Separability of different classes
### Aliases: Planarplot
### Keywords: multivariate

### ** Examples

### simple linear discrimination for the golub data:
data(golub)
golubY <- golub[,1]
golubX <- as.matrix(golub[,-1])
golubn <- nrow(golubX)
set.seed(111)
learnind <- sample(golubn, size=floor(2/3*golubn))
Planarplot(X=golubX, y=golubY, learnind=learnind, predind=c(2,4),
           classifier=ldaCMA)




cleanEx()
nameEx("classification")
### * classification

flush(stderr()); flush(stdout())

### Name: classification
### Title: General method for classification with various methods
### Aliases: classification
### Keywords: multivariate

### ** Examples

### a simple k-nearest neighbour example
### datasets
## Not run: 
##D plot(x)
##D data(golub)
##D golubY <- golub[,1]
##D golubX <- as.matrix(golub[,-1])
##D ### learningsets
##D set.seed(111)
##D lset <- GenerateLearningsets(y=golubY, method = "CV", fold=5, strat =TRUE)
##D ### 1. GeneSelection
##D selttest <- GeneSelection(golubX, golubY, learningsets = lset, method = "t.test")
##D ### 2. tuning
##D tunek <- tune(golubX, golubY, learningsets = lset, genesel = selttest, nbgene = 20, classifier = knnCMA)
##D ### 3. classification
##D knn1 <- classification(golubX, golubY, learningsets = lset, genesel = selttest,
##D                        tuneres = tunek, nbgene = 20, classifier = knnCMA)
##D ### steps 1.-3. combined into one step:
##D knn2 <- classification(golubX, golubY, learningsets = lset,
##D                        genesellist = list(method  = "t.test"), classifier = knnCMA,
##D                        tuninglist = list(grids = list(k = c(1:8))), nbgene = 20)
##D ### show and analyze results:
##D knnjoin <- join(knn2)
##D show(knn2)
##D eval <- evaluation(knn2, measure = "misclassification")
##D show(eval)
##D summary(eval)
##D boxplot(eval)
## End(Not run)



cleanEx()
nameEx("compBoostCMA")
### * compBoostCMA

flush(stderr()); flush(stdout())

### Name: compBoostCMA
### Title: Componentwise Boosting
### Aliases: compBoostCMA
### Keywords: multivariate

### ** Examples

 ### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression
golubX <- as.matrix(golub[,-1])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run componentwise (logit)-boosting (not tuned)
result <- compBoostCMA(X=golubX, y=golubY, learnind=learnind, mstop = 500)
### show results
show(result)
ftable(result)
plot(result)
### multiclass example:
### load Khan data
data(khan)
### extract class labels
khanY <- khan[,1]
### extract gene expression
khanX <- as.matrix(khan[,-1])
### select learningset
set.seed(111)
learnind <- sample(length(khanY), size=floor(ratio*length(khanY)))
### run componentwise multivariate (logit)-boosting (not tuned)
result <- compBoostCMA(X=khanX, y=khanY, learnind=learnind, mstop = 1000)
### show results
show(result)
ftable(result)
plot(result)



cleanEx()
nameEx("compare")
### * compare

flush(stderr()); flush(stdout())

### Name: compare
### Title: Compare different classifiers
### Aliases: compare
### Keywords: multivariate

### ** Examples

## Not run: 
##D ### compare the performance of several discriminant analysis methods
##D ### for the Khan dataset:
##D data(khan)
##D khanX <- as.matrix(khan[,-1])
##D khanY <- khan[,1]
##D set.seed(27611)
##D fiveCV10iter <- GenerateLearningsets(y=khanY, method = "CV", fold = 5, niter = 2, strat = TRUE)
##D ### candidate methods:  DLDA, LDA, QDA, pls_LDA, sclda
##D class_dlda <- classification(X = khanX, y=khanY, learningsets = fiveCV10iter, classifier = dldaCMA)
##D ### peform GeneSlection for LDA, FDA, QDA (using F-Tests):
##D genesel_da <- GeneSelection(X=khanX, y=khanY, learningsets = fiveCV10iter, method = "f.test")
##D ###
##D class_lda <- classification(X = khanX, y=khanY, learningsets = fiveCV10iter, classifier = ldaCMA, genesel= genesel_da, nbgene = 10)
##D 
##D class_qda <- classification(X = khanX, y=khanY, learningsets = fiveCV10iter, classifier = qdaCMA, genesel = genesel_da, nbgene = 2)
##D 
##D ### We now make a comparison concerning the performance (sev. measures):
##D ### first, collect in a list:
##D dalike <- list(class_dlda, class_lda, class_qda)
##D ### use pre-defined compare function:
##D comparison <- compare(dalike, plot = TRUE, measure = c("misclassification", "brier score", "average probability"))
##D print(comparison)
## End(Not run)



cleanEx()
nameEx("dldaCMA")
### * dldaCMA

flush(stderr()); flush(stdout())

### Name: dldaCMA
### Title: Diagonal Discriminant Analysis
### Aliases: dldaCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression
golubX <- as.matrix(golub[,-1])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run DLDA
dldaresult <- dldaCMA(X=golubX, y=golubY, learnind=learnind)
### show results
show(dldaresult)
ftable(dldaresult)
plot(dldaresult)
### multiclass example:
### load Khan data
data(khan)
### extract class labels
khanY <- khan[,1]
### extract gene expression
khanX <- as.matrix(khan[,-1])
### select learningset
set.seed(111)
learnind <- sample(length(khanY), size=floor(ratio*length(khanY)))
### run LDA
ldaresult <- dldaCMA(X=khanX, y=khanY, learnind=learnind)
### show results
show(dldaresult)
ftable(dldaresult)
plot(dldaresult)



cleanEx()
nameEx("evaluation")
### * evaluation

flush(stderr()); flush(stdout())

### Name: evaluation
### Title: Evaluation of classifiers
### Aliases: evaluation
### Keywords: multivariate

### ** Examples

### simple linear discriminant analysis example using bootstrap datasets:
### datasets:
data(golub)
golubY <- golub[,1]
### extract gene expression from first 10 genes
golubX <- as.matrix(golub[,2:11])
### generate 25 bootstrap datasets
set.seed(333)
bootds <- GenerateLearningsets(y = golubY, method = "bootstrap", ntrain = 30, niter = 10, strat = TRUE)
### run classification()
ldalist <- classification(X=golubX, y=golubY, learningsets = bootds, classifier=ldaCMA)
### Evaluation:
eval_iter <- evaluation(ldalist, scheme = "iter")
eval_obs <- evaluation(ldalist, scheme = "obs")
show(eval_iter)
show(eval_obs)
summary(eval_iter)
summary(eval_obs)
### auc with boxplot
eval_auc <- evaluation(ldalist, scheme = "iter", measure = "auc")
boxplot(eval_auc)
### which observations have often been misclassified ?
obsinfo(eval_obs, threshold = 0.75)



cleanEx()
nameEx("fdaCMA")
### * fdaCMA

flush(stderr()); flush(stdout())

### Name: fdaCMA
### Title: Fisher's Linear Discriminant Analysis
### Aliases: fdaCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression from first 10 genes
golubX <- as.matrix(golub[,2:11])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run FDA
fdaresult <- fdaCMA(X=golubX, y=golubY, learnind=learnind, comp = 1, plot = TRUE)
### show results
show(fdaresult)
ftable(fdaresult)
plot(fdaresult)
### multiclass example:
### load Khan data
data(khan)
### extract class labels
khanY <- khan[,1]
### extract gene expression from first 10 genes
khanX <- as.matrix(khan[,2:11])
### select learningset
set.seed(111)
learnind <- sample(length(khanY), size=floor(ratio*length(khanY)))
### run FDA
fdaresult <- fdaCMA(X=khanX, y=khanY, learnind=learnind, comp = 2, plot = TRUE)
### show results
show(fdaresult)
ftable(fdaresult)
plot(fdaresult)



cleanEx()
nameEx("flexdaCMA")
### * flexdaCMA

flush(stderr()); flush(stdout())

### Name: flexdaCMA
### Title: Flexible Discriminant Analysis
### Aliases: flexdaCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression from first 5 genes
golubX <- as.matrix(golub[,2:6])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run flexible Discriminant Analysis
result <- flexdaCMA(X=golubX, y=golubY, learnind=learnind, comp = 1)
### show results
show(result)
ftable(result)
plot(result)



cleanEx()
nameEx("gbmCMA")
### * gbmCMA

flush(stderr()); flush(stdout())

### Name: gbmCMA
### Title: Tree-based Gradient Boosting
### Aliases: gbmCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression
golubX <- as.matrix(golub[,-1])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run tree-based gradient boosting (no tuning)
gbmresult <- gbmCMA(X=golubX, y=golubY, learnind=learnind, n.trees = 500)
show(gbmresult)
ftable(gbmresult)
plot(gbmresult)



cleanEx()
nameEx("golub")
### * golub

flush(stderr()); flush(stdout())

### Name: golub
### Title: ALL/AML dataset of Golub et al. (1999)
### Aliases: golub
### Keywords: datasets

### ** Examples

data(golub)



cleanEx()
nameEx("khan")
### * khan

flush(stderr()); flush(stdout())

### Name: khan
### Title: Small blue round cell tumor dataset of Khan et al. (2001)
### Aliases: khan
### Keywords: datasets

### ** Examples

data(khan)


cleanEx()
nameEx("knnCMA")
### * knnCMA

flush(stderr()); flush(stdout())

### Name: knnCMA
### Title: Nearest Neighbours
### Aliases: knnCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression from first 10 genes
golubX <- as.matrix(golub[,-1])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run k-nearest neighbours
result <- knnCMA(X=golubX, y=golubY, learnind=learnind, k = 3)
### show results
show(result)
ftable(result)
### multiclass example:
### load Khan data
data(khan)
### extract class labels
khanY <- khan[,1]
### extract gene expression
khanX <- as.matrix(khan[,-1])
### select learningset
set.seed(111)
learnind <- sample(length(khanY), size=floor(ratio*length(khanY)))
### run knn
result <- knnCMA(X=khanX, y=khanY, learnind=learnind, k = 5)
### show results
show(result)
ftable(result)



cleanEx()
nameEx("ldaCMA")
### * ldaCMA

flush(stderr()); flush(stdout())

### Name: ldaCMA
### Title: Linear Discriminant Analysis
### Aliases: ldaCMA
### Keywords: multivariate

### ** Examples

## Not run: 
##D ### load Golub AML/ALL data
##D data(golub)
##D ### extract class labels
##D golubY <- golub[,1]
##D ### extract gene expression from first 10 genes
##D golubX <- as.matrix(golub[,2:11])
##D ### select learningset
##D ratio <- 2/3
##D set.seed(111)
##D learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
##D ### run LDA
##D ldaresult <- ldaCMA(X=golubX, y=golubY, learnind=learnind)
##D ### show results
##D show(ldaresult)
##D ftable(ldaresult)
##D plot(ldaresult)
##D ### multiclass example:
##D ### load Khan data
##D data(khan)
##D ### extract class labels
##D khanY <- khan[,1]
##D ### extract gene expression from first 10 genes
##D khanX <- as.matrix(khan[,2:11])
##D ### select learningset
##D set.seed(111)
##D learnind <- sample(length(khanY), size=floor(ratio*length(khanY)))
##D ### run LDA
##D ldaresult <- ldaCMA(X=khanX, y=khanY, learnind=learnind)
##D ### show results
##D show(ldaresult)
##D ftable(ldaresult)
##D plot(ldaresult)
## End(Not run)


cleanEx()
nameEx("nnetCMA")
### * nnetCMA

flush(stderr()); flush(stdout())

### Name: nnetCMA
### Title: Feed-forward Neural Networks
### Aliases: nnetCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression from first 10 genes
golubX <- as.matrix(golub[,2:11])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run nnet (not tuned)
nnetresult <- nnetCMA(X=golubX, y=golubY, learnind=learnind, size = 3, decay = 0.01)
### show results
show(nnetresult)
ftable(nnetresult)
plot(nnetresult)
### in the space of eigengenes (not tuned)
golubXfull <-  as.matrix(golubX[,-1])
nnetresult <- nnetCMA(X=golubXfull, y=golubY, learnind = learnind, eigengenes = TRUE,
                      size = 3, decay = 0.01)
### show results
show(nnetresult)
ftable(nnetresult)
plot(nnetresult)



cleanEx()
nameEx("pknnCMA")
### * pknnCMA

flush(stderr()); flush(stdout())

### Name: pknnCMA
### Title: Probabilistic Nearest Neighbours
### Aliases: pknnCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression from first 10 genes
golubX <- as.matrix(golub[,-1])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run probabilistic k-nearest neighbours
result <- pknnCMA(X=golubX, y=golubY, learnind=learnind, k = 3)
### show results
show(result)
ftable(result)
plot(result)



cleanEx()
nameEx("plrCMA")
### * plrCMA

flush(stderr()); flush(stdout())

### Name: plrCMA
### Title: L2 penalized logistic regression
### Aliases: plrCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression from first 10 genes
golubX <- as.matrix(golub[,-1])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run penalized logistic regression (no tuning)
plrresult <- plrCMA(X=golubX, y=golubY, learnind=learnind)
### show results
show(plrresult)
ftable(plrresult)
plot(plrresult)
### multiclass example:
### load Khan data
data(khan)
### extract class labels
khanY <- khan[,1]
### extract gene expression from first 10 genes
khanX <- as.matrix(khan[,-1])
### select learningset
set.seed(111)
learnind <- sample(length(khanY), size=floor(ratio*length(khanY)))
### run penalized logistic regression (no tuning)
plrresult <- plrCMA(X=khanX, y=khanY, learnind=learnind)
### show results
show(plrresult)
ftable(plrresult)
plot(plrresult)



cleanEx()
nameEx("pls_ldaCMA")
### * pls_ldaCMA

flush(stderr()); flush(stdout())

### Name: pls_ldaCMA
### Title: Partial Least Squares combined with Linear Discriminant Analysis
### Aliases: pls_ldaCMA
### Keywords: multivariate

### ** Examples

## Not run: 
##D ### load Khan data
##D data(khan)
##D ### extract class labels
##D khanY <- khan[,1]
##D ### extract gene expression
##D khanX <- as.matrix(khan[,-1])
##D ### select learningset
##D set.seed(111)
##D learnind <- sample(length(khanY), size=floor(2/3*length(khanY)))
##D ### run Shrunken Centroids classfier, without tuning
##D plsresult <- pls_ldaCMA(X=khanX, y=khanY, learnind=learnind, comp = 4)
##D ### show results
##D show(plsresult)
##D ftable(plsresult)
##D plot(plsresult)
## End(Not run)


cleanEx()
nameEx("pls_lrCMA")
### * pls_lrCMA

flush(stderr()); flush(stdout())

### Name: pls_lrCMA
### Title: Partial Least Squares followed by logistic regression
### Aliases: pls_lrCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression
golubX <- as.matrix(golub[,-1])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run PLS, combined with logistic regression
result <- pls_lrCMA(X=golubX, y=golubY, learnind=learnind)
### show results
show(result)
ftable(result)
plot(result)



cleanEx()
nameEx("pls_rfCMA")
### * pls_rfCMA

flush(stderr()); flush(stdout())

### Name: pls_rfCMA
### Title: Partial Least Squares followed by random forests
### Aliases: pls_rfCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression
golubX <- as.matrix(golub[,-1])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run PLS, combined with Random Forest
#result <- pls_rfCMA(X=golubX, y=golubY, learnind=learnind)
### show results
#show(result)
#ftable(result)
#plot(result)



cleanEx()
nameEx("pnnCMA")
### * pnnCMA

flush(stderr()); flush(stdout())

### Name: pnnCMA
### Title: Probabilistic Neural Networks
### Aliases: pnnCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression from first 10 genes
golubX <- as.matrix(golub[,2:11])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run PNN
pnnresult <- pnnCMA(X=golubX, y=golubY, learnind=learnind, sigma = 3)
### show results
show(pnnresult)
ftable(pnnresult)
plot(pnnresult)



cleanEx()
nameEx("prediction")
### * prediction

flush(stderr()); flush(stdout())

### Name: prediction
### Title: General method for predicting classes of new observations
### Aliases: prediction
### Keywords: multivariate

### ** Examples

### a simple k-nearest neighbour example
### datasets
## Not run: 
##D plot(x)
##D data(golub)
##D golubY <- golub[,1]
##D golubX <- as.matrix(golub[,-1])
##D ###Splitting data into training and test set
##D X.tr<-golubX[1:30]
##D X.new<-golubX[31:39]
##D y.tr<-golubY[1:30]
##D ### 1. GeneSelection
##D selttest <- GeneSelection(X=X.tr, y=y.tr, method = "t.test")
##D ### 2. tuning
##D tunek <- tune(X.tr, y.tr, genesel = selttest, nbgene = 20, classifier = knnCMA)
##D ### 3. classification
##D pred <- prediction(X.tr=X.tr,y.tr=y.tr,X.new=X.new, genesel = selttest,
##D                        tuneres = tunek, nbgene = 20, classifier = knnCMA)
##D ### show and analyze results:
##D show(pred)
##D 
## End(Not run)



cleanEx()
nameEx("qdaCMA")
### * qdaCMA

flush(stderr()); flush(stdout())

### Name: qdaCMA
### Title: Quadratic Discriminant Analysis
### Aliases: qdaCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression from first 3 genes
golubX <- as.matrix(golub[,2:4])
### select learningset
ratio <- 2/3
set.seed(112)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run QDA
qdaresult <- qdaCMA(X=golubX, y=golubY, learnind=learnind)
### show results
show(qdaresult)
ftable(qdaresult)
plot(qdaresult)
### multiclass example:
### load Khan data
data(khan)
### extract class labels
khanY <- khan[,1]
### extract gene expression from first 4 genes
khanX <- as.matrix(khan[,2:5])
### select learningset
set.seed(111)
learnind <- sample(length(khanY), size=floor(ratio*length(khanY)))
### run QDA
qdaresult <- qdaCMA(X=khanX, y=khanY, learnind=learnind)
### show results
show(qdaresult)
ftable(qdaresult)
plot(qdaresult)



cleanEx()
nameEx("rfCMA")
### * rfCMA

flush(stderr()); flush(stdout())

### Name: rfCMA
### Title: Classification based on Random Forests
### Aliases: rfCMA
### Keywords: multivariate

### ** Examples

 ### load Khan data
data(khan)
### extract class labels
khanY <- khan[,1]
### extract gene expression
khanX <- as.matrix(khan[,-1])
### select learningset
set.seed(111)
learnind <- sample(length(khanY), size=floor(2/3*length(khanY)))
### run random Forest
#rfresult <- rfCMA(X=khanX, y=khanY, learnind=learnind, varimp = FALSE)
### show results
#show(rfresult)
#ftable(rfresult)
#plot(rfresult)


cleanEx()
nameEx("scdaCMA")
### * scdaCMA

flush(stderr()); flush(stdout())

### Name: scdaCMA
### Title: Shrunken Centroids Discriminant Analysis
### Aliases: scdaCMA
### Keywords: multivariate

### ** Examples

### load Khan data
data(khan)
### extract class labels
khanY <- khan[,1]
### extract gene expression
khanX <- as.matrix(khan[,-1])
### select learningset
set.seed(111)
learnind <- sample(length(khanY), size=floor(2/3*length(khanY)))
### run Shrunken Centroids classfier, without tuning
scdaresult <- scdaCMA(X=khanX, y=khanY, learnind=learnind)
### show results
show(scdaresult)
ftable(scdaresult)
plot(scdaresult)


cleanEx()
nameEx("shrinkldaCMA")
### * shrinkldaCMA

flush(stderr()); flush(stdout())

### Name: shrinkldaCMA
### Title: Shrinkage linear discriminant analysis
### Aliases: shrinkldaCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression
golubX <- as.matrix(golub[,-1])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run  shrinkage-LDA
result <- shrinkldaCMA(X=golubX, y=golubY, learnind=learnind)
### show results
show(result)
ftable(result)
plot(result)


cleanEx()
nameEx("svmCMA")
### * svmCMA

flush(stderr()); flush(stdout())

### Name: svmCMA
### Title: Support Vector Machine
### Aliases: svmCMA
### Keywords: multivariate

### ** Examples

### load Golub AML/ALL data
data(golub)
### extract class labels
golubY <- golub[,1]
### extract gene expression
golubX <- as.matrix(golub[,-1])
### select learningset
ratio <- 2/3
set.seed(111)
learnind <- sample(length(golubY), size=floor(ratio*length(golubY)))
### run _untuned_linear SVM
svmresult <- svmCMA(X=golubX, y=golubY, learnind=learnind,probability=TRUE)
### show results
show(svmresult)
ftable(svmresult)
plot(svmresult)


cleanEx()
nameEx("tune")
### * tune

flush(stderr()); flush(stdout())

### Name: tune
### Title: Hyperparameter tuning for classifiers
### Aliases: tune
### Keywords: multivariate

### ** Examples

## Not run: 
##D ### simple example for a one-dimensional grid, using compBoostCMA.
##D ### dataset
##D data(golub)
##D golubY <- golub[,1]
##D golubX <- as.matrix(golub[,-1])
##D ### learningsets
##D set.seed(111)
##D lset <- GenerateLearningsets(y=golubY, method = "CV", fold=5, strat =TRUE)
##D ### tuning after gene selection with the t.test
##D tuneres <- tune(X = golubX, y = golubY, learningsets = lset,
##D               genesellist = list(method = "t.test"),
##D               classifier=compBoostCMA, nbgene = 100,
##D               grids = list(mstop = c(50, 100, 250, 500, 1000)))
##D ### inspect results
##D show(tuneres)
##D best(tuneres)
##D plot(tuneres, iter = 3)
## End(Not run)



cleanEx()
nameEx("weighted_mcr")
### * weighted_mcr

flush(stderr()); flush(stdout())

### Name: weighted.mcr
### Title: Tuning / Selection bias correction
### Aliases: weighted.mcr
### Keywords: tuning bias, selection bias, corrected misclassification rate

### ** Examples

#inputs
classifiers<-rep('knnCMA',7)
nbgenes<-rep(50,7)
parameters<-c('k=1','k=3','k=5','k=7','k=9','k=11','k=13')
portion<-0.8
niter<-100
data(golub)
X<-as.matrix(golub[,-1])         
y<-golub[,1]
sel.method<-'t.test'
#function call
wmcr<-weighted.mcr(classifiers=classifiers,parameters=parameters,nbgenes=nbgenes,sel.method=sel.method,X=X,y=y,portion=portion,niter=niter)



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
