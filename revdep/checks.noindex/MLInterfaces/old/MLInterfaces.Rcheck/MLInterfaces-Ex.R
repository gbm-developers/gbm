pkgname <- "MLInterfaces"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('MLInterfaces')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("MLearn-new")
### * MLearn-new

flush(stderr()); flush(stdout())

### Name: MLearn
### Title: revised MLearn interface for machine learning
### Aliases: MLearn_new MLearn baggingI dlda glmI.logistic knnI knn.cvI
###   ksvmI ldaI lvqI naiveBayesI nnetI qdaI RABI randomForestI rpartI svmI
###   svm2 ksvm2 plsda2 plsdaI dlda2 dldaI sldaI blackboostI knn2 knn.cv2
###   ldaI.predParms lvq rab adaI
###   MLearn,formula,ExpressionSet,character,numeric-method
###   MLearn,formula,ExpressionSet,learnerSchema,numeric-method
###   MLearn,formula,data.frame,learnerSchema,numeric-method
###   MLearn,formula,data.frame,learnerSchema,xvalSpec-method
###   MLearn,formula,ExpressionSet,learnerSchema,xvalSpec-method
###   MLearn,formula,data.frame,clusteringSchema,ANY-method plotXvalRDA
###   rdacvI rdaI BgbmI gbm2 rdaML rdacvML hclustI kmeansI pamI
###   makeLearnerSchema standardMLIConverter
### Keywords: models

### ** Examples

library("MASS")
data(crabs)
set.seed(1234)
kp = sample(1:200, size=120)
rf1 = MLearn(sp~CW+RW, data=crabs, randomForestI, kp, ntree=600 )
rf1
nn1 = MLearn(sp~CW+RW, data=crabs, nnetI, kp, size=3, decay=.01,
    trace=FALSE )
nn1
RObject(nn1)
knn1 = MLearn(sp~CW+RW, data=crabs, knnI(k=3,l=2), kp)
knn1
names(RObject(knn1))
dlda1 = MLearn(sp~CW+RW, data=crabs, dldaI, kp )
dlda1
names(RObject(dlda1))
lda1 = MLearn(sp~CW+RW, data=crabs, ldaI, kp )
lda1
names(RObject(lda1))
slda1 = MLearn(sp~CW+RW, data=crabs, sldaI, kp )
slda1
names(RObject(slda1))
svm1 = MLearn(sp~CW+RW, data=crabs, svmI, kp )
svm1
names(RObject(svm1))
ldapp1 = MLearn(sp~CW+RW, data=crabs, ldaI.predParms(method="debiased"), kp )
ldapp1
names(RObject(ldapp1))
qda1 = MLearn(sp~CW+RW, data=crabs, qdaI, kp )
qda1
names(RObject(qda1))
logi = MLearn(sp~CW+RW, data=crabs, glmI.logistic(threshold=0.5), kp, family=binomial ) # need family
logi
names(RObject(logi))
rp2 = MLearn(sp~CW+RW, data=crabs, rpartI, kp)
rp2
## recode data for RAB
#nsp = ifelse(crabs$sp=="O", -1, 1)
#nsp = factor(nsp)
#ncrabs = cbind(nsp,crabs)
#rab1 = MLearn(nsp~CW+RW, data=ncrabs, RABI, kp, maxiter=10)
#rab1
#
# new approach to adaboost
#
ada1 = MLearn(sp ~ CW+RW, data = crabs, .method = adaI, 
    trainInd = kp, type = "discrete", iter = 200)
ada1
confuMat(ada1)
#
lvq.1 = MLearn(sp~CW+RW, data=crabs, lvqI, kp )
lvq.1
nb.1 = MLearn(sp~CW+RW, data=crabs, naiveBayesI, kp )
confuMat(nb.1)
bb.1 = MLearn(sp~CW+RW, data=crabs, baggingI, kp )
confuMat(bb.1)
#
# new mboost interface -- you MUST supply family for nonGaussian response
#
require(party)  # trafo ... killing cmd check
blb.1 = MLearn(sp~CW+RW+FL, data=crabs, blackboostI, kp, family=mboost::Binomial() )
confuMat(blb.1)
#
# ExpressionSet illustration
# 
data(sample.ExpressionSet)
#  needed to increase training set size to avoid a new randomForest condition
# on empty class
set.seed(1234)
X = MLearn(type~., sample.ExpressionSet[100:250,], randomForestI, 1:19, importance=TRUE )
library(randomForest)
library(hgu95av2.db)
opar = par(no.readonly=TRUE)
par(las=2)
plot(getVarImp(X), n=10, plat="hgu95av2", toktype="SYMBOL")
par(opar)
#
# demonstrate cross validation
#
nn1cv = MLearn(sp~CW+RW, data=crabs[c(1:20,101:120),], 
   nnetI, xvalSpec("LOO"), size=3, decay=.01, trace=FALSE )
confuMat(nn1cv)
nn2cv = MLearn(sp~CW+RW, data=crabs[c(1:20,101:120),], nnetI, 
   xvalSpec("LOG",5, balKfold.xvspec(5)), size=3, decay=.01,
   trace=FALSE )
confuMat(nn2cv)
nn3cv = MLearn(sp~CW+RW+CL+BD+FL, data=crabs[c(1:20,101:120),], nnetI, 
   xvalSpec("LOG",5, balKfold.xvspec(5), fsFun=fs.absT(2)), size=3, decay=.01,
   trace=FALSE )
confuMat(nn3cv)
nn4cv = MLearn(sp~.-index-sex, data=crabs[c(1:20,101:120),], nnetI, 
   xvalSpec("LOG",5, balKfold.xvspec(5), fsFun=fs.absT(2)), size=3, decay=.01,
   trace=FALSE )
confuMat(nn4cv)
#
# try with expression data
#
library(golubEsets)
data(Golub_Train)
litg = Golub_Train[ 100:150, ]
g1 = MLearn(ALL.AML~. , litg, nnetI, 
   xvalSpec("LOG",5, balKfold.xvspec(5), 
   fsFun=fs.probT(.75)), size=3, decay=.01, trace=FALSE )
confuMat(g1)
#
# illustrate rda.cv interface from package rda (requiring local bridge)
#
library(ALL)
data(ALL)
#
# restrict to BCR/ABL or NEG
#
bio <- which( ALL$mol.biol %in% c("BCR/ABL", "NEG"))
#
# restrict to B-cell
#
isb <- grep("^B", as.character(ALL$BT))
kp <- intersect(bio,isb)
all2 <- ALL[,kp]
mads = apply(exprs(all2),1,mad)
kp = which(mads>1)  # get around 250 genes
vall2 = all2[kp, ]
vall2$mol.biol = factor(vall2$mol.biol) # drop unused levels

if (requireNamespace("rda", quietly=TRUE)) {
 library("rda")
 r1 = MLearn(mol.biol~., vall2, MLInterfaces:::rdacvI, 1:40)
 confuMat(r1)
 RObject(r1)
 MLInterfaces:::plotXvalRDA(r1)  # special interface to plots of parameter space
}

# illustrate clustering support

cl1 = MLearn(~CW+RW+CL+FL+BD, data=crabs, hclustI(distFun=dist, cutParm=list(k=4)))
plot(cl1)

cl1a = MLearn(~CW+RW+CL+FL+BD, data=crabs, hclustI(distFun=dist, cutParm=list(k=4)), 
   method="complete")
plot(cl1a)

cl2 = MLearn(~CW+RW+CL+FL+BD, data=crabs, kmeansI, centers=5, algorithm="Hartigan-Wong")
plot(cl2, crabs[,-c(1:3)])

c3 = MLearn(~CL+CW+RW, crabs, pamI(dist), k=5)
c3
plot(c3, data=crabs[,c("CL", "CW", "RW")])


#  new interfaces to PLS thanks to Laurent Gatto

set.seed(1234)
kp = sample(1:200, size=120)

#plsda.1 = MLearn(sp~CW+RW, data=crabs, plsdaI, kp, probMethod="Bayes")
#plsda.1
#confuMat(plsda.1)
#confuMat(plsda.1,t=.65) ## requires at least 0.65 post error prob to assign species
#
#plsda.2 = MLearn(type~., data=sample.ExpressionSet[100:250,], plsdaI, 1:16)
#plsda.2
#confuMat(plsda.2)
#confuMat(plsda.2,t=.65) ## requires at least 0.65 post error prob to assign outcome

## examples for predict
#clout <- MLearn(type~., sample.ExpressionSet[100:250,], svmI , 1:16)
#predict(clout, sample.ExpressionSet[100:250,17:26])




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("RAB")
### * RAB

flush(stderr()); flush(stdout())

### Name: RAB
### Title: real adaboost (Friedman et al)
### Aliases: RAB RAB4es DAB Predict tonp mkfmla Predict,raboostCont-method
###   Predict,daboostCont-method
### Keywords: models

### ** Examples

library(MASS)
library(rpart)
data(Pima.tr)
data(Pima.te)
Pima.all = rbind(Pima.tr, Pima.te)
tonp = ifelse(Pima.all$type == "Yes", 1, -1)
tonp = factor(tonp)
Pima.all = data.frame(Pima.all[,1:7], mtype=tonp)
fit1 = RAB(mtype~ped+glu+npreg+bmi+age, data=Pima.all[1:200,], maxiter=10, maxdepth=5)
pfit1 = Predict(fit1, newdata=Pima.tr)
table(Pima.tr$type, pfit1)



cleanEx()
nameEx("balKfold.xvspec")
### * balKfold.xvspec

flush(stderr()); flush(stdout())

### Name: balKfold.xvspec
### Title: generate a partition function for cross-validation, where the
###   partitions are approximately balanced with respect to the
###   distribution of a response variable
### Aliases: balKfold.xvspec
### Keywords: models manip

### ** Examples

## The function is currently defined as
function (K) 
function(data, clab, iternum) {
    clabs <- data[[clab]]
    narr <- nrow(data)
    cnames <- unique(clabs)
    ilist <- list()
    for (i in 1:length(cnames)) ilist[[cnames[i]]] <- which(clabs == 
        cnames[i])
    clens <- lapply(ilist, length)
    nrep <- lapply(clens, function(x) ceiling(x/K))
    grpinds <- list()
    for (i in 1:length(nrep)) grpinds[[i]] <- rep(1:K, nrep[[i]])[1:clens[[i]]]
    (1:narr)[-which(unlist(grpinds) == iternum)]
  }
# try it out
library("MASS")
data(crabs)
p1c = balKfold.xvspec(5)
inds = p1c( crabs, "sp", 3 )
table(crabs$sp[inds] )
inds2 = p1c( crabs, "sp", 4 )
table(crabs$sp[inds2] )
allc = 1:200
# are test sets disjoint?
intersect(setdiff(allc,inds), setdiff(allc,inds2))



cleanEx()
nameEx("classifierOutput-class")
### * classifierOutput-class

flush(stderr()); flush(stdout())

### Name: classifierOutput-class
### Title: Class "classifierOutput"
### Aliases: classifierOutput-class RObject,classifierOutput-method RObject
###   trainInd,classifierOutput-method trainInd
###   show,classifierOutput-method testScores,classifierOutput-method
###   trainScores,classifierOutput-method
###   predictions,classifierOutput-method predictions
###   predScores,classifierOutput-method predScores
###   predScore,classifierOutput-method predScore
###   testPredictions,classifierOutput-method testPredictions
###   trainPredictions trainPredictions,classifierOutput-method
###   fsHistory,classifierOutput-method testScores trainScores
###   testPredictions
### Keywords: classes

### ** Examples

showClass("classifierOutput")
library(golubEsets)
data(Golub_Train) # now cross-validate a neural net
set.seed(1234)
xv5 = xvalSpec("LOG", 5, balKfold.xvspec(5))
m2 = MLearn(ALL.AML~., Golub_Train[1000:1050,], nnetI, xv5, 
   size=5, decay=.01, maxit=1900 )
testScores(RObject(m2)[[1]]$mlans)
alls = lapply(RObject(m2), function(x) testScores(x$mlans))



cleanEx()
nameEx("clusteringOutput-class")
### * clusteringOutput-class

flush(stderr()); flush(stdout())

### Name: clusteringOutput-class
### Title: container for clustering outputs in uniform structure
### Aliases: clusteringOutput-class RObject,clusteringOutput-method
###   plot,clusteringOutput,ANY-method show,clusteringOutput-method
###   show,clusteringSchema-method getConverter,clusteringSchema-method
###   getDist,clusteringSchema-method getConverter getDist
###   clusteringSchema-class prcompObj-class silhouette-class prcomp-class
### Keywords: classes

### ** Examples

showClass("clusteringOutput")



cleanEx()
nameEx("confuMat-methods")
### * confuMat-methods

flush(stderr()); flush(stdout())

### Name: confuMat-methods
### Title: Compute the confusion matrix for a classifier.
### Aliases: confuMat confuMat-methods confuMat,classifierOutput-method
###   confuMat,classifierOutput,character-method
###   confuMat,classifierOutput,missing-method
###   confuMat,classifierOutput,numeric-method
### Keywords: methods classif

### ** Examples

library(golubEsets)
data(Golub_Merge)
smallG <- Golub_Merge[101:150,]
k1 <- MLearn(ALL.AML~., smallG, knnI(k=1), 1:30)
confuMat(k1)
confuMat(k1, "train")



cleanEx()
nameEx("confuTab")
### * confuTab

flush(stderr()); flush(stdout())

### Name: confuTab
### Title: Compute confusion tables for a confusion matrix.
### Aliases: confuTab

### ** Examples

## the confusion matrix
cm <- table(iris$Species, sample(iris$Species))
## the 3 confusion tables
(ct <- confuTab(cm))



cleanEx()
nameEx("fs.absT")
### * fs.absT

flush(stderr()); flush(stdout())

### Name: fs.absT
### Title: support for feature selection in cross-validation
### Aliases: fs.absT fs.probT fs.topVariance
### Keywords: models

### ** Examples

library("MASS")
data(crabs)
# we will demonstrate this procedure with the crabs data.
# first, create the closure to pick 3 features
demFS = fs.absT(3)
# run it on the entire dataset with features excluding sex
demFS(sp~.-sex, crabs)
# emulate cross-validation by excluding last 50 records
demFS(sp~.-sex, crabs[1:150,])
# emulate cross-validation by excluding first 50 records -- different features retained
demFS(sp~.-sex, crabs[51:200,])



cleanEx()
nameEx("fsHistory")
### * fsHistory

flush(stderr()); flush(stdout())

### Name: fsHistory
### Title: extract history of feature selection for a cross-validated
###   machine learner
### Aliases: fsHistory
### Keywords: models

### ** Examples

data(iris)
iris2 = iris[ iris$Species %in% levels(iris$Species)[1:2], ]
iris2$Species = factor(iris2$Species) # drop unused levels
x1 = MLearn(Species~., iris2, ldaI, xvalSpec("LOG", 3, 
   balKfold.xvspec(3), fs.absT(3)))
fsHistory(x1)



cleanEx()
nameEx("hclustWidget")
### * hclustWidget

flush(stderr()); flush(stdout())

### Name: hclustWidget
### Title: shiny-oriented GUI for cluster or classifier exploration
### Aliases: hclustWidget mlearnWidget
### Keywords: models

### ** Examples

# should run with example(hclustWidget, ask=FALSE)
if (interactive()) {
 library(shiny)
 library(MASS)
 data(crabs)
 cr = data.matrix(crabs[,-c(1:3)])
 au = crabs[,1:3]
 show(hclustWidget(cr, auxdf=au))
## must use stop widget button to proceed
  library(ALL)
  library(hgu95av2.db)
  data(ALL)
  show(mlearnWidget(ALL[1:500,], mol.biol~.))
 }



cleanEx()
nameEx("learnerSchema-class")
### * learnerSchema-class

flush(stderr()); flush(stdout())

### Name: learnerSchema-class
### Title: Class "learnerSchema" - convey information on a machine learning
###   function to the MLearn wrapper
### Aliases: learnerSchema-class nonstandardLearnerSchema-class
###   show,learnerSchema-method
### Keywords: classes

### ** Examples

showClass("learnerSchema")



cleanEx()
nameEx("performance-analytics")
### * performance-analytics

flush(stderr()); flush(stdout())

### Name: performance-analytics
### Title: Assessing classifier performance
### Aliases: precision-methods precision,classifierOutput,character-method
###   precision,classifierOutput,missing-method
###   precision,classifierOutput,numeric-method
###   precision,table,missing-method precision recall-methods
###   recall,classifierOutput,character-method
###   recall,classifierOutput,missing-method
###   recall,classifierOutput,numeric-method recall,table,missing-method
###   recall sensitivity-methods
###   sensitivity,classifierOutput,character-method
###   sensitivity,classifierOutput,missing-method
###   sensitivity,classifierOutput,numeric-method
###   sensitivity,table,missing-method sensitivity macroF1-methods
###   macroF1,classifierOutput,character-method
###   macroF1,classifierOutput,missing-method
###   macroF1,classifierOutput,numeric-method macroF1,table,missing-method
###   macroF1,numeric,numeric-method macroF1 acc,table-method acc
###   specificity,table-method specificity tp,table-method tp
###   tn,table-method tn fp,table-method fp fn,table-method fn
###   F1,table-method F1
### Keywords: methods

### ** Examples

## the confusion matrix
cm <- table(iris$Species, sample(iris$Species))
tp(cm)
tn(cm)
fp(cm)
fn(cm)
acc(cm)
precision(cm)
recall(cm)
F1(cm)
macroF1(cm)



cleanEx()
nameEx("planarPlot-methods")
### * planarPlot-methods

flush(stderr()); flush(stdout())

### Name: planarPlot-methods
### Title: Methods for Function planarPlot in Package 'MLInterfaces'
### Aliases: planarPlot planarPlot-methods
###   planarPlot,classifierOutput,ExpressionSet,character-method
###   planarPlot,classifierOutput,data.frame,character-method
### Keywords: methods

### ** Examples

library(ALL)
library(hgu95av2.db)
data(ALL)
#
# restrict to BCR/ABL or NEG
#
bio <- which( ALL$mol.biol %in% c("BCR/ABL", "NEG"))
#
# restrict to B-cell
#
isb <- grep("^B", as.character(ALL$BT))
kp <- intersect(bio,isb)
all2 <- ALL[,kp]
#
# sample 2 genes at random
#
set.seed(1234)
ng <- nrow(exprs(all2)) # pick 5 in case any NAs come back
pick <- sample(1:ng, size=5, replace=FALSE)
gg <- all2[pick,]
sym <- unlist(mget(featureNames(gg), hgu95av2SYMBOL))
bad = which(is.na(sym))
if (length(bad)>0) {
  gg = gg[-bad,]
  sym = sym[-bad]
  }
gg = gg[1:2,]
sym = sym[1:2]
featureNames(gg) <- sym
gg$class = factor(ifelse(all2$mol.biol=="NEG", "NEG", "POS"))

cl1 <- which( gg$class == "NEG" )
cl2 <- which( gg$class != "NEG" )
#
# create balanced training sample
#
trainInds <- c( sample(cl1, size=floor(length(cl1)/2) ),
      sample(cl2, size=floor(length(cl2)/2)) )
#
# run rpart
#
tgg <- MLearn(class~., gg, rpartI, trainInds, minsplit=4 )
opar <- par(no.readonly=TRUE)
par(mfrow=c(2,2))
planarPlot( tgg, gg, "class" )
title("rpart")
points(exprs(gg)[1,trainInds], exprs(gg)[2,trainInds], col=ifelse(gg$class[trainInds]=="NEG", "yellow", "black"), pch=16)
#
# run nnet
#
ngg <- MLearn( class~., gg, nnetI, trainInds, size=8 )
planarPlot( ngg, gg, "class" )
points(exprs(gg)[1,trainInds], exprs(gg)[2,trainInds], col=ifelse(gg$class[trainInds]=="NEG", "yellow", "black"), pch=16)
title("nnet")
#
# run knn
#
kgg <- MLearn( class~.,  gg, knnI(k=3,l=1), trainInds)
planarPlot( kgg, gg, "class" )
points(exprs(gg)[1,trainInds], exprs(gg)[2,trainInds], col=ifelse(gg$class[trainInds]=="NEG", "yellow", "black"), pch=16)
title("3-nn")
#
# run svm
#
sgg <- MLearn( class~., gg, svmI, trainInds )
planarPlot( sgg, gg, "class" )
points(exprs(gg)[1,trainInds], exprs(gg)[2,trainInds], col=ifelse(gg$class[trainInds]=="NEG", "yellow", "black"), pch=16)
title("svm")
par(opar)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plspinHcube")
### * plspinHcube

flush(stderr()); flush(stdout())

### Name: plspinHcube
### Title: shiny app for interactive 3D visualization of mlbench hypercube
### Aliases: plspinHcube
### Keywords: models

### ** Examples

if (interactive()) plspinHcube()



cleanEx()
nameEx("predict.classifierOutput")
### * predict.classifierOutput

flush(stderr()); flush(stdout())

### Name: predict.classifierOutput
### Title: Predict method for 'classifierOutput' objects
### Aliases: predict.classifierOutput

### ** Examples

## Not run: 
##D set.seed(1234)
##D data(sample.ExpressionSet)
##D trainInd <- 1:16
##D 
##D clout.svm <- MLearn(type~., sample.ExpressionSet[100:250,], svmI, trainInd)
##D predict(clout.svm, sample.ExpressionSet[100:250,-trainInd])
##D 
##D clout.ksvm <- MLearn(type~., sample.ExpressionSet[100:250,], ksvmI, trainInd)
##D predict(clout.ksvm, sample.ExpressionSet[100:250,-trainInd])
##D 
##D clout.nnet <- MLearn(type~., sample.ExpressionSet[100:250,], nnetI, trainInd, size=3, decay=.01 )
##D predict(clout.nnet, sample.ExpressionSet[100:250,-trainInd])
##D 
##D clout.knn <- MLearn(type~., sample.ExpressionSet[100:250,], knnI(k=3), trainInd)
##D predict(clout.knn, sample.ExpressionSet[100:250,-trainInd],k=1)
##D predict(clout.knn, sample.ExpressionSet[100:250,-trainInd],k=3)
##D 
##D #clout.plsda <- MLearn(type~., sample.ExpressionSet[100:250,], plsdaI, trainInd)
##D #predict(clout.plsda, sample.ExpressionSet[100:250,-trainInd])
##D 
##D clout.nb <- MLearn(type~., sample.ExpressionSet[100:250,], naiveBayesI, trainInd)
##D predict(clout.nb, sample.ExpressionSet[100:250,-trainInd])
##D 
##D # this can fail if training set does not yield sufficient diversity in response vector;
##D # setting seed seems to help with this example, but other applications may have problems
##D #
##D clout.rf <- MLearn(type~., sample.ExpressionSet[100:250,], randomForestI, trainInd)
##D predict(clout.rf, sample.ExpressionSet[100:250,-trainInd])
## End(Not run) # end of dontrun



cleanEx()
nameEx("projectLearnerToGrid")
### * projectLearnerToGrid

flush(stderr()); flush(stdout())

### Name: projectLearnerToGrid
### Title: create learned tesselation of feature space after PC
###   transformation
### Aliases: projectLearnerToGrid
### Keywords: models

### ** Examples

library(mlbench)
# demostrate with 3 dimensional hypercube problem
kk = mlbench.hypercube()
colnames(kk$x) = c("f1", "f2", "f3")
hcu = data.frame(cl=kk$classes, kk$x)
set.seed(1234)
sam = sample(1:nrow(kk$x), size=nrow(kk$x)/2)
ldap = projectLearnerToGrid(cl~., data=hcu, ldaI, 
   sam, predWrapper=function(x)x$class)
plot(ldap)
confuMat(ldap@fittedLearner)
nnetp = projectLearnerToGrid(cl~., data=hcu, nnetI, sam, size=2,
   decay=.01, predExtras=list(type="class"))
plot(nnetp)
confuMat(nnetp@fittedLearner)
#if (requireNamespace("rgl") && interactive()) {
#    learnerIn3D(nnetp)
#    ## customising the rgl plot
#    learnerIn3D(nnetp, size = 10, alpha = 0.1)
#}



cleanEx()
nameEx("projectedLearner-class")
### * projectedLearner-class

flush(stderr()); flush(stdout())

### Name: projectedLearner-class
### Title: Class '"projectedLearner"'
### Aliases: projectedLearner-class learnerIn3D,projectedLearner-method
###   plot,projectedLearner,ANY-method plotOne,projectedLearner-method
###   show,projectedLearner-method plotOne learnerIn3D
### Keywords: classes

### ** Examples

showClass("projectedLearner")



cleanEx()
nameEx("raboostCont-class")
### * raboostCont-class

flush(stderr()); flush(stdout())

### Name: raboostCont-class
### Title: Class "raboostCont" ~~~
### Aliases: raboostCont-class daboostCont-class show,raboostCont-method
### Keywords: classes

### ** Examples

showClass("raboostCont")



cleanEx()
nameEx("varImpStruct-class")
### * varImpStruct-class

flush(stderr()); flush(stdout())

### Name: varImpStruct-class
### Title: Class "varImpStruct" - collect data on variable importance from
###   various machine learning methods
### Aliases: varImpStruct-class plot plot,varImpStruct-method
###   plot,varImpStruct,ANY-method show,varImpStruct-method
###   report,varImpStruct-method report getVarImp
###   getVarImp,classifOutput,logical-method
###   getVarImp,classifierOutput,logical-method
###   getVarImp,classifierOutput,missing-method
### Keywords: classes

### ** Examples

library(golubEsets)
data(Golub_Merge)
library(hu6800.db)
smallG <- Golub_Merge[1001:1060,]
set.seed(1234)
opar=par(no.readonly=TRUE)
par(las=2, mar=c(10,11,5,5))
rf2 <- MLearn(ALL.AML~., smallG, randomForestI, 1:40, importance=TRUE,
 sampsize=table(smallG$ALL.AML[1:40]), mtry=sqrt(ncol(exprs(smallG))))
plot( getVarImp( rf2, FALSE ), n=10, plat="hu6800", toktype="SYMBOL")
par(opar)
report( getVarImp( rf2, FALSE ), n=10, plat="hu6800", toktype="SYMBOL")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("xvalLoop")
### * xvalLoop

flush(stderr()); flush(stdout())

### Name: xvalLoop
### Title: Cross-validation in clustered computing environments
### Aliases: xvalLoop xvalLoop,ANY-method
### Keywords: methods

### ** Examples

## Not run: 
##D library(golubEsets)
##D data(Golub_Merge)
##D smallG <- Golub_Merge[200:250,]
##D 
##D # Evaluation on one node
##D 
##D lk1 <- xval(smallG, "ALL.AML", knnB, xvalMethod="LOO", group=as.integer(0))
##D table(lk1,smallG$ALL.AML)
##D 
##D # Evaluation on several nodes -- a cluster programmer might write the following...
##D 
##D library(snow)
##D setOldClass("spawnedMPIcluster")
##D 
##D setMethod("xvalLoop", signature( cluster = "spawnedMPIcluster"),
##D ## use the function returned below to evalutae
##D ## the central cross-validation loop in xval
##D function( cluster, ... ) {
##D     clusterExportEnv <- function (cl, env = .GlobalEnv)
##D     {
##D         unpackEnv <- function(env) {
##D             for ( name in ls(env) ) assign(name, get(name, env), .GlobalEnv )
##D             NULL
##D         }
##D         clusterCall(cl, unpackEnv, env)
##D     }
##D     function(X, FUN, ...) { # this gets returned to xval
##D         ## send all visible variables from the parent (i.e., xval) frame
##D         clusterExportEnv( cluster, parent.frame(1) )
##D         parLapply( cluster, X, FUN, ... )
##D     }
##D })
##D 
##D # ... and use the cluster like this...
##D 
##D cl <- makeCluster(2, "MPI")
##D clusterEvalQ(cl, library(MLInterfaces))
##D 
##D lk1 <- xval(smallG, "ALL.AML", knnB, xvalMethod="LOO", group=as.integer(0), cluster = cl)
##D table(lk1,smallG$ALL.AML)
## End(Not run)


cleanEx()
nameEx("xvalSpec")
### * xvalSpec

flush(stderr()); flush(stdout())

### Name: xvalSpec
### Title: container for information specifying a cross-validated machine
###   learning exercise
### Aliases: xvalSpec xvalSpec-class
### Keywords: models

### ** Examples

library("MASS")
data(crabs)
set.seed(1234)
#
# demonstrate cross validation
#
nn1cv = MLearn(sp~CW+RW, data=crabs, nnetI, xvalSpec("LOG",
   5, balKfold.xvspec(5)), size=3, decay=.01 )
nn1cv
confuMat(nn1cv)
names(RObject(nn1cv)[[1]])
RObject(RObject(nn1cv)[[1]]$mlans)



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
