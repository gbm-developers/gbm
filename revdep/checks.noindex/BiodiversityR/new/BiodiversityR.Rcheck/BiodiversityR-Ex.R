pkgname <- "BiodiversityR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('BiodiversityR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("BCI.env")
### * BCI.env

flush(stderr()); flush(stdout())

### Name: BCI.env
### Title: Barro Colorado Island Quadrat Descriptions
### Aliases: BCI.env
### Keywords: datasets

### ** Examples

data(BCI.env)



cleanEx()
nameEx("CAPdiscrim")
### * CAPdiscrim

flush(stderr()); flush(stdout())

### Name: CAPdiscrim
### Title: Canonical Analysis of Principal Coordinates based on
###   Discriminant Analysis
### Aliases: CAPdiscrim
### Keywords: multivariate

### ** Examples

library(vegan)
library(MASS)
data(dune)
data(dune.env)
Ordination.model1 <- CAPdiscrim(dune~Management, data=dune.env,
    dist="bray", axes=2, m=0, add=FALSE)
Ordination.model1
plot1 <- ordiplot(Ordination.model1, type="none")
ordisymbol(plot1, dune.env, "Management", legend=TRUE)

# plot change in classification success against m
plot(seq(1:14), rep(-1000, 14), xlim=c(1, 14), ylim=c(0, 100), xlab="m", 
    ylab="classification success (percent)", type="n")
for (mseq in 1:14) {
    CAPdiscrim.result <- CAPdiscrim(dune~Management, data=dune.env, 
        dist="bray", axes=2, m=mseq)
    points(mseq, CAPdiscrim.result$percent)
}

#




cleanEx()
nameEx("NMSrandom")
### * NMSrandom

flush(stderr()); flush(stdout())

### Name: NMSrandom
### Title: Calculate the NMS Result with the Smallest Stress from Various
###   Random Starts
### Aliases: NMSrandom
### Keywords: multivariate

### ** Examples

library(vegan)
library(MASS)
data(dune)
distmatrix <- vegdist(dune)
Ordination.model1 <- NMSrandom(distmatrix,perm=100,k=2)
Ordination.model1 <- add.spec.scores(Ordination.model1,dune, 
    method='wa.scores')
Ordination.model1



cleanEx()
nameEx("PCAsignificance")
### * PCAsignificance

flush(stderr()); flush(stdout())

### Name: PCAsignificance
### Title: PCA Significance
### Aliases: PCAsignificance ordiequilibriumcircle
### Keywords: multivariate

### ** Examples

library(vegan)
data(dune)
Ordination.model1 <- rda(dune)
PCAsignificance(Ordination.model1)
plot1 <- ordiplot(Ordination.model1, choices=c(1,2), scaling=1)
ordiequilibriumcircle(Ordination.model1,plot1)



cleanEx()
nameEx("accumresult")
### * accumresult

flush(stderr()); flush(stdout())

### Name: accumresult
### Title: Alternative Species Accumulation Curve Results
### Aliases: accumresult accumplot accumcomp
### Keywords: multivariate

### ** Examples

library(vegan)
data(dune.env)
data(dune)
dune.env$site.totals <- apply(dune,1,sum)
Accum.1 <- accumresult(dune, y=dune.env, scale='site.totals', method='exact', conditioned=TRUE)
Accum.1
accumplot(Accum.1)
accumcomp(dune, y=dune.env, factor='Management', method='exact', legend=FALSE, conditioned=TRUE)
## CLICK IN THE GRAPH TO INDICATE WHERE THE LEGEND NEEDS TO BE PLACED FOR
## OPTION WHERE LEGEND=TRUE (DEFAULT).



cleanEx()
nameEx("add.spec.scores")
### * add.spec.scores

flush(stderr()); flush(stdout())

### Name: add.spec.scores
### Title: Add Species Scores to Unconstrained Ordination Results
### Aliases: add.spec.scores
### Keywords: multivariate

### ** Examples

library(vegan)
data(dune)
distmatrix <-vegdist(dune, method="euc")
# Principal coordinates analysis with 19 axes to estimate total variance
Ordination.model1 <- cmdscale (distmatrix, k=19, eig=TRUE, add=FALSE)
# Change scores for second axis
Ordination.model1$points[,2] <- -1.0 * Ordination.model1$points[,2]
Ordination.model1 <- add.spec.scores(Ordination.model1, dune, 
    method='pcoa.scores', Rscale=TRUE, scaling=1, multi=1)
# Compare Ordination.model1 with PCA
Ordination.model2 <- rda(dune, scale=FALSE)
#
par(mfrow=c(1,2))
ordiplot(Ordination.model1, type="text")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
plot(Ordination.model2, type="text", scaling=1)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("balanced.specaccum")
### * balanced.specaccum

flush(stderr()); flush(stdout())

### Name: balanced.specaccum
### Title: Balanced Species Accumulation Curves
### Aliases: balanced.specaccum
### Keywords: multivariate

### ** Examples

library(vegan)
data(dune.env)
data(dune)

# not balancing species accumulation
Accum.orig <- specaccum(dune)
Accum.orig

# randomly sample 3 quadrats from each stratum of Management
Accum.1 <- balanced.specaccum(dune, strata=dune.env$Management, reps=3)
Accum.1

# scale results by number of trees per quadrat
dune.env$site.totals <- apply(dune,1,sum)
Accum.2 <- balanced.specaccum(dune, strata=dune.env$Management, reps=3, scale=dune.env$site.totals)
Accum.2



cleanEx()
nameEx("caprescale")
### * caprescale

flush(stderr()); flush(stdout())

### Name: caprescale
### Title: Rescaling of Capscale Results to Reflect Total Sums of Squares
###   Of Distance Matrix
### Aliases: caprescale
### Keywords: multivariate

### ** Examples

library(vegan)
library(MASS)
data(dune)
data(dune.env)
Distmatrix.1 <- vegdist(dune,method='bray')
Ordination.model1 <- cmdscale(Distmatrix.1, k=19, eig=TRUE, add=FALSE)
# Sum of all eigenvalues
sum(Ordination.model1$eig)
# [1] 4.395807541512926
sum(Ordination.model1$eig[1:14])
# [1] 4.593946896588808
Distmatrix.2 <- as.matrix(vegdist(Ordination.model1$points[,1:14],method='euc'))
totalsumsquares1 <- sum(Distmatrix.2^2)/(2*20)
# Sum of distances among sites in principal coordinates analysis on axes
# corresponding to positive eigenvalues
totalsumsquares1
# [1] 4.593946896588808
Ordination.model2 <- capscale(dune ~ Management,dune.env,dist='bray', add=FALSE)
# Total sums of positive eigenvalues of the distance-based redundancy analysis
Ordination.model2$CA$tot.chi+Ordination.model2$CCA$tot.chi
# [1] 4.593946896588808
Ordination.model3 <- caprescale(Ordination.model2, verbose=TRUE)
sum1 <- summary(Ordination.model3,axes=17,scaling=1)$constraints
Distmatrix.3 <- as.matrix(vegdist(sum1 ,method='euc'))
totalsumsquares2 <- sum((Distmatrix.3)^2)/(2*20)/19
totalsumsquares2
# [1] 4.593946896588808




cleanEx()
nameEx("crosstabanalysis")
### * crosstabanalysis

flush(stderr()); flush(stdout())

### Name: crosstabanalysis
### Title: Presence-absence Analysis by Cross Tabulation
### Aliases: crosstabanalysis
### Keywords: multivariate

### ** Examples

library(vegan)
data(dune.env)
crosstabanalysis(dune.env,"Manure","Management")



cleanEx()
nameEx("deviancepercentage")
### * deviancepercentage

flush(stderr()); flush(stdout())

### Name: deviancepercentage
### Title: Calculate Percentage and Significance of Deviance Explained by a
###   GLM
### Aliases: deviancepercentage
### Keywords: multivariate

### ** Examples

library(vegan)
data(dune)
data(dune.env)
dune.env$Agrostol <- dune$Agrostol
Count.model1 <- glm(Agrostol ~ Management + A1, family=quasipoisson(link=log), 
    data=dune.env, na.action=na.omit)
summary(Count.model1)
deviancepercentage(Count.model1, dune.env, digits=3)



cleanEx()
nameEx("dist.eval")
### * dist.eval

flush(stderr()); flush(stdout())

### Name: dist.eval
### Title: Distance Matrix Evaluation
### Aliases: dist.eval prepare.bioenv
### Keywords: multivariate

### ** Examples

library(vegan)
data(dune)
dist.eval(dune,"euclidean")
dist.eval(dune,"bray")

## Not run: 
##D data(dune.env)
##D dune.env2 <- dune.env[,c('A1', 'Moisture', 'Manure')]
##D dune.env2$Moisture <- as.numeric(dune.env2$Moisture)
##D dune.env2$Manure <- as.numeric(dune.env2$Manure)
##D sol <- bioenv(dune ~ A1 + Moisture + Manure, dune.env2)
##D sol
##D summary(sol)
##D dune.env3 <- prepare.bioenv(dune.env, as.numeric=c('Moisture', 'Manure'))
##D bioenv(dune, dune.env3)
## End(Not run)





cleanEx()
nameEx("dist.zeroes")
### * dist.zeroes

flush(stderr()); flush(stdout())

### Name: dist.zeroes
### Title: Distance Matrix Transformation
### Aliases: dist.zeroes
### Keywords: multivariate

### ** Examples

library(vegan)
matrix <- array(0,dim=c(5,3))
matrix[4,] <- c(1,2,3)
matrix[5,] <- c(1,0,0)
dist1 <- vegdist(matrix,method="kulc")
dist1
dist2 <- dist.zeroes(matrix,dist1)
dist2



cleanEx()
nameEx("distdisplayed")
### * distdisplayed

flush(stderr()); flush(stdout())

### Name: distdisplayed
### Title: Compare Distance Displayed in Ordination Diagram with Distances
###   of Distance Matrix
### Aliases: distdisplayed
### Keywords: multivariate

### ** Examples

library(vegan)
library(mgcv)
data(dune)
distmatrix <- vegdist(dune,method="kulc")
ordination.model1 <- cmdscale(distmatrix,k=2)
ordiplot1 <- ordiplot(ordination.model1)
distdisplayed(dune,ordiplot=ordiplot1,distx="kulc",plotit=TRUE,
    method="spearman",permutations=100,gam=TRUE)



cleanEx()
nameEx("disttransform")
### * disttransform

flush(stderr()); flush(stdout())

### Name: disttransform
### Title: Community Matrix Transformation
### Aliases: disttransform
### Keywords: multivariate

### ** Examples

library(vegan)
data(dune)
Community.1 <- disttransform(dune, method='hellinger')
Distmatrix.1 <- vegdist(Community.1,method='euclidean')
Distmatrix.1



cleanEx()
nameEx("diversityresult")
### * diversityresult

flush(stderr()); flush(stdout())

### Name: diversityresult
### Title: Alternative Diversity Results
### Aliases: diversityresult diversitycomp diversityvariables
### Keywords: multivariate

### ** Examples


library(vegan)
data(dune.env)
data(dune)

diversityresult(dune, y=NULL, index="Shannon", method="each site", 
    sortit=TRUE, digits=5)
diversityresult(dune, y=dune.env, factor="Management", level="NM", 
    index="Shannon", method="each site", 
    sortit=TRUE, digits=5)
diversityresult(dune, y=NULL, index="Shannon", method="pooled", digits=5)
diversityresult(dune, y=dune.env, factor="Management", level="NM", 
    index="Shannon", method="pooled", digits=5)
diversityresult(dune, y=NULL, index="Shannon", method="mean", 
    digits=5)
diversityresult(dune, y=NULL, index="Shannon", method="sd", 
    digits=5)
diversityresult(dune, y=NULL, index="Shannon", method="jackknife", 
    digits=5)
diversityresult(dune, y=dune.env, factor="Management", level="NM", 
    index="Shannon", method="jackknife", digits=5)

diversitycomp(dune, y=dune.env, factor1="Moisture", index="Shannon",
    method="pooled", sortit=TRUE)
diversitycomp(dune, y=dune.env, factor1="Moisture", index="Shannon",
    method="mean", sortit=TRUE)
diversitycomp(dune, y=dune.env, factor1="Management", index="Shannon",
    method="jackknife", sortit=TRUE)

diversitycomp(dune, y=dune.env, factor1="Management", factor2="Moisture", 
    index="Shannon", method="pooled", digits=6)
diversitycomp(dune, y=dune.env, factor1="Management", factor2="Moisture", 
    index="Shannon", method="mean", digits=6)



cleanEx()
nameEx("ensemble")
### * ensemble

flush(stderr()); flush(stdout())

### Name: ensemble.calibrate.models
### Title: Suitability mapping based on ensembles of modelling algorithms:
###   calibration of models and weights
### Aliases: ensemble.calibrate.models ensemble.calibrate.weights
###   ensemble.calibrate.models.gbm ensemble.calibrate.models.nnet
###   ensemble.drop1 ensemble.formulae ensemble.weights ensemble.strategy
###   ensemble.threshold ensemble.VIF ensemble.pairs

### ** Examples

## Not run: 
##D # based on examples in the dismo package
##D 
##D # get predictor variables
##D library(dismo)
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D # subset based on Variance Inflation Factors
##D predictors <- subset(predictors, subset=c("bio5", "bio6", 
##D     "bio16", "bio17", "biome"))
##D predictors
##D predictors@title <- "predictors"
##D 
##D # presence points
##D presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
##D pres <- read.table(presence_file, header=TRUE, sep=',')[,-1]
##D 
##D # the kfold function randomly assigns data to groups; 
##D # groups are used as calibration (1/4) and training (3/4) data
##D groupp <- kfold(pres, 4)
##D pres_train <- pres[groupp !=  1, ]
##D pres_test <- pres[groupp ==  1, ]
##D 
##D # choose background points
##D background <- randomPoints(predictors, n=1000, extf=1.00)
##D colnames(background)=c('lon', 'lat')
##D groupa <- kfold(background, 4)
##D backg_train <- background[groupa != 1, ]
##D backg_test <- background[groupa == 1, ]
##D 
##D # formulae for random forest and generalized linear model
##D # compare with: ensemble.formulae(predictors, factors=c("biome"))
##D 
##D rfformula <- as.formula(pb ~ bio5+bio6+bio16+bio17)
##D 
##D glmformula <- as.formula(pb ~ bio5 + I(bio5^2) + I(bio5^3) + 
##D     bio6 + I(bio6^2) + I(bio6^3) + bio16 + I(bio16^2) + I(bio16^3) + 
##D     bio17 + I(bio17^2) + I(bio17^3) )
##D 
##D # fit four ensemble models (RF, GLM, BIOCLIM, DOMAIN)
##D ensemble.nofactors <- ensemble.calibrate.models(x=predictors, p=pres_train, a=backg_train, 
##D     pt=pres_test, at=backg_test,
##D     species.name="Bradypus",
##D     ENSEMBLE.tune=TRUE,
##D     ENSEMBLE.min = 0.65,
##D     MAXENT=0, MAXLIKE=0, GBM=0, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=0, GAM=0, 
##D     GAMSTEP=0, MGCV=0, MGCVFIX=0, EARTH=0, RPART=0, NNET=0, FDA=0, 
##D     SVM=0, SVME=0, GLMNET=0,
##D     BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=0,
##D     Yweights="BIOMOD",
##D     evaluations.keep=TRUE, models.keep=TRUE,
##D     RF.formula=rfformula,
##D     GLM.formula=glmformula)
##D 
##D # with option models.keep, all model objects are saved in ensemble object
##D # the same slots can be used to replace model objects with new calibrations
##D ensemble.nofactors$models$RF
##D summary(ensemble.nofactors$models$GLM)
##D ensemble.nofactors$models$BIOCLIM
##D ensemble.nofactors$models$DOMAIN
##D ensemble.nofactors$models$formulae
##D 
##D # evaluations are kept in different slot
##D attributes(ensemble.nofactors$evaluations)
##D plot(ensemble.nofactors$evaluations$RF.T, "ROC")
##D 
##D # fit four ensemble models (RF, GLM, BIOCLIM, DOMAIN) using default formulae
##D # variable 'biome' is not included as explanatory variable
##D # results are provided in a file in the 'outputs' subfolder of the working
##D # directory
##D ensemble.nofactors <- ensemble.calibrate.models(x=predictors,
##D     p=pres_train, a=backg_train, 
##D     pt=pres_test, at=backg_test,
##D     layer.drops="biome",
##D     species.name="Bradypus",
##D     ENSEMBLE.tune=TRUE,
##D     ENSEMBLE.min=0.65,
##D     SINK=TRUE,
##D     MAXENT=0, MAXLIKE=0, GBM=0, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=0, GAM=0, 
##D     GAMSTEP=0, MGCV=0, MGCVFIX=0, EARTH=0, RPART=0, NNET=0, FDA=0, 
##D     SVM=0, SVME=0, GLMNET=0,
##D     BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=0,
##D     Yweights="BIOMOD", 
##D     evaluations.keep=TRUE,
##D     formulae.defaults=TRUE)    
##D 
##D # after fitting the individual algorithms (submodels),
##D # transform predictions with a probit link.
##D ensemble.nofactors <- ensemble.calibrate.models(x=predictors,
##D     p=pres_train, a=backg_train, 
##D     pt=pres_test, at=backg_test,
##D     layer.drops="biome",
##D     species.name="Bradypus",
##D     SINK=TRUE,
##D     ENSEMBLE.tune=TRUE,
##D     ENSEMBLE.min=0.65,
##D     MAXENT=0, MAXLIKE=0, GBM=0, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=0, GAM=0, 
##D     GAMSTEP=0, MGCV=0, MGCVFIX=0, EARTH=0, RPART=0, NNET=0, FDA=0, 
##D     SVM=0, SVME=0, GLMNET=0,
##D     BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=0,
##D     PROBIT=TRUE,
##D     Yweights="BIOMOD", factors="biome",
##D     evaluations.keep=TRUE,
##D     formulae.defaults=TRUE)    
##D 
##D # Instead of providing presence and background locations, provide data.frames.
##D # Because 'biome' is a factor, RasterStack needs to be provided
##D # to check for levels in the Training and Testing data set.
##D TrainData1 <- prepareData(x=predictors, p=pres_train, b=backg_train, 
##D     factors=c("biome"), xy=FALSE)
##D TestData1 <- prepareData(x=predictors, p=pres_test, b=backg_test, 
##D     factors=c("biome"), xy=FALSE)
##D ensemble.factors1 <- ensemble.calibrate.models(x=predictors, 
##D     TrainData=TrainData1, TestData=TestData1,
##D     p=pres_train, a=backg_train, 
##D     pt=pres_test, at=backg_test,
##D     species.name="Bradypus",
##D     SINK=TRUE,
##D     ENSEMBLE.tune=TRUE,
##D     ENSEMBLE.min=0.65,
##D     MAXENT=0, MAXLIKE=1, GBM=1, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=1, GAM=1, 
##D     GAMSTEP=1, MGCV=1, MGCVFIX=0, EARTH=1, RPART=1, NNET=1, FDA=1, 
##D     SVM=1, SVME=1, GLMNET=1,
##D     BIOCLIM.O=1, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
##D     Yweights="BIOMOD", factors="biome",
##D     evaluations.keep=TRUE)
##D 
##D # compare different methods of calculating ensembles
##D ensemble.factors2 <- ensemble.calibrate.models(x=predictors, 
##D     TrainData=TrainData1, TestData=TestData1,
##D     p=pres_train, a=backg_train, 
##D     pt=pres_test, at=backg_test,
##D     species.name="Bradypus",
##D     SINK=TRUE,
##D     ENSEMBLE.tune=TRUE,
##D     MAXENT=1, MAXLIKE=1, GBM=1, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=1, GAM=1, 
##D     GAMSTEP=1, MGCV=1, MGCVFIX=1, EARTH=1, RPART=1, NNET=1, FDA=1, 
##D     SVM=1, SVME=1, BIOCLIM.O=1, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
##D     ENSEMBLE.best=c(4:10), ENSEMBLE.exponent=c(1, 2, 3),
##D     Yweights="BIOMOD", factors="biome",
##D     evaluations.keep=TRUE)
##D 
##D # test performance of different suitability models
##D # data are split in 4 subsets, each used once for evaluation
##D ensemble.nofactors2 <- ensemble.calibrate.weights(x=predictors, 
##D     p=pres, a=background, k=4, 
##D     species.name="Bradypus",
##D     SINK=TRUE, PROBIT=TRUE,
##D     MAXENT=1, MAXLIKE=1, GBM=1, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=1, GAM=1, 
##D     GAMSTEP=1, MGCV=1, MGCVFIX=1, EARTH=1, RPART=1, NNET=1, FDA=1, 
##D     SVM=1, SVME=1, BIOCLIM.O=1, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
##D     ENSEMBLE.tune=TRUE,
##D     ENSEMBLE.best=0, ENSEMBLE.exponent=c(1, 2, 3),
##D     ENSEMBLE.min=0.7,
##D     Yweights="BIOMOD", 
##D     formulae.defaults=TRUE)
##D ensemble.nofactors2$AUC.table
##D 
##D # test the result of leaving out one of the variables from the model
##D # note that positive differences indicate that the model without the variable 
##D # has higher AUC than the full model
##D ensemble.variables <- ensemble.drop1(x=predictors, 
##D     p=pres, a=background, k=4,
##D     species.name="Bradypus",
##D     SINK=TRUE,
##D     difference=TRUE,
##D     VIF=TRUE, PROBIT=TRUE,
##D     MAXENT=0, MAXLIKE=0, GBM=1, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=1, GAM=1, 
##D     GAMSTEP=1, MGCV=1, MGCVFIX=1, EARTH=1, RPART=1, NNET=1, FDA=1, 
##D     SVM=1, SVME=1, GLMNET=1,
##D     BIOCLIM.O=0, BIOCLIM=0, DOMAIN=1, MAHAL=0, MAHAL01=1,
##D     ENSEMBLE.tune=TRUE,
##D     ENSEMBLE.best=0, ENSEMBLE.exponent=c(1, 2, 3),
##D     ENSEMBLE.min=0.7,
##D     Yweights="BIOMOD", factors="biome")
##D ensemble.variables
##D 
##D # use function ensemble.VIF to select a subset of variables
##D # factor variables are not handled well by the function
##D # and therefore factors are removed
##D # however, one can check for factors with car::vif through
##D # the ensemble.calibrate.models function
##D # VIF.analysis$var.drops can be used as input for ensemble.calibrate.models or
##D # ensemble.calibrate.weights
##D 
##D predictors <- stack(predictor.files)
##D predictors <- subset(predictors, subset=c("bio1", "bio5", "bio6", "bio8", 
##D     "bio12", "bio16", "bio17", "biome"))
##D 
##D ensemble.pairs(predictors)
##D 
##D VIF.analysis <- ensemble.VIF(predictors, factors="biome")
##D VIF.analysis
##D # alternative solution where bio1 and bio12 are kept
##D VIF.analysis <- ensemble.VIF(predictors, factors="biome", 
##D     keep=c("bio1", "bio12"))
##D VIF.analysis
## End(Not run)



cleanEx()
nameEx("ensemble.analogue")
### * ensemble.analogue

flush(stderr()); flush(stdout())

### Name: ensemble.analogue
### Title: Climate analogues from climatic distance raster layers.
### Aliases: ensemble.analogue ensemble.analogue.object

### ** Examples


## Not run: 
##D # get predictor variables
##D library(dismo)
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D predictors <- subset(predictors, subset=c("bio1", "bio5", "bio6", "bio7", "bio8", 
##D     "bio12", "bio16", "bio17"))
##D predictors
##D predictors@title <- "base"
##D 
##D # instead of searching for current analogue of future climate conditions,
##D # search for analogue in southern hemisphere
##D future.stack <- stack(crop(predictors, y=extent(-125, -32, 0, 40)))
##D future.stack@title <- "north"
##D current.stack <- stack(crop(predictors, y=extent(-125, -32, -56, 0)))
##D current.stack@title <- "south"
##D 
##D # reference location in Florida
##D # in this case future.stack and current.stack are both current
##D ref.loc <- data.frame(t(c(-80.19, 25.76)))
##D names(ref.loc) <- c("lon", "lat")
##D 
##D # climate analogue analysis based on the Mahalanobis distance
##D Florida.object.mahal <- ensemble.analogue.object(ref.location=ref.loc, 
##D     future.stack=future.stack, current.stack=current.stack, 
##D     name="FloridaMahal", method="mahal", an=10000)
##D Florida.object.mahal
##D 
##D Florida.analogue.mahal <- ensemble.analogue(x=current.stack, 
##D     analogue.object=Florida.object.mahal, analogues=50)
##D Florida.analogue.mahal
##D 
##D # climate analogue analysis based on the Euclidean distance and dividing each variable by the sd
##D Florida.object.sd <- ensemble.analogue.object(ref.location=ref.loc, 
##D     future.stack=future.stack, current.stack=current.stack, 
##D     name="FloridaSD", method="sd", z=2)
##D Florida.object.sd
##D 
##D Florida.analogue.sd <- ensemble.analogue(x=current.stack, 
##D     analogue.object=Florida.object.sd, analogues=50)
##D Florida.analogue.sd
##D 
##D # plot analogues on climatic distance maps
##D par(mfrow=c(1,2))
##D analogue.file <- paste(getwd(), "//ensembles//analogue//FloridaMahal_south_analogue.grd", sep="")
##D plot(raster(analogue.file), main="Mahalanobis climatic distance")
##D points(Florida.analogue.sd[3:50, "lat"] ~ Florida.analogue.sd[3:50, "lon"], 
##D     pch=1, col="red", cex=1)
##D points(Florida.analogue.mahal[3:50, "lat"] ~ Florida.analogue.mahal[3:50, "lon"], 
##D     pch=3, col="black", cex=1)
##D points(Florida.analogue.mahal[2, "lat"] ~ Florida.analogue.mahal[2, "lon"], 
##D     pch=22, col="blue", cex=2)
##D legend(x="topright", legend=c("closest", "Mahalanobis", "SD"), pch=c(22, 3 , 1), 
##D     col=c("blue" , "black", "red"))
##D 
##D analogue.file <- paste(getwd(), "//ensembles//analogue//FloridaSD_south_analogue.grd", sep="")
##D plot(raster(analogue.file), main="Climatic distance normalized by standard deviation")
##D points(Florida.analogue.mahal[3:50, "lat"] ~ Florida.analogue.mahal[3:50, "lon"], 
##D     pch=3, col="black", cex=1)
##D points(Florida.analogue.sd[3:50, "lat"] ~ Florida.analogue.sd[3:50, "lon"], 
##D     pch=1, col="red", cex=1)
##D points(Florida.analogue.sd[2, "lat"] ~ Florida.analogue.sd[2, "lon"], 
##D     pch=22, col="blue", cex=2)
##D legend(x="topright", legend=c("closest", "Mahalanobis", "SD"), pch=c(22, 3 , 1), 
##D     col=c("blue" , "black", "red"))
##D par(mfrow=c(1,1))
## End(Not run)




cleanEx()
nameEx("ensemble.batch")
### * ensemble.batch

flush(stderr()); flush(stdout())

### Name: ensemble.batch
### Title: Suitability mapping based on ensembles of modelling algorithms:
###   batch processing
### Aliases: ensemble.batch ensemble.mean ensemble.plot

### ** Examples

## Not run: 
##D # based on examples in the dismo package
##D 
##D # get predictor variables
##D library(dismo)
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D # subset based on Variance Inflation Factors
##D predictors <- subset(predictors, subset=c("bio5", "bio6", 
##D     "bio16", "bio17", "biome"))
##D predictors
##D predictors@title <- "base"
##D 
##D # presence points
##D presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
##D pres <- read.table(presence_file, header=TRUE, sep=',')
##D pres[,1] <- rep("Bradypus", nrow(pres))
##D 
##D # choose background points
##D background <- randomPoints(predictors, n=1000, extf = 1.00)
##D 
##D # north and south for new predictions (as if new climates)
##D ext2 <- extent(-90, -32, 0, 23)
##D predictors2 <- crop(predictors, y=ext2)
##D predictors2 <- stack(predictors2)
##D predictors2@title <- "north"
##D 
##D ext3 <- extent(-90, -32, -33, 0)
##D predictors3 <- crop(predictors, y=ext3)
##D predictors3 <- stack(predictors3)
##D predictors3@title <- "south"
##D 
##D # fit 3 ensembles with batch processing, choosing the best ensemble model based on the 
##D # average weights of 4-fold split of calibration and testing data
##D # final models use all available presence data and average weights determined by the 
##D # ensemble.calibrate.weights function (called internally)
##D # batch processing can handle several species by using 3-column species.presence and 
##D # species.absence data sets
##D # note that these calculations can take a while
##D 
##D ensemble.nofactors <- ensemble.batch(x=predictors, 
##D     xn=c(predictors, predictors2, predictors3),
##D     species.presence=pres, 
##D     species.absence=background, 
##D     k.splits=4, k.test=0, 
##D     n.ensembles=3, 
##D     SINK=TRUE, 
##D     layer.drops=c("biome"),
##D     ENSEMBLE.best=0, ENSEMBLE.exponent=c(1, 2, 3), 
##D     ENSEMBLE.min=0.7,
##D     MAXENT=1, MAXLIKE=1, GBM=1, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=1, GAM=1, 
##D     GAMSTEP=1, MGCV=1, MGCVFIX=1, EARTH=1, RPART=1, NNET=1, FDA=1, 
##D     SVM=1, SVME=1, GLMNET=1,
##D     BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
##D     PROBIT=TRUE,
##D     Yweights="BIOMOD",
##D     formulae.defaults=TRUE)
##D 
##D # summaries for the 3 ensembles for the species
##D # summaries are based on files in folders ensemble/suitability, 
##D # ensemble/presence and ensemble/count
##D # ensemble.mean is used internally in ensemble.batch
##D 
##D ensemble.mean(RASTER.species.name="Bradypus", RASTER.stack.name="base",
##D     p=pres, a=background)
##D 
##D # plot mean suitability without specifying colours
##D plot1 <- ensemble.plot(RASTER.species.name="Bradypus", RASTER.stack.name="base",
##D     plot.method="consensussuitability",
##D     p=pres, a=background, abs.breaks=4, pres.breaks=9)
##D plot1
##D 
##D # only colour the areas where species is predicted to be present
##D # option is invoked by having no absence breaks
##D # same colourscheme as \url{http://www.worldagroforestry.org/atlas-central-america}
##D LAatlascols <- grDevices::colorRampPalette(c("#FFFF80", "#38E009","#1A93AB", "#0C1078"))
##D plot2 <- ensemble.plot(RASTER.species.name="Bradypus", RASTER.stack.name="base",
##D     plot.method="consensussuitability",
##D     p=pres, a=background, abs.breaks=0, pres.breaks=9, pres.col=LAatlascols(8))
##D plot2
##D 
##D # only colour the areas where species is predicted to be present
##D # option is invoked by only setting one colour for absence-presence
##D plot3 <- ensemble.plot(RASTER.species.name="Bradypus", RASTER.stack.name="base",
##D     plot.method="consensuspresence",
##D     absencePresence.col=c("#90EE90"))
##D 
##D # only colour presence area by specifying colours > 0
##D plot4 <- ensemble.plot(RASTER.species.name="Bradypus", RASTER.stack.name="base",
##D     plot.method="consensuscount",
##D     count.col=LAatlascols(3))
##D 
##D 
##D 
## End(Not run)




cleanEx()
nameEx("ensemble.bioclim")
### * ensemble.bioclim

flush(stderr()); flush(stdout())

### Name: ensemble.bioclim
### Title: Suitability mapping based on the BIOCLIM algorithm
### Aliases: ensemble.bioclim ensemble.bioclim.object

### ** Examples


## Not run: 
##D # get predictor variables
##D library(dismo)
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D # subset based on Variance Inflation Factors
##D predictors <- subset(predictors, subset=c("bio5", "bio6", 
##D     "bio16", "bio17", "biome"))
##D predictors
##D predictors@title <- "base"
##D 
##D # presence points
##D presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
##D pres <- read.table(presence_file, header=TRUE, sep=',')[,-1]
##D 
##D background <- dismo::randomPoints(predictors, n=100)
##D colnames(background)=c('lon', 'lat')
##D 
##D pres.dataset <- data.frame(extract(predictors, y=pres))
##D names(pres.dataset) <- names(predictors)
##D pres.dataset$biome <- as.factor(pres.dataset$biome)
##D 
##D Bradypus.bioclim <- ensemble.bioclim.object(predictors, quantiles=T, 
##D     p=pres, factors="biome", species.name="Bradypus")
##D Bradypus.bioclim
##D # obtain the same results with a data.frame
##D Bradypus.bioclim2 <- ensemble.bioclim.object(pres.dataset, quantiles=T, 
##D     species.name="Bradypus")
##D Bradypus.bioclim2
##D # obtain results for entire rasterStack
##D Bradypus.bioclim3 <- ensemble.bioclim.object(predictors, p=NULL, quantiles=T, 
##D     factors="biome", species.name="America")
##D Bradypus.bioclim3
##D 
##D ensemble.bioclim(x=predictors, bioclim.object=Bradypus.bioclim, KML.out=T)
##D ensemble.bioclim(x=predictors, bioclim.object=Bradypus.bioclim3, KML.out=T)
##D 
##D par.old <- graphics::par(no.readonly=T)
##D graphics::par(mfrow=c(1,2))
##D 
##D rasterfull1 <- paste("ensembles//Bradypus_base_BIOCLIM_orig", sep="")
##D raster::plot(raster(rasterfull1), breaks=c(-0.1, 0, 0.5, 1), 
##D     col=c("grey", "blue", "green"), main="original method")
##D rasterfull2 <- paste("ensembles//America_base_BIOCLIM_orig", sep="")
##D raster::plot(raster(rasterfull2), breaks=c(-0.1, 0, 0.5, 1), 
##D     col=c("grey", "blue", "green"), main="America")
##D 
##D graphics::par(par.old)
##D 
##D # compare with implementation bioclim in dismo
##D bioclim.dismo <- bioclim(predictors, p=pres)
##D rasterfull2 <- paste("ensembles//Bradypus_base_BIOCLIM_dismo", sep="")
##D raster::predict(object=predictors, model=bioclim.dismo, na.rm=TRUE, 
##D     filename=rasterfull2, progress='text', overwrite=TRUE)
##D 
##D par.old <- graphics::par(no.readonly=T)
##D graphics::par(mfrow=c(1,2))
##D 
##D raster::plot(raster(rasterfull1), breaks=c(-0.1, 0, 0.5, 1), 
##D     col=c("grey", "blue", "green"), main="original method")
##D raster::plot(raster(rasterfull2), main="dismo method")
##D 
##D graphics::par(par.old)
##D 
##D # use dummy variables to deal with factors
##D predictors <- stack(predictor.files)
##D biome.layer <- predictors[["biome"]]
##D biome.layer
##D ensemble.dummy.variables(xcat=biome.layer, most.frequent=0, freq.min=1,
##D     overwrite=TRUE)
##D 
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D predictors.dummy <- subset(predictors, subset=c("biome_1", "biome_2",  "biome_3",  
##D     "biome_4", "biome_5", "biome_7",  "biome_8",  "biome_9", "biome_10", 
##D     "biome_12", "biome_13", "biome_14"))
##D predictors.dummy
##D predictors.dummy@title <- "base_dummy"
##D 
##D Bradypus.dummy <- ensemble.bioclim.object(predictors.dummy, quantiles=T, 
##D     p=pres, species.name="Bradypus")
##D Bradypus.dummy
##D ensemble.bioclim(x=predictors.dummy, bioclim.object=Bradypus.dummy, KML.out=F)
##D 
##D par.old <- graphics::par(no.readonly=T)
##D graphics::par(mfrow=c(1,2))
##D 
##D rasterfull3 <- paste("ensembles//Bradypus_base_dummy_BIOCLIM_orig", sep="")
##D raster::plot(raster(rasterfull1), breaks=c(-0.1, 0, 0.5, 1), col=c("grey", "blue", "green"), 
##D     main="numeric predictors")
##D raster::plot(raster(rasterfull3), breaks=c(-0.1, 0, 0.5, 1), col=c("grey", "blue", "green"), 
##D     main="dummy predictors")
##D 
##D graphics::par(par.old)
## End(Not run)




cleanEx()
nameEx("ensemble.bioclim.graph")
### * ensemble.bioclim.graph

flush(stderr()); flush(stdout())

### Name: ensemble.bioclim.graph
### Title: Graphs of bioclimatic ranges of species and climates
### Aliases: ensemble.bioclim.graph ensemble.bioclim.graph.data

### ** Examples

## Not run: 
##D 
##D # get predictor variables
##D library(dismo)
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D # subset based on Variance Inflation Factors
##D predictors <- subset(predictors, subset=c("bio5", "bio6", 
##D     "bio16", "bio17", "biome"))
##D predictors
##D predictors@title <- "base"
##D 
##D # presence points
##D presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
##D pres <- read.table(presence_file, header=TRUE, sep=',')[,-1]
##D 
##D # climates for north and south (use same process for future climates)
##D ext2 <- extent(-90, -32, 0, 23)
##D predictors2 <- crop(predictors, y=ext2)
##D predictors2 <- stack(predictors2)
##D predictors2@title <- "north"
##D 
##D ext3 <- extent(-90, -32, -33, 0)
##D predictors3 <- crop(predictors, y=ext3)
##D predictors3 <- stack(predictors3)
##D predictors3@title <- "south"
##D 
##D graph.data1 <- ensemble.bioclim.graph.data(predictors, p=pres, 
##D     factors="biome", species.climate.name="Bradypus")
##D graph.data2 <- ensemble.bioclim.graph.data(predictors, p=NULL, 
##D     factors="biome", species.climate.name="baseline")
##D graph.data3 <- ensemble.bioclim.graph.data(predictors2, p=NULL, 
##D     factors="biome", species.climate.name="north")
##D graph.data4 <- ensemble.bioclim.graph.data(predictors3, p=NULL, 
##D     factors="biome", species.climate.name="south")
##D graph.data.all <- rbind(graph.data1, graph.data2, graph.data3, graph.data4)
##D 
##D par.old <- graphics::par(no.readonly=T)
##D graphics::par(mfrow=c(2, 2))
##D 
##D ensemble.bioclim.graph(graph.data.all, focal.var="bio5", 
##D     var.multiply=0.1, cols=c("black", rep("blue", 3)))
##D ensemble.bioclim.graph(graph.data.all, focal.var="bio6", 
##D     var.multiply=0.1, cols=c("black", rep("blue", 3)))
##D ensemble.bioclim.graph(graph.data.all, focal.var="bio16", 
##D     var.multiply=1.0, cols=c("black", rep("blue", 3)))
##D ensemble.bioclim.graph(graph.data.all, focal.var="bio17", 
##D     var.multiply=1.0, cols=c("black", rep("blue", 3)))
##D 
##D graphics::par(par.old)
##D 
## End(Not run)



cleanEx()
nameEx("ensemble.dummy.variables")
### * ensemble.dummy.variables

flush(stderr()); flush(stdout())

### Name: ensemble.dummy.variables
### Title: Suitability mapping based on ensembles of modelling algorithms:
###   handling of categorical data
### Aliases: ensemble.dummy.variables ensemble.accepted.categories
###   ensemble.simplified.categories

### ** Examples

## Not run: 
##D 
##D # get predictor variables
##D library(dismo)
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D biome.layer <- predictors[["biome"]]
##D biome.layer
##D 
##D # create dummy layers for the 5 most frequent factor levels
##D 
##D ensemble.dummy.variables(xcat=biome.layer, most.frequent=5,
##D     overwrite=TRUE)
##D 
##D # check whether dummy variables were created
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D predictors
##D names(predictors)
##D 
##D # once dummy variables were created, avoid using the original categorical data layer
##D predictors <- subset(predictors, subset=c("bio5", "bio6", "bio16", "bio17", 
##D     "biome_1", "biome_2", "biome_7", "biome_8", "biome_13"))
##D predictors
##D predictors@title <- "base"
##D 
##D # presence points
##D presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
##D pres <- read.table(presence_file, header=TRUE, sep=',')[,-1]
##D 
##D # the kfold function randomly assigns data to groups; 
##D # groups are used as calibration (1/5) and training (4/5) data
##D groupp <- kfold(pres, 5)
##D pres_train <- pres[groupp !=  1, ]
##D pres_test <- pres[groupp ==  1, ]
##D 
##D # choose background points
##D background <- randomPoints(predictors, n=1000, extf=1.00)
##D colnames(background)=c('lon', 'lat')
##D groupa <- kfold(background, 5)
##D backg_train <- background[groupa != 1, ]
##D backg_test <- background[groupa == 1, ]
##D 
##D # note that dummy variables with no variation are not used by DOMAIN
##D # note that dummy variables are not used by MAHAL and MAHAL01
##D # (neither are categorical variables)
##D ensemble.nofactors <- ensemble.calibrate.models(x=predictors, p=pres_train, a=backg_train, 
##D     pt=pres_test, at=backg_test,
##D     species.name="Bradypus",
##D     VIF=T,
##D     MAXENT=1, MAXLIKE=1, GBM=1, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=0, GAM=1, 
##D     GAMSTEP=0, MGCV=1, MGCVFIX=0, EARTH=1, RPART=1, NNET=1, FDA=1, 
##D     SVM=1, SVME=1, BIOCLIM.O=1, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
##D     Yweights="BIOMOD", 
##D     dummy.vars=c("biome_1", "biome_2", "biome_7", "biome_8", "biome_13"),
##D     PLOTS=FALSE, evaluations.keep=TRUE)
## End(Not run)




cleanEx()
nameEx("ensemble.ecocrop")
### * ensemble.ecocrop

flush(stderr()); flush(stdout())

### Name: ensemble.ecocrop
### Title: Mapping of novel environmental conditions (areas where some of
###   the environmental conditions are outside the range of environmental
###   conditions of a reference area).
### Aliases: ensemble.ecocrop ensemble.ecocrop.object

### ** Examples


## Not run: 
##D #test with Brazil nut (limits from FAO ecocrop)
##D #temperature: (12) 20-36 (40)
##D #annnual rainfall: (1400) 2400-2800 (3500)
##D 
##D # get predictor variables
##D library(dismo)
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D # subset based on Variance Inflation Factors
##D predictors <- subset(predictors, subset=c("bio5", "bio6", "bio12"))
##D predictors
##D predictors@title <- "base"
##D 
##D Brazil.ecocrop <- ensemble.ecocrop.object(temp.thresholds=c(20, 36, 12, 40), 
##D     rain.thresholds=c(2400, 2800, 1400, 3500), 
##D     annual.temps=FALSE, name="Bertholletia_excelsa")
##D Brazil.ecocrop
##D ensemble.ecocrop(predictors, ecocrop.object=Brazil.ecocrop)
##D 
##D dev.new()
##D par.old <- graphics::par(no.readonly=T)
##D graphics::par(mfrow=c(1,2))
##D 
##D 
##D rasterfull1 <- paste("ensembles//ecocrop//Bertholletia_excelsa_base.grd", sep="")
##D rasterfull1 <- raster(rasterfull1)
##D # raster file saved probabilities as integer values between 0 and 1000
##D rasterfull1 <- rasterfull1/1000
##D raster::plot(rasterfull1, main="Ecocrop suitability")
##D 
##D GBIFloc <- gbif(genus="Bertholletia", species="excelsa", geo=TRUE)
##D GBIFpres <- cbind(GBIFloc$lon, GBIFloc$lat)
##D GBIFpres <- GBIFpres[complete.cases(GBIFpres), ]
##D GBIFpres <- GBIFpres[duplicated(GBIFpres) == FALSE, ]
##D point.suitability <- extract(rasterfull1, y=GBIFpres)
##D point.suitability[is.na(point.suitability)] <- -1
##D 
##D GBIFpres.optimal <- GBIFpres[point.suitability == 1, ]
##D GBIFpres.suboptimal <- GBIFpres[point.suitability < 1 & point.suitability > 0, ]
##D GBIFpres.not <- GBIFpres[point.suitability == 0, ]
##D 
##D raster::plot(rasterfull1, main="GBIF locations", 
##D     sub="blue: optimal, cyan: suboptimal, red: not suitable")
##D bg.legend <- c("blue", "cyan", "red")
##D 
##D points(GBIFpres.suboptimal, pch=21, cex=1.2, bg=bg.legend[2])
##D points(GBIFpres.optimal, pch=21, cex=1.2, bg=bg.legend[1])
##D points(GBIFpres.not, pch=21, cex=1.2, bg=bg.legend[3])
##D 
##D graphics::par(par.old)
## End(Not run)




cleanEx()
nameEx("ensemble.novel")
### * ensemble.novel

flush(stderr()); flush(stdout())

### Name: ensemble.novel
### Title: Mapping of novel environmental conditions (areas where some of
###   the environmental conditions are outside the range of environmental
###   conditions of a reference area).
### Aliases: ensemble.novel ensemble.novel.object

### ** Examples


## Not run: 
##D # get predictor variables
##D library(dismo)
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D predictors <- subset(predictors, subset=c("bio1", "bio5", "bio6", "bio7", "bio8", 
##D     "bio12", "bio16", "bio17"))
##D predictors
##D predictors@title <- "base"
##D 
##D # reference area to calculate environmental ranges
##D ext <- extent(-70, -50, -10, 10)
##D extent.values2 <- c(-70, -50, -10, 10)
##D predictors.current <- crop(predictors, y=ext)
##D predictors.current <- stack(predictors.current)
##D 
##D novel.test <- ensemble.novel.object(predictors.current, name="noveltest")
##D novel.test
##D novel.raster <- ensemble.novel(x=predictors, novel.object=novel.test, KML.out=T)
##D novel.raster
##D 
##D plot(novel.raster)
##D # no novel conditions within reference area
##D rect(extent.values2[1], extent.values2[3], extent.values2[2], extent.values2[4])
##D 
##D # use novel conditions as a simple species suitability mapping method
##D # presence points
##D presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
##D pres <- read.table(presence_file, header=TRUE, sep=',')[,-1]
##D pres.data <- data.frame(extract(predictors, y=pres))
##D 
##D # ranges and maps
##D Bradypus.ranges1 <- ensemble.novel.object(pres.data, name="Bradypus", quantiles=F)
##D Bradypus.ranges1
##D Bradypus.novel1 <- ensemble.novel(x=predictors, novel.object=Bradypus.ranges1, KML.out=T)
##D Bradypus.novel1
##D 
##D par.old <- graphics::par(no.readonly=T)
##D graphics::par(mfrow=c(1,2))
##D 
##D # suitable where there are no novel conditions
##D raster::plot(Bradypus.novel1, breaks=c(-0.1, 0, 1), col=c("green", "grey"), 
##D     main="Suitability mapping using minimum to maximum range")
##D points(pres[, 2] ~ pres[, 1], pch=1, col="red", cex=0.8)
##D 
##D # use 90 percent intervals similar to BIOCLIM methodology
##D Bradypus.ranges2 <- ensemble.novel.object(pres.data, name="BradypusQuantiles", quantiles=T)
##D Bradypus.ranges2
##D Bradypus.novel2 <- ensemble.novel(x=predictors, novel.object=Bradypus.ranges2, KML.out=T)
##D Bradypus.novel2
##D raster::plot(Bradypus.novel2, breaks=c(-0.1, 0, 1), col=c("green", "grey"), 
##D     main="Suitability mapping using quantile range")
##D points(pres[, 2] ~ pres[, 1], pch=1, col="red", cex=0.8)
##D 
##D graphics::par(par.old)
##D 
##D # deal with novel factor levels through dummy variables
##D predictors <- stack(predictor.files)
##D biome.layer <- predictors[["biome"]]
##D biome.layer
##D ensemble.dummy.variables(xcat=biome.layer, most.frequent=0, freq.min=1,
##D     overwrite=TRUE)
##D 
##D predictors.dummy <- stack(predictor.files)
##D predictors.dummy <- subset(predictors.dummy, subset=c("biome_1", "biome_2",  "biome_3",  
##D     "biome_4", "biome_5", "biome_7",  "biome_8",  "biome_9", 
##D     "biome_10", "biome_12", "biome_13", "biome_14"))
##D predictors.dummy
##D predictors.dummy@title <- "base_dummy"
##D 
##D predictors.dummy.current <- crop(predictors.dummy, y=ext)
##D predictors.dummy.current <- stack(predictors.dummy.current)
##D 
##D novel.levels <- ensemble.novel.object(predictors.dummy.current, name="novellevels")
##D novel.levels
##D novel.levels.raster <- ensemble.novel(x=predictors.dummy, novel.object=novel.levels, 
##D     KML.out=T)
##D novel.levels.raster
##D 
##D novel.levels.quantiles <- ensemble.novel.object(predictors.dummy.current, quantiles=TRUE,
##D     name="novellevels_quantiles")
##D novel.levels.quantiles
##D novel.levels.quantiles.raster <- ensemble.novel(x=predictors.dummy, 
##D     novel.object=novel.levels.quantiles, KML.out=T)
##D novel.levels.quantiles.raster
##D 
##D # difference in ranges for variables with low frequencies
##D background <- dismo::randomPoints(predictors.dummy.current, n=10000, p=NULL, excludep=F)
##D extract.data <- extract(predictors.dummy.current, y=background)
##D colSums(extract.data)/sum(extract.data)*100
##D novel.levels
##D novel.levels.quantiles
##D 
##D par.old <- graphics::par(no.readonly=T)
##D graphics::par(mfrow=c(1,2))
##D raster::plot(novel.levels.raster, breaks=c(-0.1, 0, 1), col=c("grey", "green"), 
##D     main="novel outside minimum to maximum range")
##D rect(extent.values2[1], extent.values2[3], extent.values2[2], extent.values2[4])
##D raster::plot(novel.levels.quantiles.raster, breaks=c(-0.1, 0, 1), col=c("grey", "green"), 
##D     main="novel outside quantile range")
##D rect(extent.values2[1], extent.values2[3], extent.values2[2], extent.values2[4])
##D graphics::par(par.old)
##D 
## End(Not run)



cleanEx()
nameEx("ensemble.raster")
### * ensemble.raster

flush(stderr()); flush(stdout())

### Name: ensemble.raster
### Title: Suitability mapping based on ensembles of modelling algorithms:
###   consensus mapping
### Aliases: ensemble.raster ensemble.habitat.change ensemble.area

### ** Examples

## Not run: 
##D # based on examples in the dismo package
##D 
##D # get predictor variables
##D library(dismo)
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D # subset based on Variance Inflation Factors
##D predictors <- subset(predictors, subset=c("bio5", "bio6", 
##D     "bio16", "bio17"))
##D predictors
##D predictors@title <- "base"
##D 
##D # presence points
##D # presence points
##D presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
##D pres <- read.table(presence_file, header=TRUE, sep=',')[,-1]
##D 
##D # choose background points
##D background <- randomPoints(predictors, n=1000, extf = 1.00)
##D 
##D # if desired, change working directory where subfolders of "models" and 
##D # "ensembles" will be created
##D # raster layers will be saved in subfolders of /models and /ensembles:
##D getwd()
##D 
##D # first calibrate the ensemble
##D # calibration is done in two steps
##D # in step 1, a k-fold procedure is used to determine the weights
##D # in step 2, models are calibrated for all presence and background locations
##D # factor is not used as it is not certain whether correct levels will be used
##D # it may therefore be better to use dummy variables
##D 
##D # step 1: determine weights through 4-fold cross-validation
##D ensemble.calibrate.step1 <- ensemble.calibrate.weights(
##D     x=predictors, p=pres, a=background, k=4, 
##D     SINK=TRUE, species.name="Bradypus",
##D     MAXENT=1, MAXLIKE=0, GBM=1, GBMSTEP=0, RF=0, GLM=1, GLMSTEP=0, 
##D     GAM=1, GAMSTEP=0, MGCV=1, MGCVFIX=0,
##D     EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, GLMNET=1,
##D     BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
##D     ENSEMBLE.tune=TRUE, PROBIT=TRUE,
##D     ENSEMBLE.best=0, ENSEMBLE.exponent=c(1, 2, 3),
##D     ENSEMBLE.min=c(0.65, 0.7),
##D     Yweights="BIOMOD",
##D     PLOTS=FALSE, formulae.defaults=TRUE)
##D 
##D # step 1 generated the weights for each algorithm
##D model.weights <- ensemble.calibrate.step1$output.weights
##D x.batch <- ensemble.calibrate.step1$x
##D p.batch <- ensemble.calibrate.step1$p
##D a.batch <- ensemble.calibrate.step1$a
##D MAXENT.a.batch <- ensemble.calibrate.step1$MAXENT.a
##D factors.batch <- ensemble.calibrate.step1$factors
##D dummy.vars.batch <- ensemble.calibrate.step1$dummy.vars
##D 
##D # step 2: calibrate models with all available presence locations
##D # weights determined in step 1 calculate ensemble in step 2
##D ensemble.calibrate.step2 <- ensemble.calibrate.models(
##D     x=x.batch, p=p.batch, a=a.batch, MAXENT.a=MAXENT.a.batch, 
##D     factors=factors.batch, dummy.vars=dummy.vars.batch, 
##D     SINK=TRUE, species.name="Bradypus",
##D     models.keep=TRUE,
##D     input.weights=model.weights,
##D     ENSEMBLE.tune=FALSE, PROBIT=TRUE,
##D     Yweights="BIOMOD",
##D     PLOTS=FALSE, formulae.defaults=TRUE)
##D 
##D # step 3: use previously calibrated models to create ensemble raster layers
##D # re-evaluate the created maps at presence and background locations
##D # (note that re-evaluation will be different due to truncation of raster layers
##D # as they wered saved as integer values ranged 0 to 1000)
##D ensemble.raster.results <- ensemble.raster(xn=predictors, 
##D     models.list=ensemble.calibrate.step2$models, 
##D     input.weights=model.weights,
##D     SINK=TRUE, evaluate=TRUE,
##D     RASTER.species.name="Bradypus", RASTER.stack.name="base")
##D 
##D # use the base map to check for changes in suitable habitat
##D # this type of analysis is typically done with different predictor layers
##D # (for example, predictor layers representing different possible future climates)
##D # In this example, changes from a previous model (ensemble.raster.results)
##D # are contrasted with a newly calibrated model (ensemble.raster.results2)
##D # step 1: 4-fold cross-validation
##D ensemble.calibrate2.step1 <- ensemble.calibrate.weights(
##D     x=x.batch, p=p.batch, a=a.batch, MAXENT.a=MAXENT.a.batch, 
##D     factors=factors.batch, dummy.vars=dummy.vars.batch, 
##D     k=4, 
##D     SINK=TRUE, species.name="Bradypus",
##D     MAXENT=1, MAXLIKE=0, GBM=1, GBMSTEP=0, RF=0, GLM=1, GLMSTEP=0, 
##D     GAM=1, GAMSTEP=0, MGCV=1, MGCVFIX=0,
##D     EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, GLMNET=1,
##D     BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
##D     ENSEMBLE.tune=TRUE, PROBIT=TRUE,
##D     ENSEMBLE.best=0, ENSEMBLE.exponent=c(1, 2, 3),
##D     ENSEMBLE.min=c(0.65, 0.7),
##D     Yweights="BIOMOD",
##D     PLOTS=FALSE, formulae.defaults=TRUE)
##D 
##D model.weights2 <- ensemble.calibrate2.step1$output.weights
##D 
##D ensemble.calibrate2.step2 <- ensemble.calibrate.models(
##D     x=x.batch, p=p.batch, a=a.batch, MAXENT.a=MAXENT.a.batch, 
##D     factors=factors.batch, dummy.vars=dummy.vars.batch, 
##D     SINK=TRUE, species.name="Bradypus",
##D     models.keep=TRUE,
##D     input.weights=model.weights2,
##D     ENSEMBLE.tune=FALSE, PROBIT=TRUE,
##D     Yweights="BIOMOD",
##D     PLOTS=FALSE, formulae.defaults=TRUE)
##D 
##D ensemble.raster.results2 <- ensemble.raster(
##D     xn=predictors, 
##D     models.list=ensemble.calibrate2.step2$models, 
##D     input.weights=model.weights2,
##D     SINK=TRUE, evaluate=TRUE,
##D     RASTER.species.name="Bradypus", RASTER.stack.name="recalibrated")
##D 
##D base.file <- paste(getwd(), "/ensembles/presence/Bradypus_base.grd", sep="")
##D other.file <- paste(getwd(), "/ensembles/presence/Bradypus_recalibrated.grd", sep="")
##D 
##D changed.habitat <- ensemble.habitat.change(base.map=base.file, 
##D     other.maps=c(other.file),
##D     change.folder="ensembles/change")
##D 
##D change.file <- paste(getwd(), "/ensembles/change/Bradypus_recalibrated_presence.grd", sep="")
##D 
##D par.old <- graphics::par(no.readonly=T)
##D dev.new()
##D par(mfrow=c(2,2))
##D raster::plot(raster(base.file), breaks=c(-1, 0, 1), col=c("grey", "green"), 
##D     legend.shrink=0.8, main="base presence")
##D raster::plot(raster(other.file), breaks=c(-1, 0, 1), col=c("grey", "green"), 
##D     legend.shrink=0.8, main="other presence")
##D raster::plot(raster(change.file), breaks=c(-1, 0, 1, 10, 11), 
##D     col=c("grey", "blue", "red", "green"), 
##D     legend.shrink=0.8, main="habitat change", sub="11 remaining, 10 lost, 1 new")
##D graphics::par(par.old)
##D 
##D areas <- ensemble.area(raster(change.file))
##D areas
## End(Not run)



cleanEx()
nameEx("ensemble.red")
### * ensemble.red

flush(stderr()); flush(stdout())

### Name: ensemble.red
### Title: Area of Occupancy (AOO) and Extent of Occurrence (EOO) via the
###   'red' library.
### Aliases: ensemble.red ensemble.chull.create ensemble.chull.apply

### ** Examples


## Not run: 
##D 
##D ## Not run: 
##D # based on examples in the dismo package
##D 
##D # get predictor variables
##D library(dismo)
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D # subset based on Variance Inflation Factors
##D predictors <- subset(predictors, subset=c("bio5", "bio6", 
##D     "bio16", "bio17"))
##D predictors
##D predictors@title <- "red"
##D 
##D # presence points
##D presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
##D pres <- read.table(presence_file, header=TRUE, sep=',')
##D 
##D # fit 5 ensemble models (could take some time!)
##D # (examples for the red package use 100 models)
##D ensembles <- ensemble.batch(x=predictors, 
##D     xn=c(predictors),
##D     species.presence=pres, 
##D     thin.km=100,
##D     k.splits=4, k.test=0, 
##D     n.ensembles=5, 
##D     SINK=TRUE, 
##D     ENSEMBLE.best=10, ENSEMBLE.exponent=c(1, 2, 3), 
##D     ENSEMBLE.min=0.6,
##D     MAXENT=0, MAXLIKE=1, GBM=1, GBMSTEP=0, RF=1, GLM=1, GLMSTEP=1, GAM=1, 
##D     GAMSTEP=1, MGCV=1, MGCVFIX=1, EARTH=1, RPART=1, NNET=1, FDA=1, 
##D     SVM=1, SVME=1, GLMNET=1,
##D     BIOCLIM.O=0, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
##D     PROBIT=TRUE,
##D     Yweights="BIOMOD",
##D     formulae.defaults=TRUE)
##D 
##D # first application of ensemble.red before applying the convex hull mask
##D # AOO and EOO are determined for each count level
##D library(red)
##D count.file <- paste(getwd(), "/ensembles/consensuscount/Bradypus variegatus_red.grd", sep="")
##D count.raster <- raster(count.file)
##D ensemble.red(count.raster)
##D 
##D # do not predict presence in polygons completely outside convex hull
##D # of known presence locations
##D pres.file <- paste(getwd(), "/ensembles/consensuspresence/Bradypus variegatus_red.grd", sep="")
##D pres.raster <- raster(pres.file)
##D pres1 <- pres[, -1]
##D chull.created <- ensemble.chull.create(x.pres=pres.raster, p=pres1)
##D 
##D mask.raster <- chull.created$mask.layer
##D mask.poly <- chull.created$convex.hull
##D par.old <- graphics::par(no.readonly=T)
##D par(mfrow=c(1,2))
##D plot(pres.raster, breaks=c(-1, 0, 1), col=c("grey", "green"),
##D     main="before convex hull")
##D points(pres1, col="blue")
##D 
##D pres.chull <- ensemble.chull.apply(pres.raster, mask=mask.raster, keep.old=T)
##D # load new
##D pres.file <- paste(getwd(), "/ensembles/consensuspresence/Bradypus variegatus_red.grd", sep="")
##D pres.raster <- raster(pres.file)
##D plot(pres.raster, breaks=c(-1, 0, 1), col=c("grey", "green"),
##D     main="after convex hull")
##D plot(mask.poly, add=T, border="blue")
##D 
##D # new application of ensemble.red
##D dev.new()
##D plot(count.raster, main="before convex hull")
##D ensemble.red(count.raster)
##D # all cells where species is predicted not to be present according to the mask layer
##D # will be modified to a count of zero
##D count.chull <- ensemble.chull.apply(count.raster, mask=mask.raster, keep.old=T)
##D # load new
##D count.file <- paste(getwd(), "/ensembles/consensuscount/Bradypus variegatus_red.grd", sep="")
##D count.raster <- raster(count.file)
##D ensemble.red(count.raster)
##D dev.new()
##D plot(count.raster, main="after convex hull")
##D par.old <- graphics::par(no.readonly=T)
##D 
## End(Not run)



cleanEx()
nameEx("ensemble.spatialThin")
### * ensemble.spatialThin

flush(stderr()); flush(stdout())

### Name: ensemble.spatialThin
### Title: Spatial thinning of presence point locations using the highly
###   accurate geodesic estimates from the geosphere package
### Aliases: ensemble.spatialThin

### ** Examples

## Not run: 
##D # get predictor variables, only needed for plotting
##D library(dismo)
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D # subset based on Variance Inflation Factors
##D predictors <- subset(predictors, subset=c("bio5", "bio6", 
##D     "bio16", "bio17", "biome"))
##D predictors
##D predictors@title <- "base"
##D 
##D # presence points
##D presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
##D pres <- read.table(presence_file, header=TRUE, sep=',')[, -1]
##D 
##D # number of locations
##D nrow(pres)
##D 
##D par.old <- graphics::par(no.readonly=T)
##D par(mfrow=c(2,2))
##D 
##D pres.thin1 <- ensemble.spatialThin(pres, thin.km=100, runs=10, verbose=T)
##D plot(predictors[[1]], main="10 runs", ext=extent(SpatialPoints(pres.thin1)))
##D points(pres.thin1, pch=20, col="red")
##D 
##D pres.thin2 <- ensemble.spatialThin(pres, thin.km=100, runs=10, verbose=T)
##D plot(predictors[[1]], main="10 runs", ext=extent(SpatialPoints(pres.thin2)))
##D points(pres.thin2, pch=20, col="red")
##D 
##D pres.thin3 <- ensemble.spatialThin(pres, thin.km=100, runs=100, verbose=T)
##D plot(predictors[[1]], main="100 runs", ext=extent(SpatialPoints(pres.thin3)))
##D points(pres.thin3, pch=20, col="red")
##D 
##D pres.thin4 <- ensemble.spatialThin(pres, thin.km=100, runs=100, verbose=T)
##D plot(predictors[[1]], main="100 runs", ext=extent(SpatialPoints(pres.thin4)))
##D points(pres.thin4, pch=20, col="red")
##D 
##D graphics::par(par.old)
##D 
## End(Not run)




cleanEx()
nameEx("ensemble.zones")
### * ensemble.zones

flush(stderr()); flush(stdout())

### Name: ensemble.zones
### Title: Mapping of environmental zones based on the Mahalanobis distance
###   from centroids in environmental space.
### Aliases: ensemble.zones ensemble.centroids

### ** Examples


## Not run: 
##D # get predictor variables
##D library(dismo)
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D predictors <- subset(predictors, subset=c("bio1", "bio5", "bio6", "bio7", "bio8", 
##D     "bio12", "bio16", "bio17"))
##D predictors
##D predictors@title <- "base"
##D 
##D # choose background points
##D background <- randomPoints(predictors, n=1000, extf=1.00)
##D 
##D # predicted presence from GLM
##D ensemble.calibrate.step1 <- ensemble.calibrate.models(
##D     x=predictors, p=pres, a=background,
##D     species.name="Bradypus",
##D     MAXENT=0, MAXLIKE=0, GBM=0, GBMSTEP=0, RF=0, GLM=1, GLMSTEP=0, 
##D     GAM=0, GAMSTEP=0, MGCV=0, MGCVFIX=0,
##D     EARTH=0, RPART=0, NNET=0, FDA=0, SVM=0, SVME=0, GLMNET=0,
##D     BIOCLIM.O=0, BIOCLIM=0, DOMAIN=0, MAHAL=0, MAHAL01=0,
##D     Yweights="BIOMOD",
##D     models.keep=TRUE)
##D 
##D ensemble.raster.results <- ensemble.raster(xn=predictors, 
##D     models.list=ensemble.calibrate.step1$models, 
##D     RASTER.species.name="Bradypus", RASTER.stack.name="base")
##D 
##D # get presence map as for example created with ensemble.raster in subfolder 'ensemble/presence'
##D # presence values are values equal to 1
##D presence.file <- paste("ensembles//presence//Bradypus_base.grd", sep="")
##D presence.raster <- raster(presence.file)
##D 
##D # let cascadeKM decide on the number of clusters
##D dev.new()
##D centroids <- ensemble.centroids(presence.raster=presence.raster, 
##D     x=predictors, an=1000, plotit=T)
##D ensemble.zones(presence.raster=presence.raster, centroid.object=centroids, 
##D     x=predictors, RASTER.species.name="Bradypus", KML.out=T)
##D 
##D dev.new()
##D zones.file <- paste("ensembles//zones//Bradypus_base.grd", sep="")
##D zones.raster <- raster(zones.file)
##D max.zones <- maxValue(zones.raster)
##D plot(zones.raster, breaks=c(0, c(1:max.zones)), 
##D     col = grDevices::rainbow(n=max.zones), main="zones")
##D ensemble.zones(presence.raster=presence.raster, centroid.object=centroids, 
##D     x=predictors, RASTER.species.name="Bradypus", KML.out=T)
##D 
##D # manually choose 6 zones
##D dev.new()
##D centroids6 <- ensemble.centroids(presence.raster=presence.raster, 
##D     x=predictors, an=1000, plotit=T, centers=6)
##D ensemble.zones(presence.raster=presence.raster, centroid.object=centroids6, 
##D     x=predictors, RASTER.species.name="Bradypus6", KML.out=T)
##D 
##D dev.new()
##D zones.file <- paste("ensembles//zones//Bradypus6_base.grd", sep="")
##D zones.raster <- raster(zones.file)
##D max.zones <- maxValue(zones.raster)
##D plot(zones.raster, breaks=c(0, c(1:max.zones)), 
##D     col = grDevices::rainbow(n=max.zones), main="six zones")
##D 
## End(Not run)



cleanEx()
nameEx("evaluation.strip")
### * evaluation.strip

flush(stderr()); flush(stdout())

### Name: evaluation.strip.data
### Title: Evaluation strips for ensemble suitability mapping
### Aliases: evaluation.strip.data evaluation.strip.plot

### ** Examples

## Not run: 
##D 
##D # get predictor variables
##D library(dismo)
##D predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
##D     pattern='grd', full.names=TRUE)
##D predictors <- stack(predictor.files)
##D # subset based on Variance Inflation Factors
##D predictors <- subset(predictors, subset=c("bio5", "bio6", 
##D     "bio16", "bio17"))
##D predictors <- stack(predictors)
##D predictors
##D predictors@title <- "base"
##D 
##D # presence points
##D presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
##D pres <- read.table(presence_file, header=TRUE, sep=',')[,-1]
##D 
##D # the kfold function randomly assigns data to groups; 
##D # groups are used as calibration (1/5) and training (4/5) data
##D groupp <- kfold(pres, 5)
##D pres_train <- pres[groupp !=  1, ]
##D pres_test <- pres[groupp ==  1, ]
##D 
##D # choose background points
##D background <- randomPoints(predictors, n=1000, extf=1.00)
##D colnames(background)=c('lon', 'lat')
##D groupa <- kfold(background, 5)
##D backg_train <- background[groupa != 1, ]
##D backg_test <- background[groupa == 1, ]
##D 
##D # calibrate the models
##D # MAXLIKE not included as does not allow predictions for data.frames
##D # ENSEMBLE.min and ENSEMBLE.weight.min set very low to explore all
##D # algorithms.
##D # If focus is on actual ensemble, then set ENSEMBLE.min and 
##D # ENSEMBLE.weight.min to more usual values
##D ensemble.calibrate <- ensemble.calibrate.models(x=predictors, 
##D     p=pres_train, a=backg_train, 
##D     pt=pres_test, at=backg_test,
##D     ENSEMBLE.min=0.5, ENSEMBLE.weight.min = 0.001,
##D     MAXENT=1, MAXLIKE=0, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, 
##D     GAMSTEP=1, MGCV=1, MGCVFIX=1, EARTH=1, RPART=1, NNET=1, FDA=1, 
##D     SVM=1, SVME=1, GLMNET=1,
##D     BIOCLIM.O=1, BIOCLIM=1, DOMAIN=1, MAHAL=0, MAHAL01=1,
##D     Yweights="BIOMOD", 
##D     PLOTS=FALSE, models.keep=TRUE)
##D 
##D # obtain data for plotting the evaluation strip
##D strip.data <- evaluation.strip.data(xn=predictors, steps=500,
##D     models.list=ensemble.calibrate$models)
##D 
##D # in case predictions for DOMAIN failed
##D # however, ENSEMBLE should also be recalculated
##D DOMAIN.model <- ensemble.calibrate$models$DOMAIN
##D strip.data$plot.data[, "DOMAIN"] <- dismo::predict(object=DOMAIN.model, 
##D     x=strip.data$plot.data)
##D 
##D # in case predictions for MAHAL01 failed
##D predict.MAHAL01 <- function(model, newdata, MAHAL.shape) {
##D     p <- dismo::predict(object=model, x=newdata)
##D     p <- p - 1 - MAHAL.shape
##D     p <- abs(p)
##D     p <- MAHAL.shape / p
##D     return(as.numeric(p))
##D }
##D 
##D MAHAL01.model <- ensemble.calibrate$models$MAHAL01
##D MAHAL.shape1 <- ensemble.calibrate$models$formulae$MAHAL.shape
##D strip.data$plot.data[, "MAHAL01"] <- predict.MAHAL01(model=MAHAL01.model, 
##D     newdata=strip.data$plot.data, MAHAL.shape=MAHAL.shape1)
##D 
##D # create graphs
##D evaluation.strip.plot(data=strip.data$plot.data, variable.focal="bio6",
##D     TrainData=strip.data$TrainData,
##D     type="o", col="red")
##D evaluation.strip.plot(data=strip.data$plot.data, model.focal="ENSEMBLE",
##D     TrainData=strip.data$TrainData,
##D     type="o", col="red")
##D 
## End(Not run)



cleanEx()
nameEx("faramea")
### * faramea

flush(stderr()); flush(stdout())

### Name: faramea
### Title: Faramea occidentalis abundance in Panama
### Aliases: faramea
### Keywords: datasets

### ** Examples

data(faramea)



cleanEx()
nameEx("ifri")
### * ifri

flush(stderr()); flush(stdout())

### Name: ifri
### Title: Example data from the International Forestry Resources and
###   Institutions (IFRI) research network
### Aliases: ifri
### Keywords: datasets

### ** Examples

data(ifri)



cleanEx()
nameEx("importancevalue")
### * importancevalue

flush(stderr()); flush(stdout())

### Name: importancevalue
### Title: Importance Value
### Aliases: importancevalue importancevalue.comp
### Keywords: multivariate

### ** Examples

data(ifri)
importancevalue(ifri, site='plotID', species='species', count='count', 
    basal='basal', factor='forest', level='YSF')
importancevalue.comp(ifri, site='plotID', species='species', count='count', 
    basal='basal', factor='forest')

# When all survey plots are the same size, importance value
# is not affected. Counts and basal areas now calculated per square metre
ifri$count <- ifri$count/314.16
ifri$basal <- ifri$basal/314.16

importancevalue(ifri, site='plotID', species='species', count='count', 
    basal='basal', factor='forest', level='YSF')
importancevalue.comp(ifri, site='plotID', species='species', count='count', 
    basal='basal', factor='forest')




cleanEx()
nameEx("makecommunitydataset")
### * makecommunitydataset

flush(stderr()); flush(stdout())

### Name: makecommunitydataset
### Title: Make a Community Dataset from a Stacked Dataset
### Aliases: makecommunitydataset stackcommunitydataset
### Keywords: multivariate

### ** Examples


## Not run: 
##D dune.file <- normalizePath(paste(system.file(package="BiodiversityR"), 
##D     '/etc/dunestacked.csv', sep=''))
##D dune.stacked <- read.csv(dune.file)
##D 
##D # dune.stacked has different variables for sites, species and abundance
##D head(dune.stacked)
##D dune.comm2 <- makecommunitydataset(dune.stacked, row='sites', column='species', 
##D     value='abundance')
##D 
##D # recreate the original stack
##D dune.stacked2 <- stackcommunitydataset(dune.comm2, remove.zeroes=T)
##D 
## End(Not run)




cleanEx()
nameEx("multiconstrained")
### * multiconstrained

flush(stderr()); flush(stdout())

### Name: multiconstrained
### Title: Pairwise Comparisons for All Levels of a Categorical Variable by
###   RDA, CCA or Capscale
### Aliases: multiconstrained
### Keywords: multivariate

### ** Examples

## Not run: 
##D library(vegan)
##D library(MASS)
##D data(dune)
##D data(dune.env)
##D multiconstrained(method="capscale", dune~Management, data=dune.env,
##D     distance="bray",add=TRUE)
##D multiconstrained(method="capscale", dune~Management, data=dune.env, 
##D     distance="bray", add=TRUE, contrast=3)
## End(Not run)



cleanEx()
nameEx("nested.anova.dbrda")
### * nested.anova.dbrda

flush(stderr()); flush(stdout())

### Name: nested.anova.dbrda
### Title: Nested Analysis of Variance via Distance-based Redundancy
###   Analysis or Non-parametric Multivariate Analysis of Variance
### Aliases: nested.anova.dbrda nested.npmanova
### Keywords: multivariate

### ** Examples

## Not run: 
##D library(vegan)
##D data(warcom)
##D data(warenv)
##D # use larger number of permutations for real studies
##D nested.npmanova(warcom~rift.valley+popshort, data=warenv, method="jac", 
##D     permutations=5)
##D nested.anova.dbrda(warcom~rift.valley+popshort, data=warenv, method="jac", 
##D     permutations=5)
## End(Not run)



cleanEx()
nameEx("nnetrandom")
### * nnetrandom

flush(stderr()); flush(stdout())

### Name: nnetrandom
### Title: Calculate the NNET Result with the Smallest Value from Various
###   Random Starts
### Aliases: nnetrandom
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(faramea)
##D faramea <- na.omit(faramea)
##D faramea$presence <- as.numeric(faramea$Faramea.occidentalis > 0)
##D attach(faramea)
##D library(nnet)
##D result <- nnetrandom(presence ~ Elevation, data=faramea, size=2, 
##D     skip=FALSE, entropy=TRUE, trace=FALSE, maxit=1000, tries=100, 
##D     leave.one.out=FALSE)
##D summary(result)
##D result$fitted.values
##D result$value
##D result2 <- nnetrandom(presence ~ Elevation, data=faramea, size=2, 
##D     skip=FALSE, entropy=TRUE, trace=FALSE, maxit=1000, tries=50, 
##D     leave.one.out=TRUE)
##D result2$range
##D result2$CV
##D result2$successful
## End(Not run)



cleanEx()
nameEx("ordicoeno")
### * ordicoeno

flush(stderr()); flush(stdout())

### Name: ordicoeno
### Title: Coenoclines for an Ordination Axis
### Aliases: ordicoeno
### Keywords: multivariate

### ** Examples

library(vegan)
library(mgcv)
data(dune)
Ordination.model1 <- rda(dune)
plot1 <- ordiplot(Ordination.model1, choices=c(1,2), scaling=1)
ordicoeno(dune, ordiplot=plot1, legend=TRUE)



cleanEx()
nameEx("ordisymbol")
### * ordisymbol

flush(stderr()); flush(stdout())

### Name: ordisymbol
### Title: Add Other Graphical Items to Ordination Diagrams
### Aliases: ordisymbol ordibubble ordicluster2 ordinearest ordivector
### Keywords: multivariate

### ** Examples

library(vegan)
data(dune)
data(dune.env)
Ordination.model1 <- rda(dune)
plot1 <- ordiplot(Ordination.model1, choices=c(1,2), scaling=2)
ordisymbol(plot1, dune.env, "Management", legend=TRUE, 
    legend.x="topleft", legend.ncol=1)
plot2 <- ordiplot(Ordination.model1, choices=c(1,2), scaling=1)
distmatrix <- vegdist(dune, method='bray')
cluster <- hclust(distmatrix, method='single')
ordicluster2(plot2, cluster)
ordinearest(plot2, distmatrix, col=2)
ordivector(plot2, "Agrostol", lty=2)



cleanEx()
nameEx("radfitresult")
### * radfitresult

flush(stderr()); flush(stdout())

### Name: radfitresult
### Title: Alternative Rank Abundance Fitting Results
### Aliases: radfitresult
### Keywords: multivariate

### ** Examples

library(vegan)
data(BCI)
BCIall <- t(as.matrix(colSums(BCI)))
radfitresult(BCIall)



cleanEx()
nameEx("rankabundance")
### * rankabundance

flush(stderr()); flush(stdout())

### Name: rankabundance
### Title: Rank Abundance Curves
### Aliases: rankabundance rankabunplot rankabuncomp
### Keywords: multivariate

### ** Examples

library(vegan)
data(dune.env)
data(dune)
RankAbun.1 <- rankabundance(dune)
RankAbun.1
rankabunplot(RankAbun.1, scale='abundance', addit=FALSE, specnames=c(1,2,3))
rankabunplot(RankAbun.1, scale='logabun', addit=FALSE, specnames=c(1:30), 
    srt=45, ylim=c(1,100))
rankabuncomp(dune, y=dune.env, factor='Management', 
    scale='proportion', legend=FALSE)
## CLICK IN THE GRAPH TO INDICATE WHERE THE LEGEND NEEDS TO BE PLACED
## IF YOU OPT FOR LEGEND=TRUE.



cleanEx()
nameEx("removeNAcomm")
### * removeNAcomm

flush(stderr()); flush(stdout())

### Name: removeNAcomm
### Title: Synchronize Community and Environmental Datasets
### Aliases: removeNAcomm removeNAenv same.sites check.datasets
###   check.ordiscores replaceNAcomm removezerospecies subsetcomm
### Keywords: multivariate

### ** Examples

library(vegan)
data(dune.env)
data(dune)
dune.env2 <- dune.env
dune.env2[1:4,"Moisture"] <- NA
dune2 <- removeNAcomm(dune,dune.env2,"Moisture")
dune.env2 <- removeNAenv(dune.env2,"Moisture")
dune3 <- same.sites(dune,dune.env2)
check.datasets(dune,dune.env2)
check.datasets(dune2,dune.env2)
check.datasets(dune3,dune.env2)
dune4 <- subsetcomm(dune,dune.env,"Management","NM",returncomm=TRUE)
dune.env4 <- subsetcomm(dune,dune.env,"Management","NM",returncomm=FALSE)
dune5 <- same.sites(dune,dune.env4)
check.datasets(dune4,dune5)



cleanEx()
nameEx("renyiresult")
### * renyiresult

flush(stderr()); flush(stdout())

### Name: renyiresult
### Title: Alternative Renyi Diversity Results
### Aliases: renyiresult renyiplot renyiaccumresult renyicomp
### Keywords: multivariate

### ** Examples

library(vegan)
data(dune.env)
data(dune)
Renyi.1 <- renyiresult(dune, y=dune.env, factor='Management', level='NM', 
    method='s')
Renyi.1
renyiplot(Renyi.1, evenness=FALSE, addit=FALSE, pch=1,col='1', cex=1, 
    legend=FALSE)
## CLICK IN THE GRAPH TO INDICATE WHERE THE LEGEND NEEDS TO BE PLACED
## IN CASE THAT YOU OPT FOR LEGEND=TRUE



cleanEx()
nameEx("residualssurface")
### * residualssurface

flush(stderr()); flush(stdout())

### Name: residualssurface
### Title: Show and Interpolate Two Dimensional Distribution of Residuals
### Aliases: residualssurface
### Keywords: multivariate

### ** Examples

library(vegan)
library(mgcv)
library(akima)
data(faramea)
Count.model1 <- lm(Faramea.occidentalis ~ Precipitation,
    data=faramea, na.action=na.exclude)
surface.1 <- residualssurface(Count.model1, na.omit(faramea),
    'UTM.EW', 'UTM.NS', gam=TRUE, plotit=TRUE, bubble=TRUE)



cleanEx()
nameEx("spatialsample")
### * spatialsample

flush(stderr()); flush(stdout())

### Name: spatialsample
### Title: Spatial Sampling within a Polygon
### Aliases: spatialsample
### Keywords: multivariate

### ** Examples

library(splancs)
area <- array(c(10,10,15,35,40,35,5,35,35,30,30,10), dim=c(6,2))
landuse1 <- array(c(10,10,15,15,30,35,35,30), dim=c(4,2))
landuse2 <- array(c(10,10,15,15,35,30,10,30,30,35,30,15), dim=c(6,2))
landuse3 <- array(c(10,10,30,35,40,35,5,10,15,30,30,10), dim=c(6,2))
plot(area[,1], area[,2], type="n", xlab="horizontal position", 
    ylab="vertical position", lwd=2, bty="l")
polygon(landuse1)
polygon(landuse2)
polygon(landuse3)
spatialsample(area, method="random", n=20, xwidth=1, ywidth=1, plotit=TRUE, 
    plothull=FALSE)
spatialsample(area, method="grid", xwidth=1, ywidth=1, plotit=TRUE, xleft=12, 
    ylower=7, xdist=4, ydist=4)
spatialsample(area, method="random grid", n=20, xwidth=1, ywidth=1, 
    plotit=TRUE, xleft=12, ylower=7, xdist=4, ydist=4)



cleanEx()
nameEx("transfgradient")
### * transfgradient

flush(stderr()); flush(stdout())

### Name: transfgradient
### Title: Gradient for Hypothetical Example of Turover of Species
###   Composition
### Aliases: transfgradient
### Keywords: datasets

### ** Examples

data(transfspecies)
data(transfgradient)
plot(transfspecies[,1]~transfgradient[,1],xlab="gradient",
    ylab="species abundance",type="n",ylim=c(0.5,8.5))
for (i in 1:9) {points(transfgradient[,1],transfspecies[,i],type="o",pch=i)}



cleanEx()
nameEx("transfspecies")
### * transfspecies

flush(stderr()); flush(stdout())

### Name: transfspecies
### Title: Hypothetical Example of Turover of Species Composition
### Aliases: transfspecies
### Keywords: datasets

### ** Examples

data(transfspecies)
data(transfgradient)
plot(transfspecies[,1]~transfgradient[,1],xlab="gradient",
    ylab="species abundance",type="n",ylim=c(0.5,8.5))
for (i in 1:9) {points(transfgradient[,1],transfspecies[,i],type="o",pch=i)}



cleanEx()
nameEx("warcom")
### * warcom

flush(stderr()); flush(stdout())

### Name: warcom
### Title: Warburgia ugandensis AFLP Scores
### Aliases: warcom
### Keywords: datasets

### ** Examples

data(warcom)



cleanEx()
nameEx("warenv")
### * warenv

flush(stderr()); flush(stdout())

### Name: warenv
### Title: Warburgia ugandensis Population Structure
### Aliases: warenv
### Keywords: datasets

### ** Examples

data(warenv)



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
