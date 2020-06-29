pkgname <- "ecospat"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ecospat')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ecospat.CCV.communityEvaluation.bin")
### * ecospat.CCV.communityEvaluation.bin

flush(stderr()); flush(stdout())

### Name: ecospat.CCV.communityEvaluation.bin
### Title: Calculates a range of community evaluation metrics based on
###   different thresholding techniques.
### Aliases: ecospat.CCV.communityEvaluation.bin

### ** Examples




cleanEx()
nameEx("ecospat.CCV.communityEvaluation.prob")
### * ecospat.CCV.communityEvaluation.prob

flush(stderr()); flush(stdout())

### Name: ecospat.CCV.communityEvaluation.prob
### Title: Evaluates community predictions directly on the probabilities
###   (i.e., threshold independent)
### Aliases: ecospat.CCV.communityEvaluation.prob

### ** Examples




cleanEx()
nameEx("ecospat.CCV.createDataSplitTable")
### * ecospat.CCV.createDataSplitTable

flush(stderr()); flush(stdout())

### Name: ecospat.CCV.createDataSplitTable
### Title: Creates a DataSplitTable for usage in ecospat.ccv.modeling.
### Aliases: ecospat.CCV.createDataSplitTable

### ** Examples





cleanEx()
nameEx("ecospat.CCV.modeling")
### * ecospat.CCV.modeling

flush(stderr()); flush(stdout())

### Name: ecospat.CCV.modeling
### Title: Runs indivudual species distribuion models with SDMs or ESMs
### Aliases: ecospat.CCV.modeling

### ** Examples




cleanEx()
nameEx("ecospat.CommunityEval")
### * ecospat.CommunityEval

flush(stderr()); flush(stdout())

### Name: ecospat.CommunityEval
### Title: Community Evaluation
### Aliases: ecospat.CommunityEval

### ** Examples



cleanEx()
nameEx("ecospat.Cscore")
### * ecospat.Cscore

flush(stderr()); flush(stdout())

### Name: ecospat.Cscore
### Title: Pairwise co-occurrence Analysis with calculation of the C-score
###   index.
### Aliases: ecospat.Cscore

### ** Examples

## Not run: 
##D data<- ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
##D nperm <- 10000
##D outpath <- getwd()
##D Cscore<-ecospat.Cscore(data, nperm, outpath)
##D 
## End(Not run)


cleanEx()
nameEx("ecospat.ESM.EnsembleModeling")
### * ecospat.ESM.EnsembleModeling

flush(stderr()); flush(stdout())

### Name: ecospat.ESM.EnsembleModeling
### Title: Ensamble of Small Models: Evaluates and Averages Simple
###   Bivariate Models To ESMs
### Aliases: ecospat.ESM.EnsembleModeling

### ** Examples



cleanEx()
nameEx("ecospat.ESM.EnsembleProjection")
### * ecospat.ESM.EnsembleProjection

flush(stderr()); flush(stdout())

### Name: ecospat.ESM.EnsembleProjection
### Title: Ensamble of Small Models: Projects Calibrated ESMs Into New
###   Space Or Time.
### Aliases: ecospat.ESM.EnsembleProjection

### ** Examples



cleanEx()
nameEx("ecospat.ESM.Modeling")
### * ecospat.ESM.Modeling

flush(stderr()); flush(stdout())

### Name: ecospat.ESM.Modeling
### Title: Ensamble of Small Models: Calibration of Simple Bivariate Models
### Aliases: ecospat.ESM.Modeling

### ** Examples



cleanEx()
nameEx("ecospat.ESM.Projection")
### * ecospat.ESM.Projection

flush(stderr()); flush(stdout())

### Name: ecospat.ESM.Projection
### Title: Ensamble of Small Models: Projects Simple Bivariate Models Into
###   New Space Or Time
### Aliases: ecospat.ESM.Projection

### ** Examples



cleanEx()
nameEx("ecospat.ESM.VarContrib")
### * ecospat.ESM.VarContrib

flush(stderr()); flush(stdout())

### Name: ecospat.ESM.VarContrib
### Title: Variable contribution in ESM
### Aliases: ecospat.ESM.VarContrib

### ** Examples




cleanEx()
nameEx("ecospat.ESM.responsePlot")
### * ecospat.ESM.responsePlot

flush(stderr()); flush(stdout())

### Name: ecospat.ESM.responsePlot
### Title: Produce response plots for ESMs
### Aliases: ecospat.ESM.responsePlot

### ** Examples



cleanEx()
nameEx("ecospat.ESM.threshold")
### * ecospat.ESM.threshold

flush(stderr()); flush(stdout())

### Name: ecospat.ESM.threshold
### Title: Thresholds for Ensamble of Small Models
### Aliases: ecospat.ESM.threshold

### ** Examples



cleanEx()
nameEx("ecospat.Epred")
### * ecospat.Epred

flush(stderr()); flush(stdout())

### Name: ecospat.Epred
### Title: Prediction Mean
### Aliases: ecospat.Epred

### ** Examples

x <- ecospat.testData[c(92,96)]
mean <- ecospat.Epred (x, w=rep(1,ncol(x)), th=0.5)



cleanEx()
nameEx("ecospat.SESAM.prr")
### * ecospat.SESAM.prr

flush(stderr()); flush(stdout())

### Name: ecospat.SESAM.prr
### Title: SESAM Probability Ranking Rule
### Aliases: ecospat.SESAM.prr

### ** Examples

proba <- ecospat.testData[,73:92]
sr <- as.data.frame(rowSums(proba))
ppr<-ecospat.SESAM.prr(proba, sr)
head(ppr)




cleanEx()
nameEx("ecospat.adj.D2")
### * ecospat.adj.D2

flush(stderr()); flush(stdout())

### Name: ecospat.adj.D2.glm
### Title: Calculate An Adjusted D2
### Aliases: ecospat.adj.D2.glm

### ** Examples


data(ecospat.testData)
glm.obj<-glm(Achillea_millefolium~ddeg+mind+srad+slp+topo, 
family = binomial, data=ecospat.testData)

ecospat.adj.D2.glm(glm.obj)




cleanEx()
nameEx("ecospat.binary.model")
### * ecospat.binary.model

flush(stderr()); flush(stdout())

### Name: ecospat.binary.model
### Title: Generate Binary Models
### Aliases: ecospat.binary.model

### ** Examples





cleanEx()
nameEx("ecospat.boyce")
### * ecospat.boyce

flush(stderr()); flush(stdout())

### Name: ecospat.boyce
### Title: Calculate Boyce Index
### Aliases: ecospat.boyce

### ** Examples

obs <- (ecospat.testData$glm_Saxifraga_oppositifolia
[which(ecospat.testData$Saxifraga_oppositifolia==1)])

ecospat.boyce (fit = ecospat.testData$glm_Saxifraga_oppositifolia , obs, nclass=0, 
window.w="default", res=100, PEplot = TRUE)



cleanEx()
nameEx("ecospat.calculate.pd")
### * ecospat.calculate.pd

flush(stderr()); flush(stdout())

### Name: ecospat.calculate.pd
### Title: Calculate Phylogenetic Diversity Measures
### Aliases: ecospat.calculate.pd

### ** Examples

fpath <- system.file("extdata", "ecospat.testTree.tre", package="ecospat")
tree <-read.tree(fpath)
data <- ecospat.testData[9:52] 

pd <- ecospat.calculate.pd(tree, data, method = "spanning", type = "species", root = FALSE, 
average = FALSE, verbose = TRUE )

plot(pd)



cleanEx()
nameEx("ecospat.caleval")
### * ecospat.caleval

flush(stderr()); flush(stdout())

### Name: ecospat.caleval
### Title: Calibration And Evaluation Dataset
### Aliases: ecospat.caleval

### ** Examples

data <- ecospat.testData
caleval <- ecospat.caleval (data = ecospat.testData[53], xy = data[2:3], row.num = 1:nrow(data), 
nrep = 2, ratio = 0.7, disaggregate = 0.2, pseudoabs = 100, npres = 10, replace = FALSE)
caleval



cleanEx()
nameEx("ecospat.climan")
### * ecospat.climan

flush(stderr()); flush(stdout())

### Name: ecospat.climan
### Title: A climate analogy setection tool for the modeling of species
###   distributions
### Aliases: ecospat.climan

### ** Examples

x <- ecospat.testData[c(4:8)]
p<- x[1:90,] #A projection dataset.
ref<- x[91:300,] #A reference dataset
ecospat.climan(ref,p)




cleanEx()
nameEx("ecospat.co_occurrences")
### * ecospat.co_occurrences

flush(stderr()); flush(stdout())

### Name: ecospat.co_occurrences
### Title: Species Co-Occurrences
### Aliases: ecospat.co_occurrences

### ** Examples



cleanEx()
nameEx("ecospat.cohen.kappa")
### * ecospat.cohen.kappa

flush(stderr()); flush(stdout())

### Name: ecospat.cohen.kappa
### Title: Cohen's Kappa
### Aliases: ecospat.cohen.kappa

### ** Examples

Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris
th <- 0.39 # threshold
xtab <- table(Pred >= th, Sp.occ)

ecospat.cohen.kappa(xtab)



cleanEx()
nameEx("ecospat.cons_Cscore")
### * ecospat.cons_Cscore

flush(stderr()); flush(stdout())

### Name: ecospat.cons_Cscore
### Title: Constrained Co-Occurrence Analysis.
### Aliases: ecospat.cons_Cscore

### ** Examples



cleanEx()
nameEx("ecospat.cor.plot")
### * ecospat.cor.plot

flush(stderr()); flush(stdout())

### Name: ecospat.cor.plot
### Title: Correlation Plot
### Aliases: ecospat.cor.plot

### ** Examples

data <- ecospat.testData[,4:8]
ecospat.cor.plot(data)



cleanEx()
nameEx("ecospat.cv.gbm")
### * ecospat.cv.gbm

flush(stderr()); flush(stdout())

### Name: ecospat.cv.gbm
### Title: GBM Cross Validation
### Aliases: ecospat.cv.gbm

### ** Examples

data('ecospat.testData')

# data for Soldanella alpina
data.Solalp<- ecospat.testData[c("Soldanella_alpina","ddeg","mind","srad","slp","topo")] 

# gbm model for Soldanella alpina
gbm.Solalp <- gbm(Soldanella_alpina ~ ., data = data.Solalp,
                  distribution = "bernoulli", cv.folds = 10, n.cores=2)

# cross-validated predictions
gbm.pred <- ecospat.cv.gbm (gbm.obj= gbm.Solalp,data.Solalp, 
                            K=10, cv.lim=10, jack.knife=FALSE)



cleanEx()
nameEx("ecospat.cv.glm")
### * ecospat.cv.glm

flush(stderr()); flush(stdout())

### Name: ecospat.cv.glm
### Title: GLM Cross Validation
### Aliases: ecospat.cv.glm

### ** Examples




cleanEx()
nameEx("ecospat.cv.me")
### * ecospat.cv.me

flush(stderr()); flush(stdout())

### Name: ecospat.cv.me
### Title: Maxent Cross Validation
### Aliases: ecospat.cv.me

### ** Examples




cleanEx()
nameEx("ecospat.cv.rf")
### * ecospat.cv.rf

flush(stderr()); flush(stdout())

### Name: ecospat.cv.rf
### Title: RandomForest Cross Validation
### Aliases: ecospat.cv.rf

### ** Examples

data('ecospat.testData')

# data for Soldanella alpina
data.Solalp<- ecospat.testData[c("Soldanella_alpina","ddeg","mind","srad","slp","topo")] 

library(randomForest)
rf.Solalp <- randomForest(x = data.Solalp[,-1], y = as.factor(data.Solalp[,1]))
rf.pred <- ecospat.cv.rf(rf.Solalp, data.Solalp, K = 10, cv.lim = 10, 
                         jack.knife = FALSE, verbose = FALSE)



cleanEx()
nameEx("ecospat.grid.clim.dyn")
### * ecospat.grid.clim.dyn

flush(stderr()); flush(stdout())

### Name: ecospat.grid.clim.dyn
### Title: Dynamic Occurrence Densities Grid
### Aliases: ecospat.grid.clim.dyn

### ** Examples



cleanEx()
nameEx("ecospat.makeDataFrame")
### * ecospat.makeDataFrame

flush(stderr()); flush(stdout())

### Name: ecospat.makeDataFrame
### Title: Make Data Frame
### Aliases: ecospat.makeDataFrame

### ** Examples




cleanEx()
nameEx("ecospat.mantel.correlogram")
### * ecospat.mantel.correlogram

flush(stderr()); flush(stdout())

### Name: ecospat.mantel.correlogram
### Title: Mantel Correlogram
### Aliases: ecospat.mantel.correlogram

### ** Examples

ecospat.mantel.correlogram(dfvar=ecospat.testData[c(2:16)],colxy=1:2, n=100, colvar=3:7, 
max=1000, nclass=10, nperm=100)



cleanEx()
nameEx("ecospat.max.kappa")
### * ecospat.max.kappa

flush(stderr()); flush(stdout())

### Name: ecospat.max.kappa
### Title: Maximum Kappa
### Aliases: ecospat.max.kappa

### ** Examples






cleanEx()
nameEx("ecospat.max.tss")
### * ecospat.max.tss

flush(stderr()); flush(stdout())

### Name: ecospat.max.tss
### Title: Maximum TSS
### Aliases: ecospat.max.tss
### Keywords: file

### ** Examples


data(ecospat.testData)
Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris
TSS100 <- ecospat.max.tss(Pred, Sp.occ)



cleanEx()
nameEx("ecospat.maxentvarimport")
### * ecospat.maxentvarimport

flush(stderr()); flush(stdout())

### Name: ecospat.maxentvarimport
### Title: Maxent Variable Importance
### Aliases: ecospat.maxentvarimport

### ** Examples

library(dismo)
data('ecospat.testData')

# data for Soldanella alpina
data.Solalp<- ecospat.testData[c("Soldanella_alpina","ddeg","mind","srad","slp","topo")]

# copy maxent.jar file in the right folder
path.from<-system.file("extdata", "maxent.txt", package="ecospat")
path.to <- paste0(system.file(package="dismo"), "/java/maxent.txt")
path.to.renamed <- paste0(system.file(package="dismo"), "/java/maxent.jar")
file.copy(path.from,path.to,overwrite = TRUE)
file.rename(path.to, path.to.renamed)

if (file.exists(path.to.renamed) & require(rJava)) {
  me <- maxent(data.Solalp[,-1],data.Solalp[,1])
  ecospat.maxentvarimport (model=me, dfvar=data.Solalp[,-1], nperm=5)
  }



cleanEx()
nameEx("ecospat.mdr")
### * ecospat.mdr

flush(stderr()); flush(stdout())

### Name: ecospat.mdr
### Title: Minimum Dispersal Routes)
### Aliases: ecospat.mdr

### ** Examples

library(maps)

data(ecospat.testMdr)
data<- ecospat.testMdr
intros<-order(data$date)[1:2] # rows corresponding to first introductions

# plot observed situation

plot(data[,2:1],pch=15,cex=0.5)
points(data[intros,2:1],pch=19,col="red")
text(data[,2]+0.5,data[,1]+0.5,data[,3],cex=0.5)
map(add=TRUE)

# calculate minimum cost arborescence (MCA) of dispersal routes

obs<-ecospat.mdr(data=data,xcol=2,ycol=1,datecol=3,mode="min",rep=100,
                  mean.date.error=1,fixed.sources.rows=intros)

# plot MCA results
# arrows' thickness indicate support for the routes

mca<-obs[[1]]
plot(mca[,3:4],type="n",xlab="longitude",ylab="latitude")
arrows(mca[,1],mca[,2],mca[,3],mca[,4],length = 0.05,lwd=mca$bootstrap.value*2)
map(add=TRUE)

# plot intros

points(data[intros,2:1],pch=19,col="red")
text(data[intros,2]+0.5,data[intros,1]+0.5,data[intros,3],cex=1,col="red")

# dispersal routes statistics

obs[[2]] # total routes length in DD
obs[[3]] # median dispersal rate in DD/yr
obs[[4]] # number of outcoming nodes



cleanEx()
nameEx("ecospat.mess")
### * ecospat.mess

flush(stderr()); flush(stdout())

### Name: ecospat.mess
### Title: MESS
### Aliases: ecospat.mess

### ** Examples

x <- ecospat.testData[c(2,3,4:8)]
proj <- x[1:90,] #A projection dataset.
cal <- x[91:300,] #A calibration dataset

#Create a MESS object 
mess.object <- ecospat.mess (proj, cal, w="default")

#Plot MESS 
ecospat.plot.mess (mess.object, cex=1, pch=15)



cleanEx()
nameEx("ecospat.meva.table")
### * ecospat.meva.table

flush(stderr()); flush(stdout())

### Name: ecospat.meva.table
### Title: Model Evaluation For A Given Threshold Value
### Aliases: ecospat.meva.table
### Keywords: file

### ** Examples


Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris

meva <- ecospat.meva.table (Pred, Sp.occ, 0.39)



cleanEx()
nameEx("ecospat.mpa")
### * ecospat.mpa

flush(stderr()); flush(stdout())

### Name: ecospat.mpa
### Title: Minimal Predicted Area
### Aliases: ecospat.mpa

### ** Examples

data(ecospat.testData)
obs <- (ecospat.testData$glm_Saxifraga_oppositifolia
[which(ecospat.testData$Saxifraga_oppositifolia==1)])

ecospat.mpa(obs)
ecospat.mpa(obs,perc=1) ## 100 percent of the presences encompassed



cleanEx()
nameEx("ecospat.niche.dynIndexProjGeo")
### * ecospat.niche.dynIndexProjGeo

flush(stderr()); flush(stdout())

### Name: ecospat.niche.dynIndexProjGeo
### Title: Projection of niche dynamic indices to the Geography
### Aliases: ecospat.niche.dynIndexProjGeo

### ** Examples



cleanEx()
nameEx("ecospat.niche.zProjGeo")
### * ecospat.niche.zProjGeo

flush(stderr()); flush(stdout())

### Name: ecospat.niche.zProjGeo
### Title: Projection of Occurrence Densities to the Geography
### Aliases: ecospat.niche.zProjGeo

### ** Examples



cleanEx()
nameEx("ecospat.npred")
### * ecospat.npred

flush(stderr()); flush(stdout())

### Name: ecospat.npred
### Title: Number Of Predictors
### Aliases: ecospat.npred

### ** Examples

colvar <- ecospat.testData[c(4:8)]
x <- cor(colvar, method="pearson")
ecospat.npred (x, th=0.75)



cleanEx()
nameEx("ecospat.occ.desaggregation")
### * ecospat.occ.desaggregation

flush(stderr()); flush(stdout())

### Name: ecospat.occ.desaggregation
### Title: Species Occurrences Desaggregation
### Aliases: ecospat.occ.desaggregation

### ** Examples





cleanEx()
nameEx("ecospat.occupied.patch")
### * ecospat.occupied.patch

flush(stderr()); flush(stdout())

### Name: ecospat.occupied.patch
### Title: Extract occupied patches of a species in geographic space.)
### Aliases: ecospat.occupied.patch
### Keywords: file

### ** Examples






cleanEx()
nameEx("ecospat.permut.glm")
### * ecospat.permut.glm

flush(stderr()); flush(stdout())

### Name: ecospat.permut.glm
### Title: GLM Permutation Function
### Aliases: ecospat.permut.glm

### ** Examples






cleanEx()
nameEx("ecospat.plot.kappa")
### * ecospat.plot.kappa

flush(stderr()); flush(stdout())

### Name: ecospat.plot.kappa
### Title: Plot Kappa
### Aliases: ecospat.plot.kappa
### Keywords: file

### ** Examples



Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris
ecospat.plot.kappa(Pred, Sp.occ)



cleanEx()
nameEx("ecospat.plot.mess")
### * ecospat.plot.mess

flush(stderr()); flush(stdout())

### Name: ecospat.plot.mess
### Title: Plot MESS
### Aliases: ecospat.plot.mess

### ** Examples




cleanEx()
nameEx("ecospat.plot.tss")
### * ecospat.plot.tss

flush(stderr()); flush(stdout())

### Name: ecospat.plot.tss
### Title: Plot True skill statistic (TSS)
### Aliases: ecospat.plot.tss
### Keywords: file

### ** Examples

Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris
ecospat.plot.tss(Pred, Sp.occ)



cleanEx()
nameEx("ecospat.rand.pseudoabsences")
### * ecospat.rand.pseudoabsences

flush(stderr()); flush(stdout())

### Name: ecospat.rand.pseudoabsences
### Title: Sample Pseudo-Absences
### Aliases: ecospat.rand.pseudoabsences

### ** Examples

glob <- ecospat.testData[2:8]
presence <- ecospat.testData[c(2:3,9)]
presence <- presence[presence[,3]==1,1:2]
ecospat.rand.pseudoabsences (nbabsences=10, glob=glob, colxyglob=1:2, colvar = "all", 
presence= presence, colxypresence=1:2, mindist=20)



cleanEx()
nameEx("ecospat.rangesize")
### * ecospat.rangesize

flush(stderr()); flush(stdout())

### Name: ecospat.rangesize
### Title: Quantification of the range size of a species using habitat
###   suitability maps and IUCN criteria)
### Aliases: ecospat.rangesize
### Keywords: file

### ** Examples




cleanEx()
nameEx("ecospat.rcls.grd")
### * ecospat.rcls.grd

flush(stderr()); flush(stdout())

### Name: ecospat.rcls.grd
### Title: Reclassifying grids function
### Aliases: ecospat.rcls.grd

### ** Examples


library(raster)
library(classInt)

bio3<- raster(system.file("external/bioclim/current/bio3.grd",package="biomod2"))
bio12<- raster(system.file("external/bioclim/current/bio12.grd",package="biomod2"))

B3.rcl<-ecospat.rcls.grd(bio3,9) 
B12.rcl<-ecospat.rcls.grd(bio12,9)
B3B12.comb <- B12.rcl+B3.rcl*10

# Plotting a histogram of the classes
hist(B3B12.comb,breaks=100,col=heat.colors(88)) 
# Plotting the new RasterLayer (9x9 classes)
plot(B3B12.comb,col=rev(rainbow(88)),main="Stratified map") 




cleanEx()
nameEx("ecospat.recstrat_prop")
### * ecospat.recstrat_prop

flush(stderr()); flush(stdout())

### Name: ecospat.recstrat_prop
### Title: Random Ecologically Stratified Sampling of propotional numbers
### Aliases: ecospat.recstrat_prop

### ** Examples

    library(raster)
    library(classInt)
    
    bio3<- raster(system.file("external/bioclim/current/bio3.grd",package="biomod2"))
    bio12<- raster(system.file("external/bioclim/current/bio12.grd",package="biomod2"))
    
    B3.rcl<-ecospat.rcls.grd(bio3,9) 
    B12.rcl<-ecospat.rcls.grd(bio12,9)
    B3B12.comb <- B12.rcl+B3.rcl*10
    
    B3B12.prop_samples <- ecospat.recstrat_prop(B3B12.comb,100)
    
    plot(B3B12.comb)
    points(B3B12.prop_samples$x,B3B12.prop_samples$y,pch=16,cex=0.6,col=B3B12.prop_samples$class)



cleanEx()
nameEx("ecospat.recstrat_regl")
### * ecospat.recstrat_regl

flush(stderr()); flush(stdout())

### Name: ecospat.recstrat_regl
### Title: Random Ecologically Stratified Sampling of equal numbers
### Aliases: ecospat.recstrat_regl

### ** Examples


  library(raster)
  library(classInt)

  bio3<- raster(system.file("external/bioclim/current/bio3.grd",package="biomod2"))
  bio12<- raster(system.file("external/bioclim/current/bio12.grd",package="biomod2"))
    
  B3.rcl<-ecospat.rcls.grd(bio3,9) 
  B12.rcl<-ecospat.rcls.grd(bio12,9)
  B3B12.comb <- B12.rcl+B3.rcl*10
    
  B3B12.regl_samples <- ecospat.recstrat_prop(B3B12.comb,100)
  
  plot(B3B12.comb)
  points(B3B12.regl_samples$x,B3B12.regl_samples$y,pch=16,cex=0.6,col=B3B12.regl_samples$class)



cleanEx()
nameEx("ecospat.sample.envar")
### * ecospat.sample.envar

flush(stderr()); flush(stdout())

### Name: ecospat.sample.envar
### Title: Sample Environmental Variables
### Aliases: ecospat.sample.envar

### ** Examples




cleanEx()
nameEx("ecospat.testData")
### * ecospat.testData

flush(stderr()); flush(stdout())

### Name: ecospat.testData
### Title: Test Data For The Ecospat package
### Aliases: ecospat.testData

### ** Examples

data(ecospat.testData)
str(ecospat.testData)
dim(ecospat.testData)
names(ecospat.testData)



cleanEx()
nameEx("ecospat.testEnvRaster")
### * ecospat.testEnvRaster

flush(stderr()); flush(stdout())

### Name: ecospat.testEnvRaster
### Title: Test Environmental Rasters for The Ecospat package
### Aliases: ecospat.testEnvRaster

### ** Examples

## Not run: 
##D fpath <- system.file("extdata", "ecospat.testEnvRaster.RData", package="ecospat")
##D load(fpath)
##D plot(env)
## End(Not run)




cleanEx()
nameEx("ecospat.testMdr")
### * ecospat.testMdr

flush(stderr()); flush(stdout())

### Name: ecospat.testMdr
### Title: Test Data For The ecospat.mdr function
### Aliases: ecospat.testMdr

### ** Examples

data(ecospat.testMdr)
str(ecospat.testMdr)
dim(ecospat.testMdr)



cleanEx()
nameEx("ecospat.testNiche")
### * ecospat.testNiche

flush(stderr()); flush(stdout())

### Name: ecospat.testNiche
### Title: Test Data For The Niche Overlap Analysis
### Aliases: ecospat.testNiche

### ** Examples

data(ecospat.testNiche)
dim(ecospat.testNiche)
names(ecospat.testNiche)



cleanEx()
nameEx("ecospat.testNiche.inv")
### * ecospat.testNiche.inv

flush(stderr()); flush(stdout())

### Name: ecospat.testNiche.inv
### Title: Test Data For The Niche Dynamics Analysis In The Invaded Range
###   Of A Hypothetical Species
### Aliases: ecospat.testNiche.inv

### ** Examples

data(ecospat.testNiche.inv)
str(ecospat.testNiche.inv)
dim(ecospat.testNiche.inv)
names(ecospat.testNiche.inv)



cleanEx()
nameEx("ecospat.testNiche.nat")
### * ecospat.testNiche.nat

flush(stderr()); flush(stdout())

### Name: ecospat.testNiche.nat
### Title: Test Data For The Niche Dynamics Analysis In The Native Range Of
###   A Hypothetical Species
### Aliases: ecospat.testNiche.nat

### ** Examples

data(ecospat.testNiche.nat)
str(ecospat.testNiche.nat)
dim(ecospat.testNiche.nat)
names(ecospat.testNiche.nat)



cleanEx()
nameEx("ecospat.testTree")
### * ecospat.testTree

flush(stderr()); flush(stdout())

### Name: ecospat.testTree
### Title: Test Tree For The Ecospat package
### Aliases: ecospat.testTree

### ** Examples

fpath <- system.file("extdata", "ecospat.testTree.tre", package="ecospat")
tree <- read.tree(fpath)
plot(tree)



cleanEx()
nameEx("ecospat.varpart")
### * ecospat.varpart

flush(stderr()); flush(stdout())

### Name: ecospat.varpart
### Title: Variation Partitioning For GLM Or GAM
### Aliases: ecospat.varpart

### ** Examples

library(rms)
data('ecospat.testData')

# data for Soldanella alpina and Achillea millefolium
data.Solalp<- ecospat.testData[c("Soldanella_alpina","ddeg","mind","srad","slp","topo")]

# glm models for Soldanella alpina

glm.Solalp1 <- glm("Soldanella_alpina ~ pol(ddeg,2) + pol(mind,2) + pol(srad,2)", 
                  data = data.Solalp, family = binomial)
glm.Solalp2 <- glm("Soldanella_alpina ~ pol(slp,2) + pol(topo,2)", 
                  data = data.Solalp, family = binomial)
                  
ecospat.varpart (model.1= glm.Solalp1, model.2= glm.Solalp2, model.12= glm.Solalp2)



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
