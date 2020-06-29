pkgname <- "statVisual"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('statVisual')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("BiAxisErrBar")
### * BiAxisErrBar

flush(stderr()); flush(stdout())

### Name: BiAxisErrBar
### Title: Compare Patterns of Two Outcomes in One Scatter Plot
### Aliases: BiAxisErrBar
### Keywords: method

### ** Examples


library(tidyverse)
library(ggplot2)

print(head(mtcars))

print(table(mtcars$gear, useNA="ifany"))

statVisual(type = "BiAxisErrBar",
  dat= mtcars,
  group = "gear",
  y.left = "mpg",
  y.right = "wt")



BiAxisErrBar(
  dat = mtcars,
  group = "gear",
  y.left = "mpg",
  y.right = "wt")




cleanEx()
nameEx("Box")
### * Box

flush(stderr()); flush(stdout())

### Name: Box
### Title: Compare Groups Based on Boxplots Across Time
### Aliases: Box
### Keywords: method

### ** Examples

library(dplyr)

data(longDat)

print(dim(longDat))
print(longDat[1:3,])

print(table(longDat$time, useNA = "ifany"))
print(table(longDat$grp, useNA = "ifany"))
print(table(longDat$sid, useNA = "ifany"))

print(table(longDat$time, longDat$grp))

statVisual(type = 'Box', 
           data = longDat, 
           x = 'time', 
           y = 'y', 
           group = 'grp',
	   title = "Boxplots across time") 

Box( 
    data = longDat, 
    x = 'time', 
    y = 'y', 
    group = 'grp',
    title = "Boxplots across time") 





cleanEx()
nameEx("BoxROC")
### * BoxROC

flush(stderr()); flush(stdout())

### Name: BoxROC
### Title: Compare Boxplots with ROC Curve
### Aliases: BoxROC
### Keywords: method

### ** Examples

library(dplyr)
library(gridExtra)

data(esSim)
print(esSim)

# expression data
dat = exprs(esSim)
print(dim(dat))
print(dat[1:2,])

# phenotype data
pDat = pData(esSim)
print(dim(pDat))
print(pDat[1:2,])

# feature data
fDat = fData(esSim)
print(dim(fDat))
print(fDat[1:2,])

# choose the first probe which is over-expressed in cases
pDat$probe1 = dat[1,]

# check histograms of probe 1 expression in cases and controls
print(table(pDat$grp, useNA = "ifany"))

statVisual(type = 'BoxROC', 
           data = pDat, 
           group = 'grp', 
           y = 'probe1', 
           point.size = 1)

BoxROC(
  data = pDat,
  group = 'grp', 
  y = 'probe1', 
  point.size = 1)




cleanEx()
nameEx("Den")
### * Den

flush(stderr()); flush(stdout())

### Name: Den
### Title: Compare Groups Based on Density Plots
### Aliases: Den
### Keywords: method

### ** Examples

data(esSim)
print(esSim)

# expression data
dat = exprs(esSim)
print(dim(dat))
print(dat[1:2,])

# phenotype data
pDat = pData(esSim)
print(dim(pDat))
print(pDat[1:2,])

# feature data
fDat = fData(esSim)
print(dim(fDat))
print(fDat[1:2,])

# choose the first probe which is over-expressed in cases
pDat$probe1 = dat[1,]

# check histograms of probe 1 expression in cases and controls
print(table(pDat$grp, useNA = "ifany"))

statVisual(type = 'Den', 
           data = pDat, 
           y = 'probe1', 
           group = 'grp') 

Den( 
    data = pDat, 
    y = 'probe1', 
    group = 'grp') 




cleanEx()
nameEx("Dendro")
### * Dendro

flush(stderr()); flush(stdout())

### Name: Dendro
### Title: Compare Groups Based on Dendrogram
### Aliases: Dendro
### Keywords: method

### ** Examples


data(esSim)
print(esSim)

# expression data
dat = exprs(esSim)
print(dim(dat))
print(dat[1:2,])

# phenotype data
pDat = pData(esSim)
print(dim(pDat))
print(pDat[1:2,])

# feature data
fDat = fData(esSim)
print(dim(fDat))
print(fDat[1:2,])

# choose the first 6 probes (3 OE probes, 2 UE probes, and 1 NE probe)
pDat$probe1 = dat[1,]
pDat$probe2 = dat[2,]
pDat$probe3 = dat[3,]
pDat$probe4 = dat[4,]
pDat$probe5 = dat[5,]
pDat$probe6 = dat[6,]

print(pDat[1:2, ])

# check histograms of probe 1 expression in cases and controls
print(table(pDat$grp, useNA = "ifany"))

pDat$grp = factor(pDat$grp)

statVisual(type = 'Dendro', 
           x = pDat[, c(3:8)], 
           group = pDat$grp)

Dendro(
       x = pDat[, c(3:8)], 
       group = pDat$grp)






cleanEx()
nameEx("ErrBar")
### * ErrBar

flush(stderr()); flush(stdout())

### Name: ErrBar
### Title: Compare Groups Based on dotplots Across Time
### Aliases: ErrBar
### Keywords: method

### ** Examples



data(longDat)

print(dim(longDat))
print(longDat[1:3,])

print(table(longDat$time, useNA = "ifany"))
print(table(longDat$grp, useNA = "ifany"))
print(table(longDat$sid, useNA = "ifany"))

print(table(longDat$time, longDat$grp))

statVisual(type = 'ErrBar', 
  data = longDat, 
  x = 'time', 
  y = 'y', 
  group = 'grp',
  title = "Dot plots across time") 


ErrBar(
  data = longDat, 
  x = 'time', 
  y = 'y', 
  group = 'grp',
  title = "Dot plots across time") 





cleanEx()
nameEx("Heat")
### * Heat

flush(stderr()); flush(stdout())

### Name: Heat
### Title: Heatmap with Row Names Colored by Group
### Aliases: Heat
### Keywords: method

### ** Examples


data(esSim)
print(esSim)

# expression data
dat = exprs(esSim)
print(dim(dat))
print(dat[1:2,])

# phenotype data
pDat = pData(esSim)
print(dim(pDat))
print(pDat[1:2,])

# feature data
fDat = fData(esSim)
print(dim(fDat))
print(fDat[1:2,])

# choose the first 6 probes (3 OE probes, 2 UE probes, and 1 NE probe)
pDat$probe1 = dat[1,]
pDat$probe2 = dat[2,]
pDat$probe3 = dat[3,]
pDat$probe4 = dat[4,]
pDat$probe5 = dat[5,]
pDat$probe6 = dat[6,]

print(pDat[1:2, ])

# check histograms of probe 1 expression in cases and controls
print(table(pDat$grp, useNA = "ifany"))

statVisual(type = 'Heat', 
           data = pDat[, c(2:8)], 
           group = 'grp')

Heat(
     data = pDat[, c(2:8)], 
     group = 'grp')




cleanEx()
nameEx("Hist")
### * Hist

flush(stderr()); flush(stdout())

### Name: Hist
### Title: Compare Groups Based on Histograms
### Aliases: Hist
### Keywords: method

### ** Examples

data(esSim)
print(esSim)

# expression data
dat = exprs(esSim)
print(dim(dat))
print(dat[1:2,])

# phenotype data
pDat = pData(esSim)
print(dim(pDat))
print(pDat[1:2,])

# feature data
fDat = fData(esSim)
print(dim(fDat))
print(fDat[1:2,])

# choose the first probe which is over-expressed in cases
pDat$probe1 = dat[1,]

# check histograms of probe 1 expression in cases and controls
print(table(pDat$grp, useNA = "ifany"))

statVisual(type = 'Hist', 
       data = pDat, 
       y = 'probe1', 
       group = 'grp') 

Hist(
     data = pDat, 
     y = 'probe1', 
     group = 'grp') 





cleanEx()
nameEx("ImpPlot")
### * ImpPlot

flush(stderr()); flush(stdout())

### Name: ImpPlot
### Title: Plot of Variable Importance
### Aliases: ImpPlot
### Keywords: method

### ** Examples


library(dplyr)
library(randomForest)
library(tibble)


data(esSim)
print(esSim)

# expression data
dat = exprs(esSim)
print(dim(dat))
print(dat[1:2,])

# phenotype data
pDat = pData(esSim)
print(dim(pDat))
print(pDat[1:2,])

# feature data
fDat = fData(esSim)
print(dim(fDat))
print(fDat[1:2,])

# choose the first 6 probes (3 OE probes, 2 UE probes, and 1 NE probe)
pDat$probe1 = dat[1,]
pDat$probe2 = dat[2,]
pDat$probe3 = dat[3,]
pDat$probe4 = dat[4,]
pDat$probe5 = dat[5,]
pDat$probe6 = dat[6,]

print(pDat[1:2, ])

# check histograms of probe 1 expression in cases and controls
print(table(pDat$grp, useNA = "ifany"))

pDat$grp = factor(pDat$grp)


rf_m = randomForest(
  x = pDat[, c(3:8)], 
  y = pDat$grp, 
  importance = TRUE, proximity = TRUE
)


statVisual(type = 'ImpPlot', model = rf_m)

ImpPlot(model = rf_m)




cleanEx()
nameEx("LinePlot")
### * LinePlot

flush(stderr()); flush(stdout())

### Name: LinePlot
### Title: Compare Groups Based on Trajectory Plots
### Aliases: LinePlot
### Keywords: method

### ** Examples

data(longDat)

print(dim(longDat))
print(longDat[1:3,])

print(table(longDat$time, useNA = "ifany"))
print(table(longDat$grp, useNA = "ifany"))
print(table(longDat$sid, useNA = "ifany"))

print(table(longDat$time, longDat$grp))

statVisual(type = "LinePlot",
  data = longDat,
  x = 'time',
  y = 'y',
  sid = 'sid',
  group = 'grp')

LinePlot(
  data = longDat,
  x = 'time',
  y = 'y',
  sid = 'sid',
  group = 'grp')




cleanEx()
nameEx("PCA_score")
### * PCA_score

flush(stderr()); flush(stdout())

### Name: PCA_score
### Title: Scatter Plot of 2 Specified Principal Components
### Aliases: PCA_score
### Keywords: method

### ** Examples

library(factoextra)

data(esSim)
print(esSim)

# expression data
dat = exprs(esSim)
print(dim(dat))
print(dat[1:2,])

# phenotype data
pDat = pData(esSim)
print(dim(pDat))
print(pDat[1:2,])

# feature data
fDat = fData(esSim)
print(dim(fDat))
print(fDat[1:2,])

# choose the first 6 probes (3 OE probes, 2 UE probes, and 1 NE probe)
pDat$probe1 = dat[1,]
pDat$probe2 = dat[2,]
pDat$probe3 = dat[3,]
pDat$probe4 = dat[4,]
pDat$probe5 = dat[5,]
pDat$probe6 = dat[6,]

print(pDat[1:2, ])

# check histograms of probe 1 expression in cases and controls
print(table(pDat$grp, useNA = "ifany"))

pDat$grp = factor(pDat$grp)

###

pca.obj = iprcomp(pDat[, c(3:8)], scale. = TRUE)

# scree plot
factoextra::fviz_eig(pca.obj, addlabels = TRUE)

# scatter plot of PC1 vs PC2
statVisual(type = 'PCA_score',
           prcomp_obj = pca.obj, 
           dims = c(1, 2),
           data = pDat, 
           color = 'grp',
           loadings = FALSE)

PCA_score(prcomp_obj = pca.obj, 
          dims = c(1, 3),
          data = pDat, 
          color = 'grp',
          loadings = FALSE)




cleanEx()
nameEx("PVCA")
### * PVCA

flush(stderr()); flush(stdout())

### Name: PVCA
### Title: Principal Variance Component Analysis (PVCA)
### Aliases: PVCA
### Keywords: method

### ** Examples

library(pvca)


data(esSim)
print(esSim)

# expression data
dat = exprs(esSim)
print(dim(dat))
print(dat[1:2,])

# create a fake Batch variable
esSim$Batch=c(rep("A", 4), rep("B", 6), rep("C", 10))
# phenotype data
pDat = pData(esSim)
print(dim(pDat))
print(pDat[1:2,])


# feature data
fDat = fData(esSim)
print(dim(fDat))
print(fDat[1:2,])


statVisual(type = 'PVCA',
           clin_data = pData(esSim), 
           clin_subjid = "sid", 
           gene_data = exprs(esSim), 
           batch.factors = c("grp", "Batch"))

PVCA( 
     clin_data = pData(esSim), 
     clin_subjid = "sid", 
     gene_data = exprs(esSim), 
     batch.factors = c("grp", "Batch"))




cleanEx()
nameEx("Volcano")
### * Volcano

flush(stderr()); flush(stdout())

### Name: Volcano
### Title: Volcano Plot
### Aliases: Volcano
### Keywords: method

### ** Examples

library(ggrepel)
library(limma)

# load the simulated dataset
data(esSim)
print(esSim)

# expression levels
y = exprs(esSim)
print(dim(y))
print(y[1:2,])

# phenotype data
pDat = pData(esSim)
print(dim(pDat))
print(pDat)

# design matrix
design = model.matrix(~grp, data = pDat)
print(design)

options(digits = 3)

# Ordinary fit
fit <- lmFit(y, design)
fit2 <- eBayes(fit)

# get result data frame
resFrame = topTable(fit2,coef = 2, number = nrow(esSim))
print(dim(resFrame))
print(resFrame[1:2,])
resFrame$sigFlag  =  resFrame$adj.P.Val < 0.05

resFrame$probe  =  rownames(resFrame)
# make sure set NA to genes non-differentially expressed
resFrame$probe[which(resFrame$sigFlag == FALSE)] = NA

print(resFrame[1:2,])
print(table(resFrame$sigFlag, useNA = "ifany"))

statVisual(type = 'Volcano',
           resFrame = resFrame, 
           stats = 'logFC', 
           p.value = 'P.Value', 
           group = 'sigFlag', 
           rowname.var = 'probe', 
           point.size = 1)

Volcano(
  resFrame = resFrame, 
  stats = 'logFC', 
  p.value = 'P.Value', 
  group = 'sigFlag', 
  rowname.var = 'probe', 
  point.size = 1)





cleanEx()
nameEx("XYscatter")
### * XYscatter

flush(stderr()); flush(stdout())

### Name: XYscatter
### Title: Compare Groups Based on Scatter Plots
### Aliases: XYscatter
### Keywords: method

### ** Examples

data(diffCorDat)

print(dim(diffCorDat))
print(diffCorDat[1:2,])

statVisual(type = 'XYscatter',
  data = diffCorDat, 
  x = 'probe1', 
  y = 'probe2', 
  group = 'grp', 
  title = 'Scatter Plot: probe1 vs probe2')

XYscatter( 
  data = diffCorDat, 
  x = 'probe1', 
  y = 'probe2', 
  group = 'grp', 
  title = 'Scatter Plot: probe1 vs probe2')




cleanEx()
nameEx("barPlot")
### * barPlot

flush(stderr()); flush(stdout())

### Name: barPlot
### Title: Compare Groups Based on Barplots Across Time
### Aliases: barPlot
### Keywords: method

### ** Examples



data(longDat)

print(dim(longDat))
print(longDat[1:3,])

print(table(longDat$time, useNA = "ifany"))
print(table(longDat$grp, useNA = "ifany"))
print(table(longDat$sid, useNA = "ifany"))

print(table(longDat$time, longDat$grp))

statVisual(type = 'barPlot', 
  data = longDat, 
  x = 'time', 
  y = 'y', 
  group = 'grp',
  title = "Bar plots across time") 


barPlot(
  data = longDat, 
  x = 'time', 
  y = 'y', 
  group = 'grp',
  title = "Bar plots across time") 





cleanEx()
nameEx("cv_glmnet_plot")
### * cv_glmnet_plot

flush(stderr()); flush(stdout())

### Name: cv_glmnet_plot
### Title: Plot the Cross-Validation Curve Produced by cv.glmnet
### Aliases: cv_glmnet_plot
### Keywords: method

### ** Examples

library(dplyr)
library(tibble)
library(glmnet)

data(esSim)
print(esSim)

# expression data
dat = exprs(esSim)
print(dim(dat))
print(dat[1:2,])

# phenotype data
pDat = pData(esSim)
print(dim(pDat))
print(pDat[1:2,])

# feature data
fDat = fData(esSim)
print(dim(fDat))
print(fDat[1:2,])

# choose the first 6 probes (3 OE probes, 2 UE probes, and 1 NE probe)
pDat$probe1 = dat[1,]
pDat$probe2 = dat[2,]
pDat$probe3 = dat[3,]
pDat$probe4 = dat[4,]
pDat$probe5 = dat[5,]
pDat$probe6 = dat[6,]

print(pDat[1:2, ])

# check histograms of probe 1 expression in cases and controls
print(table(pDat$grp, useNA = "ifany"))


statVisual(type = "cv_glmnet_plot",
           x = as.matrix(pDat[, c(3:8)]), 
           y = pDat$grp, 
           family = "binomial")

cv_glmnet_plot(x = as.matrix(pDat[, c(3:8)]), 
               y = pDat$grp, 
               family = "binomial")




cleanEx()
nameEx("diffCorDat")
### * diffCorDat

flush(stderr()); flush(stdout())

### Name: diffCorDat
### Title: A Dataset for Differential Correlation Analysis
### Aliases: diffCorDat
### Keywords: datasets

### ** Examples

data(diffCorDat)

print(dim(diffCorDat))
print(diffCorDat[1:2,])



cleanEx()
nameEx("esSim")
### * esSim

flush(stderr()); flush(stdout())

### Name: esSim
### Title: A Simulated Gene Expression Dataset
### Aliases: esSim
### Keywords: datasets

### ** Examples

data(esSim)

print(esSim)

###
dat=exprs(esSim)
print(dim(dat))
print(dat[1:2,])

###
pDat=pData(esSim)
print(dim(pDat))
print(pDat)

# subject group status
print(table(esSim$grp))

###
fDat = fData(esSim)
print(dim(fDat))
print(fDat[1:2, ])

# probe's status of differential expression
print(table(fDat$memProbes))




cleanEx()
nameEx("genoSim")
### * genoSim

flush(stderr()); flush(stdout())

### Name: genoSim
### Title: An ExpressionSet Object Storing Simulated Genotype Data
### Aliases: genoSim
### Keywords: datasets

### ** Examples

data(genoSim)

print(genoSim)



cleanEx()
nameEx("iprcomp")
### * iprcomp

flush(stderr()); flush(stdout())

### Name: iprcomp
### Title: Improved Function for Obtaining Principal Components
### Aliases: iprcomp
### Keywords: method

### ** Examples

# generate simulated data
set.seed(1234567)
dat.x = matrix(rnorm(500), nrow = 100, ncol = 5)
dat.y = matrix(rnorm(500, mean = 2), nrow = 100, ncol = 5)
dat = rbind(dat.x, dat.y)
grp = c(rep(0, 100), rep(1, 100))
print(dim(dat))

res = iprcomp(dat, center = TRUE, scale.  =  FALSE)

# for each row, set one artificial missing value
dat.na=dat
nr=nrow(dat.na)
nc=ncol(dat.na)
for(i in 1:nr)
{
  posi=sample(x=1:nc, size=1)
  dat.na[i,posi]=NA
}

res.na = iprcomp(dat.na, center = TRUE, scale.  =  FALSE)

##
# pca plot
##
par(mfrow = c(3,1))
# original data without missing values
plot(x = res$x[,1], y = res$x[,2], xlab = "PC1", ylab  =  "PC2")
# perturbed data with one NA per probe 
# the pattern of original data is captured
plot(x = res.na$x[,1], y = res.na$x[,2], xlab = "PC1", ylab  =  "PC2", main = "with missing values")
par(mfrow = c(1,1))




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("longDat")
### * longDat

flush(stderr()); flush(stdout())

### Name: longDat
### Title: A Simulated Dataset for Longitudinal Data Analysis
### Aliases: longDat
### Keywords: datasets

### ** Examples

data(longDat)

print(dim(longDat))
print(longDat[1:3,])

print(table(longDat$time, useNA = "ifany"))
print(table(longDat$grp, useNA = "ifany"))
print(table(longDat$sid, useNA = "ifany"))

print(table(longDat$time, longDat$grp))



cleanEx()
nameEx("stackedBarPlot")
### * stackedBarPlot

flush(stderr()); flush(stdout())

### Name: stackedBarPlot
### Title: Draw Stacked Bar Plots
### Aliases: stackedBarPlot
### Keywords: method

### ** Examples

data(genoSim)

pDat = pData(genoSim)
geno = exprs(genoSim)

pDat$snp1 = geno[1,]

print(table(pDat$snp1, pDat$grp, useNA="ifany"))

stackedBarPlot(dat = pDat, 
	       catVar = "snp1", 
	       group = "grp", 
               xlab = "snp1", 
	       ylab = "Count", 
	       group.lab = "grp",
               title = "Stacked barplots of counts",
               catVarLevel = NULL)




cleanEx()
nameEx("statVisual")
### * statVisual

flush(stderr()); flush(stdout())

### Name: statVisual
### Title: The Wrapper Function Incorporating All Wrapper Functions in
###   statVisual
### Aliases: statVisual
### Keywords: method

### ** Examples

data(esSim)
print(esSim)

# expression data
dat = exprs(esSim)
print(dim(dat))
print(dat[1:2,])

# phenotype data
pDat = pData(esSim)
print(dim(pDat))
print(pDat[1:2,])

# feature data
fDat = fData(esSim)
print(dim(fDat))
print(fDat[1:2,])

# choose the first probe which is over-expressed in cases
pDat$probe1 = dat[1,]

# check histograms of probe 1 expression in cases and controls
print(table(pDat$grp, useNA = "ifany"))

statVisual(type = 'Hist', 
       data = pDat, 
       y = 'probe1', 
       group = 'grp') 





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
