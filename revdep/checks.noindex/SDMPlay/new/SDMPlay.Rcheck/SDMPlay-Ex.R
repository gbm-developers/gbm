pkgname <- "SDMPlay"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SDMPlay')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("SDMdata.quality")
### * SDMdata.quality

flush(stderr()); flush(stdout())

### Name: SDMdata.quality
### Title: Evaluate dataset quality
### Aliases: SDMdata.quality

### ** Examples

#Open SDMtab object example
x <- system.file ("extdata","SDMdata1500.csv", package="SDMPlay")
SDMdata <- read.table(x,header=TRUE, sep=";")

# Evaluate the dataset
SDMPlay:::SDMdata.quality(data=SDMdata)



cleanEx()
nameEx("SDMeval")
### * SDMeval

flush(stderr()); flush(stdout())

### Name: SDMeval
### Title: Evaluate species distribution models
### Aliases: SDMeval

### ** Examples

# Model example
load(system.file('extdata','model.RData',package='SDMPlay'))
modelBRT <- model

# Evaluate modelling performance
#SDMPlay:::SDMeval(modelBRT)



cleanEx()
nameEx("SDMtab")
### * SDMtab

flush(stderr()); flush(stdout())

### Name: SDMtab
### Title: Compile species distribution dataset for modelling
### Aliases: SDMtab

### ** Examples

#Open occurrence data
data('ctenocidaris.nutrix')
occ <- ctenocidaris.nutrix

#Open environmental descriptors RasterStack
r <- raster:: stack(system.file('extdata', 'pred.grd',package='SDMPlay'))

#create the dataframe for modelling
z <- SDMPlay:::SDMtab(xydata=occ[,c('decimal.Longitude','decimal.Latitude')],predictors=r)
head(z)




cleanEx()
nameEx("brisaster.antarcticus")
### * brisaster.antarcticus

flush(stderr()); flush(stdout())

### Name: brisaster.antarcticus
### Title: Records of _Brisaster antarcticus_ echinoid presences on the
###   Kerguelen Plateau
### Aliases: brisaster.antarcticus
### Keywords: datasets

### ** Examples

data(brisaster.antarcticus)
x <- brisaster.antarcticus

# plot of the occurrences:
# selecting the species according to the campaigns
brisaster7475 <- subset(x,x$year==1974 | x$year==1975)
brisaster20102015 <- subset(x,x$campaign=='POKER II'| x$campaign=='PROTEKER')

# drawing the background (depth)
library(grDevices)
blue.palette <- colorRampPalette(c('blue','deepskyblue','azure'))(100)
data('predictors1965_1974')
depth <- raster :: subset(predictors1965_1974, 1)

raster::plot(depth, col=blue.palette,main= "Brisaster antarcticus occurrences")

# adding the occurrence data to the background
points(brisaster7475[,c('decimal.Longitude','decimal.Latitude')],
      col='orange',pch=16)
points(brisaster20102015[,c('decimal.Longitude','decimal.Latitude')],
      col='darkgreen',pch=16)
legend('bottomleft',
       legend=c('Brisaster antarcticus 1974-1975','Brisaster antarcticus 2010-2015'),
       col= c('orange','darkgreen'), pch= c(15, 15),cex=0.9)





cleanEx()
nameEx("compute.brt")
### * compute.brt

flush(stderr()); flush(stdout())

### Name: compute.brt
### Title: Compute BRT (Boosted Regression Trees) model
### Aliases: compute.brt

### ** Examples

## Not run: 
##D #Download the presence data
##D data('ctenocidaris.nutrix')
##D occ <- ctenocidaris.nutrix
##D # select longitude and latitude coordinates among all the information
##D occ <- ctenocidaris.nutrix[,c('decimal.Longitude','decimal.Latitude')]
##D 
##D #Download an example of environmental predictors
##D #restricted on geographical extent and depth (-1500m)
##D envi <- raster::stack(system.file('extdata', 'pred.grd',package='SDMPlay'))
##D envi
##D 
##D #Open SDMtab matrix
##D x <- system.file(file='extdata/SDMdata1500.csv',package='SDMPlay')
##D SDMdata <- read.table(x,header=TRUE, sep=';')
##D 
##D #Run the model
##D model <- SDMPlay:::compute.brt (x=SDMdata, proj.predictors=envi,lr=0.0005)
##D 
##D #Plot the partial dependance plots
##D dismo::gbm.plot(model$response)
##D 
##D #Get the contribution of each variable for the model
##D model$response$contributions
##D 
##D #Get the interaction between variables
##D dismo::gbm.interactions(model$response)
##D #Plot the interactions
##D int <- dismo::gbm.interactions(model$response)
##D # choose the interaction to plot
##D dismo::gbm.perspec(model$response,int$rank.list[1,1],int$rank.list[1,3])
##D 
##D #Plot the map prediction
##D library(grDevices) # add nice colors
##D palet.col <- colorRampPalette(c('deepskyblue','green','yellow', 'red'))( 80 )
##D raster::plot(model$raster.prediction, col=palet.col)
##D #add data
##D points(occ, col='black',pch=16)
## End(Not run)



cleanEx()
nameEx("compute.maxent")
### * compute.maxent

flush(stderr()); flush(stdout())

### Name: compute.maxent
### Title: Compute MaxEnt model
### Aliases: compute.maxent

### ** Examples

#Download the presence data
data('ctenocidaris.nutrix')
occ <- ctenocidaris.nutrix
# select longitude and latitude coordinates among all the information
occ <- ctenocidaris.nutrix[,c('decimal.Longitude','decimal.Latitude')]

#Download an example of environmental predictors
#restricted on geographical extent and depth (-1500m)
envi <- raster::stack(system.file('extdata', 'pred.grd',package='SDMPlay'))
envi

#Open SDMtab matrix
x <- system.file(file='extdata/SDMdata1500.csv',package='SDMPlay')
SDMdata <- read.table(x,header=TRUE, sep=';')

#only run if the maxent.jar file is available, in the right folder
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
# Check first if maxent can be run (normally not part of your script)
# (file.exists(jar) & require(rJava)) == TRUE ??
# rJava may pose a problem to load automatically within this package
# please load it manually using eventually the archives available from CRAN

# Run the model
#model <- SDMPlay:::compute.maxent(x=SDMdata , proj.predictors=envi)

# Plot the map prediction
library(grDevices) # add nice colors
#palet.col <- colorRampPalette(c('deepskyblue','green','yellow','red'))(80)
#'raster::plot(model$raster.prediction, col=palet.col)
# add data
#points(occ, col='black',pch=16)

# Get the partial dependance curves
#dismo::response(model$response)

# Get the percentage of contribution of each variable to the model
#plot(model$response)

# Get all the information provided by the model on a html document
#model$response




cleanEx()
nameEx("ctenocidaris.nutrix")
### * ctenocidaris.nutrix

flush(stderr()); flush(stdout())

### Name: ctenocidaris.nutrix
### Title: Records of _Ctenocidaris nutrix_ echinoid presences on the
###   Kerguelen Plateau
### Aliases: ctenocidaris.nutrix
### Keywords: datasets

### ** Examples

data(ctenocidaris.nutrix)
x <- ctenocidaris.nutrix
# plot of the occurrences:
# selecting the species according to the campaigns
ctenocidaris7475 <- base::subset(x,x$year==1974 | x$year==1975)
ctenocidaris20102015 <- base::subset(x,x$campaign=='POKER II' | x$campaign=='PROTEKER')

# drawing the background (depth)
library(grDevices)
blue.palette <- colorRampPalette(c('blue','deepskyblue','azure'))(100)
data('predictors1965_1974')
depth <- raster :: subset(predictors1965_1974, 1)

raster::plot(depth, col=blue.palette,main= "Ctenocidaris nutrix occurrences")

# adding the occurrences data to the background
points(ctenocidaris7475[,c('decimal.Longitude','decimal.Latitude')],
      col='orange',pch=16)
points(ctenocidaris20102015[,c('decimal.Longitude','decimal.Latitude')],
      col='darkgreen',pch=16)
legend('bottomleft',
       legend=c('Ctenocidaris nutrix 1974-1975','Ctenocidaris nutrix 2010-2015'),
       col= c('orange','darkgreen'), pch= c(15, 15),cex=0.9)





cleanEx()
nameEx("delim.area")
### * delim.area

flush(stderr()); flush(stdout())

### Name: delim.area
### Title: RasterStack preparation for modelling
### Aliases: delim.area

### ** Examples

data('predictors2005_2012')
envi <- predictors2005_2012

r <- SDMPlay:::delim.area(predictors = envi,
longmin = 70,longmax = 75, latmin = -50,latmax = -40,interval = c(0,-1000))
r

library(grDevices) # plot the result with nice colors
palet.col <- colorRampPalette(c('deepskyblue','green','yellow', 'red'))(80)
raster::plot(r, col=palet.col)



cleanEx()
nameEx("null.model")
### * null.model

flush(stderr()); flush(stdout())

### Name: null.model
### Title: Compute null model
### Aliases: null.model

### ** Examples

## Not run: 
##D library(dismo)
##D #Download the environmental predictors restricted on geographical extent and depth (-1500m)
##D envi <-raster::stack(system.file('extdata', 'pred.grd',package='SDMPlay'))
##D 
##D # Realize a null model type #2 with BRT
##D #--------------------------------------
##D # NB: the following arguments chosen for the example are not relevant,
##D # in the scope to minimize running time
##D modelN2 <- SDMPlay:::null.model(xy=NULL,predictors=envi,type=2,algorithm='brt',
##D                      nb=300,unique.data=TRUE, same=TRUE, nb.rep=2,lr=0.005)
##D 
##D # Look at the inputs used to implement the model
##D modelN2$input
##D 
##D # Get the evaluation of the models produced
##D modelN2$eval
##D 
##D # Get the evaluation of the mean of all these produced models (i.e. evaluation
##D # of the null model )
##D modelN2$eval.null
##D 
##D # Get the values of Spearman correlations between the all the prediction maps produced
##D modelN2$correlation
##D 
##D # Plot the mean null model map with nice colors
##D library(grDevices)
##D palet.col <- colorRampPalette(c('deepskyblue','green','yellow', 'red'))(80)
##D raster::plot(modelN2$pred.mean, col=palet.col)
## End(Not run)



cleanEx()
nameEx("predictors1965_1974")
### * predictors1965_1974

flush(stderr()); flush(stdout())

### Name: predictors1965_1974
### Title: Environmental descriptors for 1965-1974 on the Kerguelen Plateau
### Aliases: predictors1965_1974
### Keywords: datasets

### ** Examples

data(predictors1965_1974)
raster::plot(predictors1965_1974)




cleanEx()
nameEx("predictors2005_2012")
### * predictors2005_2012

flush(stderr()); flush(stdout())

### Name: predictors2005_2012
### Title: Environmental descriptors for 2005-2012 on the Kerguelen Plateau
### Aliases: predictors2005_2012
### Keywords: datasets

### ** Examples

data(predictors2005_2012)
raster::plot(predictors2005_2012)




cleanEx()
nameEx("predictors2200AIB")
### * predictors2200AIB

flush(stderr()); flush(stdout())

### Name: predictors2200AIB
### Title: IPCC environmental descriptors predicted for 2200 (AIB scenario)
###   on the Kerguelen Plateau
### Aliases: predictors2200AIB
### Keywords: datasets

### ** Examples

data(predictors2200AIB)
raster :: plot(predictors2200AIB)




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
