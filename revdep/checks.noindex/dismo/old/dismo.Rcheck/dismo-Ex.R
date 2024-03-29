pkgname <- "dismo"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('dismo')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("bioclim")
### * bioclim

flush(stderr()); flush(stdout())

### Name: bioclim
### Title: Bioclim
### Aliases: bioclim bioclim,Raster,matrix-method
###   bioclim,Raster,SpatialPoints-method bioclim,Raster,data.frame-method
###   bioclim,SpatialGridDataFrame,matrix-method
###   bioclim,SpatialGridDataFrame,SpatialPoints-method
###   bioclim,matrix,missing-method bioclim,data.frame,missing-method
###   Bioclim-class
### Keywords: spatial

### ** Examples

logo <- stack(system.file("external/rlogo.grd", package="raster"))
#presence data
pts <- matrix(c(48.243420, 48.243420, 47.985820, 52.880230, 49.531423, 46.182616, 54.168232, 
  69.624263, 83.792291, 85.337894, 74.261072, 83.792291, 95.126713, 84.565092, 66.275456, 41.803408,
  25.832176, 3.936132, 18.876962, 17.331359,7.048974, 13.648543, 26.093446, 28.544714, 39.104026, 
  44.572240, 51.171810, 56.262906, 46.269272, 38.161230, 30.618865, 21.945145, 34.390047, 59.656971,
  69.839163, 73.233228, 63.239594, 45.892154, 43.252326, 28.356155) , ncol=2)
bc <- bioclim(logo, pts)

#or
v <- extract(logo, pts)
bc <- bioclim(v)
p1 <- predict(logo, bc)
p2 <- predict(logo, bc, tails=c('both', 'low', 'high'))

#or
#sp <- SpatialPoints(pts)
#bc <- bioclim(logo, pts)



cleanEx()
nameEx("biovars")
### * biovars

flush(stderr()); flush(stdout())

### Name: biovars
### Title: bioclimatic variables
### Aliases: biovars biovars,matrix,matrix,matrix-method
###   biovars,Raster,Raster,Raster-method
###   biovars,vector,vector,vector-method
### Keywords: spatial

### ** Examples

tmin <- c(10,12,14,16,18,20,22,21,19,17,15,12)
tmax <- tmin + 5
prec <- c(0,2,10,30,80,160,80,20,40,60,20,0)
biovars(prec, tmin, tmax)

tmn = tmx = prc = brick(nrow=1, ncol=1)
tmn <- setValues(tmn, t(matrix(c(10,12,14,16,18,20,22,21,19,17,15,12))))
tmx <- tmn + 5
prc <- setValues(prc, t(matrix(c(0,2,10,30,80,160,80,20,40,60,20,0))))
b <- biovars(prc, tmn, tmx)
as.matrix(b)



cleanEx()
nameEx("circleHull")
### * circleHull

flush(stderr()); flush(stdout())

### Name: circleHull
### Title: Circle hull model
### Aliases: circleHull circleHull,SpatialPoints-method
###   circleHull,matrix-method circleHull,data.frame-method
###   CircleHull-class
### Keywords: spatial

### ** Examples

r <- raster(system.file("external/rlogo.grd", package="raster"))
#presence data
pts <- matrix(c(17, 42, 85, 70, 19, 53, 26, 84, 84, 46, 48, 85, 4, 95, 48, 54, 66, 
 74, 50, 48, 28, 73, 38, 56, 43, 29, 63, 22, 46, 45, 7, 60, 46, 34, 14, 51, 70, 31, 39, 26), ncol=2)
train <- pts[1:12, ]
test <- pts[13:20, ]
				 
ch <- circleHull(train)
predict(ch, test)

plot(r)
plot(ch, border='red', lwd=2, add=TRUE)
points(train, col='red', pch=20, cex=2)
points(test, col='black', pch=20, cex=2)

pr <- predict(ch, r, progress='')
plot(pr)
points(test, col='black', pch=20, cex=2)
points(train, col='red', pch=20, cex=2)

# to get the polygons:
p <- polygons(ch)
p



cleanEx()
nameEx("circles")
### * circles

flush(stderr()); flush(stdout())

### Name: CirclesRange
### Title: Circles range
### Aliases: circles circles,SpatialPoints-method circles,matrix-method
###   circles,data.frame-method CirclesRange-class
### Keywords: spatial

### ** Examples

r <- raster(system.file("external/rlogo.grd", package="raster"))
#presence data
pts <- matrix(c(17, 42, 85, 70, 19, 53, 26, 84, 84, 46, 48, 85, 4, 95, 48, 54, 66,
 74, 50, 48, 28, 73, 38, 56, 43, 29, 63, 22, 46, 45, 7, 60, 46, 34, 14, 51, 70, 31, 39, 26), ncol=2)
train <- pts[1:12, ]
test <- pts[13:20, ]
				 
cc <- circles(train, lonlat=FALSE)
predict(cc, test)

plot(r)
plot(cc, border='red', lwd=2, add=TRUE)
points(train, col='red', pch=20, cex=2)
points(test, col='black', pch=20, cex=2)

pr <- predict(cc, r, progress='')
plot(r, legend=FALSE)
plot(pr, add=TRUE, col='blue')
points(test, col='black', pch=20, cex=2)
points(train, col='red', pch=20, cex=2)


# to get the polygons:
p <- polygons(cc)
p

# to compute the Circular Area Range of a species (see Hijmans and Spooner, 2001)
pts <- train*10
ca100 <- polygons(circles(pts, d=100, lonlat=FALSE))
ca150 <- polygons(circles(pts, d=150, lonlat=FALSE))
sum(area(ca150)) / (pi*150^2)
sum(area(ca100)) / (pi*100^2)
par(mfrow=c(1,2))
plot(ca100); points(pts)
plot(ca150); points(pts)






graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("convexHull")
### * convexHull

flush(stderr()); flush(stdout())

### Name: Convex Hull
### Title: Convex hull model
### Aliases: convHull convHull,SpatialPoints-method convHull,matrix-method
###   convHull,data.frame-method ConvexHull-class
### Keywords: spatial

### ** Examples

r <- raster(system.file("external/rlogo.grd", package="raster"))
#presence data
pts <- matrix(c(17, 42, 85, 70, 19, 53, 26, 84, 84, 46, 48, 85, 4, 95, 48, 54, 66, 
 74, 50, 48, 28, 73, 38, 56, 43, 29, 63, 22, 46, 45, 7, 60, 46, 34, 14, 51, 70, 31, 39, 26), ncol=2)
train <- pts[1:12, ]
test <- pts[13:20, ]
				 
ch <- convHull(train)
predict(ch, test)

plot(r)
plot(ch, border='red', lwd=2, add=TRUE)
points(train, col='red', pch=20, cex=2)
points(test, col='black', pch=20, cex=2)

pr <- predict(ch, r, progress='')
plot(pr)
points(test, col='black', pch=20, cex=2)
points(train, col='red', pch=20, cex=2)

# to get the polygons:
p <- polygons(ch)
p



cleanEx()
nameEx("domain")
### * domain

flush(stderr()); flush(stdout())

### Name: domain
### Title: Domain
### Aliases: domain domain,Raster,SpatialPoints-method
###   domain,Raster,matrix-method domain,Raster,data.frame-method
###   domain,matrix,missing-method domain,data.frame,missing-method
###   Domain-class
### Keywords: spatial

### ** Examples

logo <- stack(system.file("external/rlogo.grd", package="raster"))
#presence data
pts <- matrix(c(48.243420, 48.243420, 47.985820, 52.880230, 49.531423, 46.182616, 54.168232, 
  69.624263, 83.792291, 85.337894, 74.261072, 83.792291, 95.126713, 84.565092, 66.275456, 
  41.803408, 25.832176, 3.936132, 18.876962, 17.331359,7.048974, 13.648543, 26.093446, 
  28.544714, 39.104026, 44.572240, 51.171810, 56.262906, 46.269272, 38.161230, 30.618865,
  21.945145, 34.390047, 59.656971, 69.839163, 73.233228, 63.239594, 45.892154, 43.252326,
  28.356155), ncol=2)
d <- domain(logo, pts)
p <- predict(d, logo)



cleanEx()
nameEx("ecocrop")
### * ecocrop

flush(stderr()); flush(stdout())

### Name: ecocrop
### Title: Ecocrop model
### Aliases: ecocrop getCrop ECOcrops ECOCROP-class ECOCROPcrop-class
### Keywords: spatial

### ** Examples

ecocrop('potato', 5:16, 15:26, runif(12)*100)
getCrop('Acacia brachystachya Benth.')
crop <- getCrop('Hot pepper')
ecocrop(crop, 5:16, 15:26, rainfed=FALSE)



cleanEx()
nameEx("ecolim")
### * ecolim

flush(stderr()); flush(stdout())

### Name: ecolim
### Title: Ecolim model
### Aliases: ecolim ecolim,matrix,matrix-method EcoLim-class
###   predict,EcoLim-method
### Keywords: spatial

### ** Examples

# get predictor variables
fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
              pattern='grd', full.names=TRUE ) 
env <- stack(fnames)

bio1 <- c(200,250,400,450)
bio12 <- c(0,1000, 3000, 4000)
r1 <- c(0, 1, 1, 0)
r2 <- c(0, 0, 1, 1)
x <- cbind(bio1, bio12)
y <- cbind(r1, r2)

e <- ecolim(x, y) 
plot(e, lwd=2, col='red')
p <- predict(e, env)
plot(p)

# no extrapolation:
ef <- ecolim(x, y, extrapolate=FALSE) 
pf <- predict(ef, env)
plot(pf)


occurence <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
occ <- read.table(occurence, header=TRUE, sep=',')[,-1]
fold <- kfold(occ, k=5)
occtest <- occ[fold == 1, ]
occtrain <- occ[fold != 1, ]
bg <- randomPoints(env, 1000)


## Not run: 
##D # An approach to optimize the values based on
##D # some known presences and (here random) absences
##D # for the same species as in the maxent example
##D 
##D # intial parameters
##D v <- c(200, 250, 400, 450, 0, 1000, 3000, 4000)
##D 
##D # function to be minimized
##D f <- function(p) {
##D 	x[] <- p
##D 	# numbers must go up
##D 	if ( any(x[-1,] < x[-nrow(x), ]) ) return(Inf)
##D 	e <- ecolim(x, y) 
##D 	# we are minimizing, hence 1-AUC
##D 	1-evaluate(e, p=occtrain, a=bg, x=env)@auc
##D }
##D 
##D # patience...
##D set.seed(0)
##D z <- optim(v, f)
##D 
##D x[] <- z$par
##D eco <- ecolim(x, y) 
##D evaluate(eco, p=occtest, a=bg, x=env)
##D 
##D set.seed(0)
##D pwd <- pwdSample(occtest,bg,occtrain)
##D ptest <- occtest[!is.na(pwd),]
##D atest <- bg[na.omit(pwd),]
##D evaluate(eco, p=ptest, a=atest, x=env)
##D 
##D p2 <- predict(eco, env)
##D plot(p2)
## End(Not run)



cleanEx()
nameEx("evaluate")
### * evaluate

flush(stderr()); flush(stdout())

### Name: evaluate
### Title: Model evaluation
### Aliases: evaluate
### Keywords: spatial

### ** Examples

## See ?maxent for an example with real data.
# this is a contrived example:
# p has the predicted values for 50 known cases (locations) 
# with presence of the phenomenon (species)
p <- rnorm(50, mean=0.7, sd=0.3)
# b has the predicted values for 50 background locations (or absence)
a <- rnorm(50, mean=0.4, sd=0.4)
e <- evaluate(p=p, a=a)

threshold(e)

plot(e, 'ROC')
plot(e, 'TPR')
boxplot(e)
density(e)

str(e)



cleanEx()
nameEx("gbif")
### * gbif

flush(stderr()); flush(stdout())

### Name: gbif
### Title: Data from GBIF
### Aliases: gbif
### Keywords: spatial

### ** Examples

## Not run: 
##D 
##D gbif('solanum', download=FALSE)
##D gbif('solanum', 'acaule', download=FALSE)
##D 
##D gbif('Batrachoseps', '' , down=FALSE)
##D gbif('Batrachoseps', 'luciae', down=FALSE)
##D g <- gbif('Batrachoseps', 'luciae', geo=TRUE)
##D plot(g$lon, g$lat)
##D 
##D gs <- gbif('Batrachoseps', 'luciae', sp=TRUE)
##D plot(gs)
## End(Not run)



cleanEx()
nameEx("gbm.step")
### * gbm.step

flush(stderr()); flush(stdout())

### Name: gbm.step
### Title: gbm step
### Aliases: gbm.step
### Keywords: spatial

### ** Examples

data(Anguilla_train)
# reduce data set to speed things up a bit
Anguilla_train = Anguilla_train[1:200,]
angaus.tc5.lr01 <- gbm.step(data=Anguilla_train, gbm.x = 3:14, gbm.y = 2, family = "bernoulli",
       tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5)



cleanEx()
nameEx("geoDist")
### * geoDist

flush(stderr()); flush(stdout())

### Name: Geographic Distance
### Title: Geographic distance model
### Aliases: geoDist geoDist,SpatialPoints-method geoDist,matrix-method
###   geoDist,data.frame-method GeographicDistance-class
### Keywords: spatial

### ** Examples

r <- raster(system.file("external/rlogo.grd", package="raster"))
#presence data
pts <- matrix(c(17, 42, 85, 70, 19, 53, 26, 84, 84, 46, 48, 85, 4, 95, 48, 54, 66, 74, 50, 48, 
        28, 73, 38, 56, 43, 29, 63, 22, 46, 45, 7, 60, 46, 34, 14, 51, 70, 31, 39, 26), ncol=2)
colnames(pts) <- c('x', 'y')

train <- pts[1:12, ]
test <- pts[13:20, ]
				 
gd <- geoDist(train, lonlat=FALSE)
p <- predict(gd, r)

## Not run: 
##D plot(p)
##D points(test, col='black', pch=20, cex=2)
##D points(train, col='red', pch=20, cex=2)
## End(Not run)



cleanEx()
nameEx("geoIDW")
### * geoIDW

flush(stderr()); flush(stdout())

### Name: InvDistW
### Title: Inverse-distance weighted model
### Aliases: geoIDW geoIDW,SpatialPoints,SpatialPoints-method
###   geoIDW,matrix,matrix-method geoIDW,data.frame,data.frame-method
###   InvDistWeightModel-class
### Keywords: spatial

### ** Examples

r <- raster(system.file("external/rlogo.grd", package="raster"))
# presence points
p <- matrix(c(17, 42, 85, 70, 19, 53, 26, 84, 84, 46, 48, 85, 4, 95, 48, 54, 66, 74, 50, 48, 
      28, 73, 38, 56, 43, 29, 63, 22, 46, 45, 7, 60, 46, 34, 14, 51, 70, 31, 39, 26), ncol=2)

# absence points
a <- matrix(c(30, 23, 5, 5, 31, 33, 91, 63, 60, 88, 93, 97, 65, 68, 85, 97, 35, 32, 29, 55,
      3, 8, 19, 71, 49, 36, 69, 41, 20, 28, 18, 9, 5, 9, 25, 71, 8, 32, 46, 60), ncol=2)

idw <- geoIDW(p, a)
prd <- predict(r, idw)

## Not run: 
##D plot(prd)
##D points(p)
##D points(a, pch='x')
## End(Not run)



cleanEx()
nameEx("geocode")
### * geocode

flush(stderr()); flush(stdout())

### Name: geocode
### Title: Georeferencing with Google
### Aliases: geocode
### Keywords: spatial

### ** Examples

## Not run: 
##D geocode(c('1600 Pennsylvania Ave NW, Washington DC', 'Luca, Italy', 'Kampala'))
##D geocode(c('San Jose', 'San Jose, Mexico'))
##D geocode(c('San Jose', 'San Jose, Mexico'), oneRecord=TRUE)
## End(Not run)



cleanEx()
nameEx("gmap")
### * gmap

flush(stderr()); flush(stdout())

### Name: gmap
### Title: Get a Google map
### Aliases: gmap Mercator
### Keywords: spatial

### ** Examples

## Not run: 
##D library(rgdal)
##D 
##D # get a map using names
##D g = gmap('Australia')
##D plot(g, inter=TRUE)
##D gs = gmap('Sydney, New South Wales, Australia', type='satellite')
##D plot(gs, inter=TRUE)
##D gs = gmap('Sydney, Australia', type='satellite', exp=3)
##D plot(gs, inter=TRUE)
##D gs = gmap('Sydney, Australia', type='hybrid', zoom=10, scale=2)
##D plot(gs, inter=TRUE)
##D 
##D # from a maxtrix with lon/lat points
##D x = runif(30)*10 + 40
##D y = runif(30)*10 - 20
##D xy = cbind(x, y)
##D g = gmap(xy, type='hybrid')
##D plot(g, inter=TRUE)
##D points(Mercator(xy) , col='red', pch=20)
##D 
##D # or from an Extent object
##D e = extent( -121.9531 , -120.3897 , 35.36 , 36.61956 )
##D # you can also get an Extent object by clicking on the map twice after using:
##D # drawExtent()
##D r = gmap(e)
##D plot(r, interpolate=TRUE)
##D 
##D # transform points to Mercator for plotting on top of map:
##D pt <- matrix(c(-121, 36), ncol=2)
##D ptm <- Mercator(pt)
##D points(ptm, cex=3, pch=20, col='blue')
##D Mercator(ptm, inverse=TRUE)
##D 
##D # transform Spatial objects to Mercator for plotting on top of map
##D # here for points, but particularly relevant for lines and polygons
##D pt <- data.frame(pt)
##D coordinates(pt) <- ~X1 + X2
##D proj4string(pt) <-"+proj=longlat +datum=WGS84 +ellps=WGS84"
##D ptm2 <- spTransform(pt, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 
##D       +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"))
##D points(ptm, col='red', pch='x', cex=3)
##D 
##D # styles:
##D g <- gmap("Brooklyn", style="feature:road.local|element:geometry|hue:0x00ff00|saturation:100
##D   &style=feature:landscape|element:geometry|lightness:-100", type='roadmap')
##D plot(g)
## End(Not run)



cleanEx()
nameEx("gridSample")
### * gridSample

flush(stderr()); flush(stdout())

### Name: gridSample
### Title: Stratified regular sample on a grid
### Aliases: gridSample
### Keywords: spatial

### ** Examples

x <- rnorm(1000, 10, 5)
y <- rnorm(1000, 50, 5)
xy <- cbind(x,y)
res <- 5
r <- raster(extent(range(xy[,1]), range(xy[,2])) + res)
res(r) <- res

samp <- gridSample(xy, r, n=1)
plot(xy, cex=0.1)
points(samp, pch='x', col='red')



cleanEx()
nameEx("kfold")
### * kfold

flush(stderr()); flush(stdout())

### Name: kfold
### Title: k-fold partitioning
### Aliases: kfold
### Keywords: spatial

### ** Examples


#library(disdat)
#data(NSWtrain)
## a single species
#srsp1 <- subset(NSWtrain, spid=='srsp1')
#kfold(srsp1, k = 5)

## all species
#k = kfold(NSWtrain, k=5, by=NSWtrain$spid)

#k[NSWtrain$spid=='srsp1']
## each group has the same number of records 
##(except for adjustments if the number of records divided by k is not an integer) 
#table(k[NSWtrain$spid=='srsp1'])
#k[NSWtrain$spid=='ousp5']




cleanEx()
nameEx("mahal")
### * mahal

flush(stderr()); flush(stdout())

### Name: mahal
### Title: Mahalanobis model
### Aliases: mahal mahal,Raster,SpatialPoints-method
###   mahal,Raster,matrix-method mahal,Raster,data.frame-method
###   mahal,matrix,missing-method mahal,data.frame,missing-method
###   Mahalanobis-class
### Keywords: spatial

### ** Examples

logo <- stack(system.file("external/rlogo.grd", package="raster"))

#presence data
pts <- matrix(c(48.243420, 48.243420, 47.985820, 52.880230, 49.531423, 46.182616, 
  54.168232, 69.624263, 83.792291, 85.337894, 74.261072, 83.792291, 95.126713, 
  84.565092, 66.275456, 41.803408, 25.832176, 3.936132, 18.876962, 17.331359, 
  7.048974, 13.648543, 26.093446, 28.544714, 39.104026, 44.572240, 51.171810, 
  56.262906, 46.269272, 38.161230, 30.618865, 21.945145, 34.390047, 59.656971, 
  69.839163, 73.233228, 63.239594, 45.892154, 43.252326, 28.356155), ncol=2)

# fit model
m <- mahal(logo, pts)

# make a prediction
predict(m, logo[1])

x <- predict(m, logo)

# or x <- predict(logo, m) via raster::predict

# plot(x > 0)



cleanEx()
nameEx("maxent")
### * maxent

flush(stderr()); flush(stdout())

### Name: maxent
### Title: Maxent
### Aliases: maxent maxent,missing,missing-method maxent,Raster,ANY-method
###   maxent,SpatialGridDataFrame,ANY-method
###   maxent,data.frame,vector-method MaxEnt-class MaxEntReplicates-class
### Keywords: spatial

### ** Examples

# only run if the maxent.jar file is available, in the right folder
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

# checking if maxent can be run (normally not part of your script)
if (file.exists(jar) & require(rJava)) {

# get predictor variables
fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
              pattern='grd', full.names=TRUE )
predictors <- stack(fnames)
#plot(predictors)

# file with presence points
occurence <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
occ <- read.table(occurence, header=TRUE, sep=',')[,-1]

# witholding a 20% sample for testing 
fold <- kfold(occ, k=5)
occtest <- occ[fold == 1, ]
occtrain <- occ[fold != 1, ]

# fit model, biome is a categorical variable
me <- maxent(predictors, occtrain, factors='biome')

# see the maxent results in a browser:
# me

# use "args"
# me2 <- maxent(predictors, occtrain, factors='biome', args=c("-J", "-P"))

# plot showing importance of each variable
plot(me)

# response curves
# response(me)

# predict to entire dataset
r <- predict(me, predictors) 

# with some options:
# r <- predict(me, predictors, args=c("outputformat=raw"), progress='text', 
#      filename='maxent_prediction.grd')

plot(r)
points(occ)

#testing
# background data
bg <- randomPoints(predictors, 1000)

#simplest way to use 'evaluate'
e1 <- evaluate(me, p=occtest, a=bg, x=predictors)

# alternative 1
# extract values
pvtest <- data.frame(extract(predictors, occtest))
avtest <- data.frame(extract(predictors, bg))

e2 <- evaluate(me, p=pvtest, a=avtest)

# alternative 2 
# predict to testing points 
testp <- predict(me, pvtest) 
head(testp)
testa <- predict(me, avtest) 

e3 <- evaluate(p=testp, a=testa)
e3
threshold(e3)

plot(e3, 'ROC')

}



cleanEx()
nameEx("mess")
### * mess

flush(stderr()); flush(stdout())

### Name: mess
### Title: Multivariate environmental similarity surfaces (MESS)
### Aliases: mess

### ** Examples


set.seed(9)
r <- raster(ncol=10, nrow=10)
r1 <- setValues(r, (1:ncell(r))/10 + rnorm(ncell(r)))
r2 <- setValues(r, (1:ncell(r))/10 + rnorm(ncell(r)))
r3 <- setValues(r, (1:ncell(r))/10 + rnorm(ncell(r)))
s <- stack(r1,r2,r3)
names(s) <- c('a', 'b', 'c')
xy <- cbind(rep(c(10,30,50), 3), rep(c(10,30,50), each=3))
refpt <- extract(s, xy)

ms <- mess(s, refpt, full=TRUE)
plot(ms)


## Not run: 
##D filename <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
##D bradypus <- read.table(filename, header=TRUE, sep=',')
##D bradypus <- bradypus[,2:3]
##D files <- list.files(path=paste(system.file(package="dismo"),'/ex', sep=''), 
##D    pattern='grd', full.names=TRUE )
##D predictors <- stack(files)
##D predictors <- dropLayer(x=predictors,i=9)
##D reference_points <- extract(predictors, bradypus)
##D mss <- mess(x=predictors, v=reference_points, full=TRUE)
##D plot(mss)
## End(Not run)




cleanEx()
nameEx("nicheOverlap")
### * nicheOverlap

flush(stderr()); flush(stdout())

### Name: nicheOverlap
### Title: Niche overlap
### Aliases: nicheOverlap

### ** Examples

r1 <- raster(nr=18, nc=36)
r2 <- raster(nr=18, nc=36)
set.seed(0)
r1[] <- runif(ncell(r1))
r2[] <- runif(ncell(r1))
nicheOverlap(r1, r2)



cleanEx()
nameEx("nullRandom")
### * nullRandom

flush(stderr()); flush(stdout())

### Name: Random null model
### Title: Random null model
### Aliases: nullRandom
### Keywords: spatial

### ** Examples

predictors <- stack(list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
               pattern='grd', full.names=TRUE ))
occurence <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
occ <- read.table(occurence, header=TRUE, sep=',')[,-1]
 
x <- extract(predictors, occ)
nr <- nullRandom(x, bioclim, n=25, rep=25, pa=FALSE)
mean(sapply(nr, function(x)x@auc))



cleanEx()
nameEx("plotEval")
### * plotEval

flush(stderr()); flush(stdout())

### Name: Evaluation plots
### Title: Plot model evaluation data
### Aliases: plot,ModelEvaluation,character-method
### Keywords: spatial

### ** Examples

# p = the predicted value for 50 known cases (locations) with presence of the phenomenon (species)
p = rnorm(50, mean=0.7, sd=0.3)
# b = the predicted value for 50 known cases (locations) with absence of the phenomenon (species)
a = rnorm(50, mean=0.4, sd=0.4)
e = evaluate(p=p, a=a)
plot(e, 'ROC')
plot(e, 'kappa')
plot(e, 'FPR')
plot(e, 'prevalence')



cleanEx()
nameEx("predict")
### * predict

flush(stderr()); flush(stdout())

### Name: predict
### Title: Distribution model predictions
### Aliases: predict predict,Bioclim-method predict,Domain-method
###   predict,Mahalanobis-method predict,MaxEnt-method
###   predict,MaxEntReplicates-method predict,ConvexHull-method
###   predict,CircleHull-method predict,RectangularHull-method
###   predict,CirclesRange-method predict,GeographicDistance-method
###   predict,InvDistWeightModel-method predict,VoronoiHull-method
### Keywords: methods spatial

### ** Examples

logo <- stack(system.file("external/rlogo.grd", package="raster"))
pts <- matrix(c(48, 48, 48, 53, 50, 46, 54, 70, 84, 85, 74, 84, 95, 85, 66, 
   42, 26, 4, 19, 17, 7, 14, 26, 29, 39, 45, 51, 56, 46, 38, 31, 22, 34, 60,
   70, 73, 63, 46, 43, 28), ncol=2)
b <- bioclim(logo, pts)
# prediction for a sub-region
e <- extent(30,90,20,60)
p <- predict(b, logo, progress='text', ext=e)
plot(p)



cleanEx()
nameEx("pwdSample")
### * pwdSample

flush(stderr()); flush(stdout())

### Name: pwdSample
### Title: Pair-wise distance sampling
### Aliases: pwdSample
### Keywords: spatial

### ** Examples

ref <- matrix(c(-54.5,-38.5, 2.5, -9.5, -45.5, 1.5, 9.5, 4.5, -10.5, -10.5), ncol=2)
fix <- matrix(c(-56.5, -30.5, -6.5, 14.5, -25.5, -48.5, 14.5, -2.5, 14.5,
               -11.5, -17.5, -11.5), ncol=2)
r <- raster()
extent(r) <- c(-110, 110, -45, 45)
r[] <- 1
set.seed(0)
sam <- randomPoints(r, n=50)

par(mfrow=c(1,2))
plot(sam, pch='x')
points(ref, col='red', pch=18, cex=2)
points(fix, col='blue', pch=20, cex=2)

i <- pwdSample(fix, sam, ref, lonlat=TRUE)
i
sfix <- fix[!is.na(i), ]
ssam <- sam[i[!is.na(i)], ]
ssam

plot(sam, pch='x', cex=0)
points(ssam, pch='x')
points(ref, col='red', pch=18, cex=2)
points(sfix, col='blue', pch=20, cex=2)

# try to get 3 pairs for each point in 'fixed'
pwdSample(fix, sam, ref, lonlat=TRUE, n=3)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("rectHull")
### * rectHull

flush(stderr()); flush(stdout())

### Name: rectHull
### Title: Rectangular hull model
### Aliases: rectHull rectHull,SpatialPoints-method rectHull,matrix-method
###   rectHull,data.frame-method RectangularHull-class
### Keywords: spatial

### ** Examples

r <- raster(system.file("external/rlogo.grd", package="raster"))
# presence data
pts <- matrix(c(17, 42, 85, 70, 19, 53, 26, 84, 84, 46, 48, 85, 4, 95, 48, 54, 66, 
 74, 50, 48, 28, 73, 38, 56, 43, 29, 63, 22, 46, 45, 7, 60, 46, 34, 14, 51, 70, 31, 39, 26), ncol=2)
train <- pts[1:12, ]
test <- pts[13:20, ]
				 
rh <- rectHull(train)
predict(rh, test)

plot(r)
plot(rh, border='red', lwd=2, add=TRUE)
points(train, col='red', pch=20, cex=2)
points(test, col='black', pch=20, cex=2)

pr <- predict(rh, r, progress='')
plot(pr)
points(test, col='black', pch=20, cex=2)
points(train, col='red', pch=20, cex=2)



cleanEx()
nameEx("ssb")
### * ssb

flush(stderr()); flush(stdout())

### Name: ssb
### Title: Spatial sorting bias
### Aliases: ssb
### Keywords: spatial

### ** Examples

ref <- matrix(c(-54.5,-38.5, 2.5, -9.5, -45.5, 1.5, 9.5, 4.5, -10.5, -10.5), ncol=2)
p <- matrix(c(-56.5, -30.5, -6.5, 14.5, -25.5, -48.5, 14.5, -2.5, 14.5, 
        -11.5, -17.5, -11.5), ncol=2)
r <- raster()
extent(r) <- c(-110, 110, -45, 45)
r[] <- 1
set.seed(0)
a <- randomPoints(r, n=50)
b <- ssb(p, a, ref)

# distances in km
b / 1000

# an index of spatial sorting bias (1 is no bias, near 0 is extreme bias)
b[1] / b[2]



cleanEx()
nameEx("threshold")
### * threshold

flush(stderr()); flush(stdout())

### Name: threshold
### Title: Find a threshold
### Aliases: threshold threshold,ModelEvaluation-method
### Keywords: spatial

### ** Examples

## See ?maxent for an example with real data.
# this is a contrived example:
# p has the predicted values for 50 known cases (locations)
# with presence of the phenomenon (species)
p <- rnorm(50, mean=0.7, sd=0.3)
# b has the predicted values for 50 background locations (or absence)
a <- rnorm(50, mean=0.4, sd=0.4)
e <- evaluate(p=p, a=a)

threshold(e)



cleanEx()
nameEx("voronoi")
### * voronoi

flush(stderr()); flush(stdout())

### Name: voronoi
### Title: Voronoi polygons
### Aliases: voronoi
### Keywords: spatial

### ** Examples

# points
p <- matrix(c(17, 42, 85, 70, 19, 53, 26, 84, 84, 46, 48, 85, 4, 95, 48, 54, 66, 74, 50, 48, 
      28, 73, 38, 56, 43, 29, 63, 22, 46, 45, 7, 60, 46, 34, 14, 51, 70, 31, 39, 26), ncol=2)
	  
v <- voronoi(p)
v



cleanEx()
nameEx("voronoiHulll")
### * voronoiHulll

flush(stderr()); flush(stdout())

### Name: Voronoi Hull
### Title: Voronoi hull model
### Aliases: voronoiHull voronoiHull,SpatialPoints,SpatialPoints-method
###   voronoiHull,matrix,matrix-method
###   voronoiHull,data.frame,data.frame-method VoronoiHull-class
### Keywords: spatial

### ** Examples

r <- raster(system.file("external/rlogo.grd", package="raster"))
# presence points
p <- matrix(c(17, 42, 85, 70, 19, 53, 26, 84, 84, 46, 48, 85, 4, 95, 48, 54, 66, 74, 50, 48, 
      28, 73, 38, 56, 43, 29, 63, 22, 46, 45, 7, 60, 46, 34, 14, 51, 70, 31, 39, 26), ncol=2)

# absence points
a <- matrix(c(30, 23, 5, 5, 31, 33, 91, 63, 60, 88, 93, 97, 65, 68, 85, 97, 35, 32, 29, 55,
      3, 8, 19, 71, 49, 36, 69, 41, 20, 28, 18, 9, 5, 9, 25, 71, 8, 32, 46, 60), ncol=2)

v <- voronoiHull(p, a)
	  
x <- predict(r, v)

## Not run: 
##D plot(x)
##D points(p, col='black', pch=20, cex=2)
##D points(a, col='red', pch=20, cex=2)
## End(Not run)



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
