### R code from vignette source 'sdm.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: foo
###################################################
options(width = 60)
set.seed(0)


###################################################
### code chunk number 2: sdm1
###################################################
# loads the dismo library
library(dismo)
file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
# this is the file we will use:
file

# read it
bradypus <- read.table(file,  header=TRUE,  sep=",")
# inspect the values of the file
# first rows
head(bradypus)
# we only need columns 2 and 3:
bradypus <- bradypus[,2:3]
head(bradypus)


###################################################
### code chunk number 3: sdm2
###################################################
# load the saved S. acaule data
data(acaule)

# how many rows and colums?
dim(acaule)

#select the records that have longitude and latitude data
colnames(acaule)
acgeo <- subset(acaule, !is.na(lon) & !is.na(lat))
dim(acgeo)

# show some values
acgeo[1:4, c(1:5,7:10)]


###################################################
### code chunk number 4: sdm3
###################################################
library(maptools)
data(wrld_simpl)
plot(wrld_simpl, xlim=c(-80,70), ylim=c(-60,10), axes=TRUE, 
        col="light yellow")

# restore the box around the map
box()

# plot points
points(acgeo$lon, acgeo$lat, col='orange', pch=20, cex=0.75)
# plot points again to add a border, for better visibility 
points(acgeo$lon, acgeo$lat, col='red', cex=0.75)


###################################################
### code chunk number 5: sdm4a
###################################################
acaule[c(303,885),1:10]


###################################################
### code chunk number 6: sdm4b
###################################################
lonzero = subset(acgeo, lon==0)
# show all records, only the first 13 columns
lonzero[, 1:13]


###################################################
### code chunk number 7: sdm5a
###################################################
# which records are duplicates (only for the first 10 columns)?
dups <- duplicated(lonzero[, 1:10])
# remove duplicates
lonzero  <-  lonzero[dups, ]
lonzero[,1:13]


###################################################
### code chunk number 8: sdm5b
###################################################
# differentiating by (sub) species
# dups2 <- duplicated(acgeo[, c('species', 'lon', 'lat')])
# ignoring (sub) species and other naming variation
dups2 <- duplicated(acgeo[, c('lon', 'lat')])
# number of duplicates
sum(dups2)
# keep the records that are _not_ duplicated
acg <- acgeo[!dups2, ]


###################################################
### code chunk number 9: sdm5c
###################################################
i <- acg$lon > 0 & acg$lat > 0
acg$lon[i] <- -1 * acg$lon[i]
acg$lat[i] <- -1 * acg$lat[i]
acg <- acg[acg$lon < -50 & acg$lat > -50, ]


###################################################
### code chunk number 10: sdm6a
###################################################
library(sp)
coordinates(acg) <- ~lon+lat
crs(acg) <- crs(wrld_simpl)
class(acg)


###################################################
### code chunk number 11: sdm6b
###################################################
class(wrld_simpl)
ovr <- over(acg, wrld_simpl)


###################################################
### code chunk number 12: sdm6c
###################################################
head(ovr)
cntr <- ovr$NAME


###################################################
### code chunk number 13: sdm6d
###################################################
i <- which(is.na(cntr))
i
j <- which(cntr != acg$country)
# for the mismatches, bind the country names of the polygons and points
cbind(cntr, acg$country)[j,]


###################################################
### code chunk number 14: sdm6d
###################################################
plot(acg)
plot(wrld_simpl, add=T, border='blue', lwd=2)
points(acg[j, ], col='red', pch=20, cex=2)


###################################################
### code chunk number 15: sdm8
###################################################
georef <- subset(acaule, (is.na(lon) | is.na(lat)) & ! is.na(locality) )
dim(georef)
georef[1:3,1:13]


###################################################
### code chunk number 16: sdm9
###################################################
georef$cloc[4]
b <- try(  geocode(georef$cloc[4]) )
b


###################################################
### code chunk number 17: sdm10
###################################################
# create a RasterLayer with the extent of acgeo
r <- raster(acg)
# set the resolution of the cells to (for example) 1 degree
res(r) <- 1

# expand (extend) the extent of the RasterLayer a little
r <- extend(r, extent(r)+1)

# sample:
acsel <- gridSample(acg, r, n=1) 

# to illustrate the method and show the result
p <- rasterToPolygons(r)
plot(p, border='gray')
points(acg)
# selected points in red
points(acsel, cex=1, col='red', pch='x')


###################################################
### code chunk number 18: sdm12
###################################################
file <- paste(system.file(package="dismo"), '/ex/acaule.csv', sep='')
acsel <- read.csv(file)


###################################################
### code chunk number 19: sdm15a
###################################################
# get the file names 
files <- list.files(path=paste(system.file(package="dismo"), '/ex', 
                       sep=''),  pattern='grd',  full.names=TRUE )

# we use the first file to create a RasterLayer
mask <- raster(files[1])

# select 500 random points
# set seed to assure that the examples will always
# have the same random sample.
set.seed(1963)
bg <- randomPoints(mask, 500 )


###################################################
### code chunk number 20: sdm15
###################################################
# set up the plotting area for two maps
par(mfrow=c(1,2))
plot(!is.na(mask), legend=FALSE)
points(bg, cex=0.5)

# now we repeat the sampling, but limit 
# the area of sampling using a spatial extent
e <- extent(-80, -53, -39, -22)
bg2 <- randomPoints(mask, 50, ext=e)
plot(!is.na(mask), legend=FALSE)
plot(e, add=TRUE, col='red')
points(bg2, cex=0.5)


###################################################
### code chunk number 21: sdm16a
###################################################
file <- paste(system.file(package="dismo"), '/ex/acaule.csv', sep='')
ac <- read.csv(file)


###################################################
### code chunk number 22: sdm16b
###################################################
coordinates(ac) <- ~lon+lat
projection(ac) <- CRS('+proj=longlat +datum=WGS84')


###################################################
### code chunk number 23: sdm17
###################################################
# circles with a radius of 50 km
x <- circles(ac, d=50000, lonlat=TRUE)
pol <- polygons(x)


###################################################
### code chunk number 24: sdm19
###################################################
# sample randomly from all circles
samp1 <- spsample(pol, 250, type='random', iter=25)
# get unique cells
cells <- cellFromXY(mask, samp1)
length(cells)
cells <- unique(cells)
length(cells)
xy <- xyFromCell(mask, cells)


###################################################
### code chunk number 25: sdm20
###################################################
plot(pol, axes=TRUE)
points(xy, cex=0.75, pch=20, col='blue')


###################################################
### code chunk number 26: sdm21a
###################################################
spxy <- SpatialPoints(xy, proj4string=CRS('+proj=longlat +datum=WGS84'))
o <- over(spxy, geometry(x))
xyInside <- xy[!is.na(o), ]


###################################################
### code chunk number 27: sdm21b
###################################################
# extract cell numbers for the circles
v <- extract(mask, x@polygons, cellnumbers=T)
# use rbind to combine the elements in list v
v <- do.call(rbind, v)

# get unique cell numbers from which you could sample
v <- unique(v[,1])
head(v)

# to display the results
m <- mask
m[] <- NA
m[v] <- 1
plot(m, ext=extent(x@polygons)+1)
plot(x@polygons, add=T)


###################################################
### code chunk number 28: sdm22
###################################################
files <- list.files(path=paste(system.file(package="dismo"), 
              '/ex', sep=''), pattern='grd', full.names=TRUE )
# The above finds all the files with extension "grd" in the
# examples ("ex") directory of the dismo package. You do not
# need such a complex statement to get your own files.
files
predictors <- stack(files)
predictors
names(predictors)
plot(predictors)


###################################################
### code chunk number 29: sdm23a
###################################################
library(maptools)
data(wrld_simpl)
file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
bradypus <- read.table(file,  header=TRUE,  sep=',')
# we do not need the first column
bradypus  <- bradypus[,-1]


###################################################
### code chunk number 30: sdm23b
###################################################
# first layer of the RasterStack
plot(predictors, 1)

# note the "add=TRUE" argument with plot
plot(wrld_simpl, add=TRUE)
# with the points function, "add" is implicit
points(bradypus, col='blue')


###################################################
### code chunk number 31: sdm24a
###################################################
presvals <- extract(predictors, bradypus)
# setting random seed to always create the same
# random set of points for this example
set.seed(0)
backgr <- randomPoints(predictors, 500)
absvals <- extract(predictors, backgr)
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))
sdmdata[,'biome'] = as.factor(sdmdata[,'biome'])
head(sdmdata)
tail(sdmdata)
summary(sdmdata)


###################################################
### code chunk number 32: sdm24b
###################################################
# pairs plot of the values of the climate data 
# at the bradypus occurrence sites.
pairs(sdmdata[,2:5], cex=0.1, fig=TRUE)


###################################################
### code chunk number 33: sdm25
###################################################
m1 <- glm(pb ~ bio1 + bio5 + bio12, data=sdmdata)
class(m1)
summary(m1)

m2 = glm(pb ~ ., data=sdmdata)
m2


###################################################
### code chunk number 34: sdm26
###################################################
bc <- bioclim(presvals[,c('bio1', 'bio5', 'bio12')])
class(bc)
bc
pairs(bc)


###################################################
### code chunk number 35: sdm27a
###################################################
bio1 = c(40, 150, 200)
bio5 = c(60, 115, 290)
bio12 = c(600, 1600, 1700)
pd = data.frame(cbind(bio1, bio5, bio12))
pd
predict(m1, pd)
predict(bc, pd)


###################################################
### code chunk number 36: sdm27b
###################################################
response(bc)


###################################################
### code chunk number 37: sdm27c
###################################################
names(predictors)
p <- predict(predictors, m1)
plot(p)


###################################################
### code chunk number 38: sdm28
###################################################
p <- rnorm(50, mean=0.7, sd=0.3)
a <- rnorm(50, mean=0.4, sd=0.4)
par(mfrow=c(1, 2))
plot(sort(p), col='red', pch=21)
points(sort(a), col='blue', pch=24)
legend(1, 0.95 * max(a,p), c('presence', 'absence'),
          pch=c(21,24), col=c('red', 'blue'))
comb = c(p,a)
group = c(rep('presence', length(p)), rep('absence', length(a)))
boxplot(comb~group, col=c('blue', 'red'))


###################################################
### code chunk number 39: sdm29
###################################################
group = c(rep(1, length(p)), rep(0, length(a))) 
cor.test(comb, group)$estimate
mv <- wilcox.test(p,a)
auc <- as.numeric(mv$statistic) / (length(p) * length(a))
auc


###################################################
### code chunk number 40: sdm40
###################################################
e <- evaluate(p=p, a=a)
class(e)
e
par(mfrow=c(1, 2))
density(e)
boxplot(e, col=c('blue', 'red'))


###################################################
### code chunk number 41: sdm41
###################################################
samp <- sample(nrow(sdmdata), round(0.75 * nrow(sdmdata)))
traindata <- sdmdata[samp,]
traindata <- traindata[traindata[,1] == 1, 2:9]
testdata <- sdmdata[-samp,]
bc <- bioclim(traindata)
e <- evaluate(testdata[testdata==1,], testdata[testdata==0,], bc)
e
plot(e, 'ROC')


###################################################
### code chunk number 42: sdm42
###################################################
pres <- sdmdata[sdmdata[,1] == 1, 2:9]
back <- sdmdata[sdmdata[,1] == 0, 2:9]


###################################################
### code chunk number 43: sdm43
###################################################
k <- 5
group <- kfold(pres, k)
group[1:10]
unique(group)


###################################################
### code chunk number 44: sdm44
###################################################
e <- list()
for (i in 1:k) {
	train <- pres[group != i,]
	test <- pres[group == i,]
	bc <- bioclim(train)
	e[[i]] <- evaluate(p=test, a=back, bc)
}	


###################################################
### code chunk number 45: sdm45a
###################################################
auc <- sapply( e, function(x){slot(x, 'auc')} )
auc
mean(auc)

 sapply( e, function(x){ x@t[which.max(x@TPR + x@TNR)] } )
# equivalent  
# sapply( e, function(x){ x@t[which.max(x@TPR + x@TNR)] } )


###################################################
### code chunk number 46: sdm45b
###################################################
nr <- nrow(bradypus)
s <- sample(nr, 0.25 * nr)
pres_train <- bradypus[-s, ]
pres_test <- bradypus[s, ]

nr <- nrow(backgr)
s <- sample(nr, 0.25 * nr)
back_train <- backgr[-s, ]
back_test <- backgr[s, ]


###################################################
### code chunk number 47: sdm45b
###################################################
sb <- ssb(pres_test, back_test, pres_train)
sb[,1] / sb[,2]


###################################################
### code chunk number 48: sdm45c
###################################################
i <- pwdSample(pres_test, back_test, pres_train, n=1, tr=0.1)
pres_test_pwd <- pres_test[!is.na(i[,1]), ]
back_test_pwd <- back_test[na.omit(as.vector(i)), ]
sb2 <- ssb(pres_test_pwd, back_test_pwd, pres_train)
sb2[1]/ sb2[2]


###################################################
### code chunk number 49: sdm45d
###################################################
bc <- bioclim(predictors, pres_train)
evaluate(bc, p=pres_test, a=back_test, x=predictors)
evaluate(bc, p=pres_test_pwd, a=back_test_pwd, x=predictors)


###################################################
### code chunk number 50: sdm46a
###################################################
files <- list.files(path=paste(system.file(package="dismo"), 
              '/ex', sep=''), pattern='grd', full.names=TRUE )
predictors <- stack(files)

file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
bradypus <- read.table(file,  header=TRUE,  sep=',')
bradypus  <- bradypus[,-1]
presvals <- extract(predictors, bradypus)
set.seed(0)
backgr <- randomPoints(predictors, 500)
absvals <- extract(predictors, backgr)
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))
sdmdata[,'biome'] = as.factor(sdmdata[,'biome'])


###################################################
### code chunk number 51: sdm46b
###################################################
pred_nf <- dropLayer(predictors, 'biome')


###################################################
### code chunk number 52: sdm47
###################################################
group <- kfold(bradypus, 5)
pres_train <- bradypus[group != 1, ]
pres_test <- bradypus[group == 1, ]


###################################################
### code chunk number 53: sdm48
###################################################
ext = extent(-90, -32, -33, 23)


###################################################
### code chunk number 54: sdm49
###################################################
backg <- randomPoints(pred_nf, n=1000, ext=ext, extf = 1.25)
colnames(backg) = c('lon', 'lat')
group <- kfold(backg, 5)
backg_train <- backg[group != 1, ]
backg_test <- backg[group == 1, ]


###################################################
### code chunk number 55: sdm50
###################################################
r = raster(pred_nf, 1)
plot(!is.na(r), col=c('white', 'light grey'), legend=FALSE)
plot(ext, add=TRUE, col='red', lwd=2)
points(backg_train, pch='-', cex=0.5, col='yellow')
points(backg_test, pch='-',  cex=0.5, col='black')
points(pres_train, pch= '+', col='green')
points(pres_test, pch='+', col='blue')


###################################################
### code chunk number 56: sdm60
###################################################
bc <- bioclim(pred_nf, pres_train)
plot(bc, a=1, b=2, p=0.85)


###################################################
### code chunk number 57: sdm61a
###################################################
e <- evaluate(pres_test, backg_test, bc, pred_nf)
e


###################################################
### code chunk number 58: sdm61b
###################################################
tr <- threshold(e, 'spec_sens')
tr


###################################################
### code chunk number 59: sdm62
###################################################
pb <- predict(pred_nf, bc, ext=ext, progress='')
pb
par(mfrow=c(1,2))
plot(pb, main='Bioclim, raw values')
plot(wrld_simpl, add=TRUE, border='dark grey')
plot(pb > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')


###################################################
### code chunk number 60: sdm63
###################################################
dm <- domain(pred_nf, pres_train)
e <- evaluate(pres_test, backg_test, dm, pred_nf)
e
pd = predict(pred_nf, dm, ext=ext, progress='')
par(mfrow=c(1,2))
plot(pd, main='Domain, raw values')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr <- threshold(e, 'spec_sens')
plot(pd > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')


###################################################
### code chunk number 61: sdm64
###################################################
mm <- mahal(pred_nf, pres_train)
e <- evaluate(pres_test, backg_test, mm, pred_nf)
e
pm = predict(pred_nf, mm, ext=ext, progress='')
par(mfrow=c(1,2))
pm[pm < -10] <- -10
plot(pm, main='Mahalanobis distance')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr <- threshold(e, 'spec_sens')
plot(pm > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')


###################################################
### code chunk number 62: sdm65
###################################################
train <- rbind(pres_train, backg_train)
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
envtrain <- extract(predictors, train)
envtrain <- data.frame( cbind(pa=pb_train, envtrain) )
envtrain[,'biome'] = factor(envtrain[,'biome'], levels=1:14)
head(envtrain)

testpres <- data.frame( extract(predictors, pres_test) )
testbackg <- data.frame( extract(predictors, backg_test) )
testpres[ ,'biome'] = factor(testpres[ ,'biome'], levels=1:14)
testbackg[ ,'biome'] = factor(testbackg[ ,'biome'], levels=1:14)


###################################################
### code chunk number 63: sdm66
###################################################
# logistic regression:
gm1 <- glm(pa ~ bio1 + bio5 + bio6 + bio7 + bio8 + bio12 + bio16 + bio17, 
            family = binomial(link = "logit"), data=envtrain)

summary(gm1)
coef(gm1)

gm2 <- glm(pa ~ bio1+bio5 + bio6 + bio7 + bio8 + bio12 + bio16 + bio17,
            family = gaussian(link = "identity"), data=envtrain)
			
evaluate(testpres, testbackg, gm1)
ge2 <- evaluate(testpres, testbackg, gm2)
ge2


###################################################
### code chunk number 64: sdm67
###################################################
pg <- predict(predictors, gm2, ext=ext)
par(mfrow=c(1,2))
plot(pg, main='GLM/gaussian, raw values')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr <- threshold(ge2, 'spec_sens')
plot(pg > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='-', cex=0.25)


###################################################
### code chunk number 65: sdm68a
###################################################
# checking if the jar file is present. If not, skip this bit
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if (file.exists(jar)) {
	xm <- maxent(predictors, pres_train, factors='biome')
	plot(xm)
} else {
    cat('cannot run this example because maxent is not available')
	plot(1)
}


###################################################
### code chunk number 66: sdm68b
###################################################
if (file.exists(jar)) {
	response(xm)
} else {
    cat('cannot run this example because maxent is not available')
	plot(1)
}


###################################################
### code chunk number 67: sdm69
###################################################
if (file.exists(jar)) {
	e <- evaluate(pres_test, backg_test, xm, predictors)
	e
	px <- predict(predictors, xm, ext=ext, progress='')
	par(mfrow=c(1,2))
	plot(px, main='Maxent, raw values')
	plot(wrld_simpl, add=TRUE, border='dark grey')
	tr <- threshold(e, 'spec_sens')
	plot(px > tr, main='presence/absence')
	plot(wrld_simpl, add=TRUE, border='dark grey')
	points(pres_train, pch='+')
} else {	
	plot(1)
}


###################################################
### code chunk number 68: sdm80
###################################################
library(randomForest)
model <- pa ~ bio1 + bio5 + bio6 + bio7 + bio8 + bio12 + bio16 + bio17
rf1 <- randomForest(model, data=envtrain)
model <- factor(pa) ~ bio1 + bio5 + bio6 + bio7 + bio8 + bio12 + bio16 + bio17
rf2 <- randomForest(model, data=envtrain)
rf3 <- randomForest(envtrain[,1:8], factor(pb_train))

erf <- evaluate(testpres, testbackg, rf1)
erf

pr <- predict(predictors, rf1, ext=ext)

par(mfrow=c(1,2))
plot(pr, main='Random Forest, regression')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr <- threshold(erf, 'spec_sens')
plot(pr > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='-', cex=0.25)


###################################################
### code chunk number 69: sdm81
###################################################
library(kernlab)
svm <- ksvm(pa ~ bio1+bio5+bio6+bio7+bio8+bio12+bio16+bio17, data=envtrain)
esv <- evaluate(testpres, testbackg, svm)
esv
ps <- predict(predictors, svm, ext=ext)

par(mfrow=c(1,2))
plot(ps, main='Support Vector Machine')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr <- threshold(esv, 'spec_sens')
plot(ps > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='-', cex=0.25)


###################################################
### code chunk number 70: sdm82
###################################################
models <- stack(pb, pd, pm, pg, pr, ps)
names(models) <- c("bioclim", "domain", "mahal", "glm", "rf", "svm")
plot(models)


###################################################
### code chunk number 71: sdm83
###################################################
m <- mean(models)
plot(m, main='average score')


###################################################
### code chunk number 72: sdm84
###################################################
auc <- sapply(list(ge2, erf, esv), function(x) x@auc)
w <- (auc-0.5)^2
m2 <- weighted.mean( models[[c("glm", "rf", "svm")]], w)
plot(m2, main='weighted mean of three models')


###################################################
### code chunk number 73: sdm100
###################################################
# first create a mask to predict to, and to use as a mask 
# to only predict to land areas
seamask <- crop(predictors[[1]], ext)

distm <- geoDist(pres_train, lonlat=TRUE)
ds <- predict(seamask, distm, mask=TRUE)

e <- evaluate(distm, p=pres_test, a=backg_test)
e

par(mfrow=c(1,2))
plot(ds, main='Geographic Distance')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr <- threshold(e, 'spec_sens')
plot(ds > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='-', cex=0.25)


###################################################
### code chunk number 74: sdm102
###################################################
hull <- convHull(pres_train, lonlat=TRUE)
e <- evaluate(hull, p=pres_test, a=backg_test)
e

h <- predict(seamask, hull, mask=TRUE)

plot(h, main='Convex Hull')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='-', cex=0.25)


###################################################
### code chunk number 75: sdm104
###################################################
circ <- circles(pres_train, lonlat=TRUE)
pc <- predict(seamask, circ, mask=TRUE)

e <- evaluate(circ, p=pres_test, a=backg_test)
e

par(mfrow=c(1,2))
plot(pc, main='Circles')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr <- threshold(e, 'spec_sens')
plot(pc > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='-', cex=0.25)


###################################################
### code chunk number 76: sdm106
###################################################
idwm <- geoIDW(p=pres_train, a=data.frame(back_train))

e <- evaluate(idwm, p=pres_test, a=backg_test)
e

iw <- predict(seamask, idwm, mask=TRUE)

par(mfrow=c(1,2))

plot(iw, main='Inv. Dist. Weighted')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr <- threshold(e, 'spec_sens')
pa <- mask(iw > tr, seamask)
plot(pa, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='-', cex=0.25)


###################################################
### code chunk number 77: sdm108
###################################################
# take a smallish sample of the background training data
va <- data.frame(back_train[sample(nrow(back_train), 100), ])
vorm <- voronoiHull(p=pres_train, a=va)

e <- evaluate(vorm, p=pres_test, a=backg_test)
e

vo <- predict(seamask, vorm, mask=T)

plot(vo, main='Voronoi Hull')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='-', cex=0.25)


