a <- 9
expect_that(a, is_less_than(10))


# test example of null model

envi <-raster::stack(system.file('extdata', 'pred.grd',package='SDMPlay'))

# Visited sites: for the example, we will take random points in the area
layer1 <- raster::subset(envi, 1)
# extract the coordinates of 500 random pixels in the area
visit <- dismo::randomPoints(layer1,n=500)
expect_equal(nrow(visit), 500)

xydata <- visit
expect_equal(nrow(xydata), 500)
predictors <- envi
background.nb <- 300
unique.data <- T
same <- T

# cleaning the xydata file and removing the data that are out of the extent of the RasterStack
xydata <- base::subset(xydata, xydata[, 1] >= predictors@extent@xmin & xydata[, 1] <= predictors@extent@xmax)
xydata <- base::subset(xydata, xydata[, 2] >= predictors@extent@ymin & xydata[, 2] <= predictors@extent@ymax)
print(paste("xydata=",nrow(xydata)))

# extracting the environmental variables associated with the presence data
presvals <- raster::extract(predictors, xydata)
expect_equal(nrow(presvals), 500)
presvals.cellnb <- raster::extract(predictors, xydata, cellnumbers = T)
expect_equal(nrow(presvals.cellnb), 500)

# remove the double occurrences
if (unique.data == TRUE) {
  presvals <- unique(presvals)
  presvals.cellnb <- unique(presvals.cellnb[,-1])
}
print(paste("presvals=",nrow(presvals)))
print(paste("presvals.cellnb=",nrow(presvals.cellnb)))
print(head(presvals))
print(head(presvals.cellnb))

#
# ####################################### try other data sets ##################################
# library(ncdf4)
# oxy_fond <-raster::raster("inst/supp_try/oxy_fond.nc")
# slope <-raster::raster("inst/supp_try/pente.nc")
# depth <-raster::raster("inst/supp_try/profondeur.nc")
# sali <-raster::raster("inst/supp_try/sali_fond_0512.nc")
#
# depth <- raster::crop(depth, extent(extent(sali)))
# raster::origin(depth) <- 0
# raster::origin(slope) <- 0
# raster::origin(sali) <- 0
# raster::origin(oxy_fond) <- 0
#
# # common extent
# mini <- depth+sali+oxy_fond+slope
#
# raster::extent(slope) <- raster::extent(mini)
# raster::extent(oxy_fond) <- raster::extent(mini)
# raster::extent(depth) <- raster::extent(mini)
# raster::extent(sali) <- raster::extent(mini)
#
# # create a stack
# pred.stack <- raster:: stack(depth,oxy_fond,slope,sali)
# names(pred.stack)<- c("a","b","c","d")
#
# # test delim.area
# stack2 <- SDMPlay:::delim.area(predictors = pred.stack, longmin=60, longmax = 80, latmin = -55, latmax = -45, interval=c(-2000,0))
# extent(stack2)
# stack2
# plot(stack2)
#
# # if the first rasterlayer is not depth
# pred.stack3 <- raster:: stack(oxy_fond,depth,slope,sali)
#
# stack3 <- SDMPlay:::delim.area(predictors = pred.stack, longmin=60, longmax = 80, latmin = -55, latmax = -45, interval=c(4,6))
# plot(stack3)
#
# # test SDMtab
# # open occurrences data
# occ <- read.csv("inst/supp_try/dataKer.csv",header=T,sep=";")
# occ <- data.frame(occ)
# occ <- occ[,c(2,16)]
# tab <- SDMPlay:::SDMtab(xydata=occ,predictors=stack2)
# head(tab)
#
# xydata <- occ
# predictors <- stack2
# unique.data <- T
# same <-T
#
# # cleaning the xydata file and removing the data that are out of the extent of the RasterStack
# xydata <- base::subset(xydata, xydata[, 1] >= predictors@extent@xmin & xydata[, 1] <= predictors@extent@xmax)
# xydata <- base::subset(xydata, xydata[, 2] >= predictors@extent@ymin & xydata[, 2] <= predictors@extent@ymax)
#
# # extracting the environmental variables associated with the presence data
# presvals <- raster::extract(predictors, xydata)
# presvals.cellnb <- raster::extract(predictors, xydata, cellnumbers = T)
#
# # remove the double occurrences
# if (unique.data == TRUE) {
#   presvals <- unique(presvals)
#   presvals.cellnb <- presvals.cellnb[!base::duplicated(presvals.cellnb[,-1]),]
# }
#
# # extract background locations (random sampling of latitude and longitudes in the area)
# if (same == TRUE) {
#   background.number <- dismo::randomPoints(predictors, n = nrow(presvals))
# } else {
#   background.number <- dismo::randomPoints(predictors, n = background.nb)
# }
# colnames(background.number) <- colnames(xydata)
#
# # extract environmental values at these locations
# pseudoabs.vals <- raster::extract(predictors, background.number)
# pseudoabs.vals.cellnb <- raster::extract(predictors, background.number, cellnumbers = T)
#
# # build the id column if
# id <- c(rep(1, nrow(presvals)), rep(0, nrow(pseudoabs.vals)))
#
# # build the final table latitudes and longitudes of the presence and background data
# xypres <- raster::xyFromCell(raster::subset(predictors, 4), presvals.cellnb[, 1]-10)
# xyback <- raster::xyFromCell(raster::subset(predictors, 1), pseudoabs.vals.cellnb[, 1])
# coord <- base::rbind(xypres, xyback)
#
#
# inter <- base::rbind(presvals, pseudoabs.vals)
# SDMtab.object <- base::data.frame(base::cbind(id, coord,inter ),check.rows=T)
# colnames(SDMtab.object) <- c("id", "latitude", "longitude", colnames(presvals))
#
#
#
# tab <- SDMPlay:::SDMtab(xydata=occ,predictors=stack2)
# nrow(tab)
# # choose a different number of background data
# tab2 <- SDMPlay:::SDMtab(xydata=occ,predictors=stack2,background.nb = NA)
# nrow(tab2)
#
# # test models functions
# model.maxent <- SDMPlay:::compute.maxent(x=tab,proj.predictors = predictors)
# model.brt <- SDMPlay:::compute.brt(x = tab,proj.predictors = predictors)
# model.brt$algorithm
# model.brt$response
# plot(model.brt$response)
# plot(model.brt$raster.prediction)
#
# # test SDMeval
# evalua <- SDMPlay:::SDMeval(model.brt)
# head(evalua)
# evalua2 <- SDMPlay:::SDMeval(model.maxent )
# head(evalua2)
#
# # test SDMdata.quality
# qual <- SDMPlay:::SDMdata.quality(tab)
# qual
#
# data <- tab
#
# table <- base::colSums(is.na(data[,-c(1:3)]))
# head(table)
#
# table <- base::data.frame(table)
# # calculate the proportion of NA values
# prop <- (table/base::nrow(data)) * 100
# base::colnames(prop) <- "NA.percent (%)"
#
# # test null model
# # map visited places
# map <- read.csv("inst/supp_try/carteVisitesALL.csv",header=T, sep=";")
# map <- base::subset(map,map$value==1 |map$value==0)
# map <- map[,1:2]
#
# nm <- SDMPlay:::null.model(predictors = predictors,xy=map,type=1,algorithm = "maxent",nb=200,unique.data = TRUE, same=T,nb.rep=4)
