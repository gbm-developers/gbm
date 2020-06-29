pkgname <- "biomod2"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('biomod2')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("BIOMOD.EnsembleModeling.out-objects")
### * BIOMOD.EnsembleModeling.out-objects

flush(stderr()); flush(stdout())

### Name: BIOMOD.EnsembleModeling.out-class
### Title: BIOMOD_EnsembleModeling() outputs objects class
### Aliases: BIOMOD.EnsembleModeling.out-class BIOMOD.EnsembleModeling.out
###   show,BIOMOD.EnsembleModeling.out-method
### Keywords: ensemble models

### ** Examples

showClass("BIOMOD.EnsembleModeling.out")




cleanEx()
nameEx("BIOMOD.Model.Options-objects")
### * BIOMOD.Model.Options-objects

flush(stderr()); flush(stdout())

### Name: BIOMOD.Model.Options-class
### Title: BIOMOD_ModelingOptions outputs objects class
### Aliases: BIOMOD.Model.Options-class show,BIOMOD.Model.Options-method
### Keywords: models options

### ** Examples

showClass("BIOMOD.Model.Options")




cleanEx()
nameEx("BIOMOD.formated.data-class")
### * BIOMOD.formated.data-class

flush(stderr()); flush(stdout())

### Name: BIOMOD.formated.data-class
### Title: BIOMOD_FormatingData() outputs objects class
### Aliases: BIOMOD.formated.data-class BIOMOD.formated.data
###   BIOMOD.formated.data.PA-class BIOMOD.formated.data.PA
###   BIOMOD.formated.data,data.frame,ANY-method
###   BIOMOD.formated.data,numeric,RasterStack-method
###   BIOMOD.formated.data,numeric,data.frame-method
###   BIOMOD.formated.data,numeric,matrix-method
###   show,BIOMOD.formated.data-method
###   plot,BIOMOD.formated.data,missing-method
###   show,BIOMOD.formated.data.PA-method
###   plot,BIOMOD.formated.data.PA,missing-method
### Keywords: models data formating

### ** Examples

showClass("BIOMOD.formated.data")



cleanEx()
nameEx("BIOMOD.models.out-class")
### * BIOMOD.models.out-class

flush(stderr()); flush(stdout())

### Name: BIOMOD.models.out-class
### Title: BIOMOD_modelling() outputs objects class
### Aliases: BIOMOD.models.out-class show,BIOMOD.models.out-method
### Keywords: models option

### ** Examples

showClass("BIOMOD.models.out")



cleanEx()
nameEx("BIOMOD.projection.out-class")
### * BIOMOD.projection.out-class

flush(stderr()); flush(stdout())

### Name: BIOMOD.projection.out-class
### Title: BIOMOD_Projection() outputs objects class
### Aliases: BIOMOD.projection.out-class BIOMOD.projection.out
###   show,BIOMOD.projection.out-method
###   plot,BIOMOD.projection.out,missing-method
### Keywords: models projection ensemble forecast

### ** Examples

showClass("BIOMOD.projection.out")



cleanEx()
nameEx("BIOMOD.stored.objects-class")
### * BIOMOD.stored.objects-class

flush(stderr()); flush(stdout())

### Name: BIOMOD.stored.objects-class
### Title: BIOMOD.stored.xxx objects class
### Aliases: BIOMOD.stored.data BIOMOD.stored.data-class
###   BIOMOD.stored.files BIOMOD.stored.files-class
###   BIOMOD.stored.data.frame BIOMOD.stored.data.frame-class
###   BIOMOD.stored.array BIOMOD.stored.array-class
###   BIOMOD.stored.formated.data BIOMOD.stored.formated.data-class
###   BIOMOD.stored.models.options BIOMOD.stored.models.options-class
###   BIOMOD.stored.models.out BIOMOD.stored.models.out-class
###   BIOMOD.stored.raster.stack BIOMOD.stored.raster.stack-class
### Keywords: models ensemble object storing

### ** Examples

showClass("BIOMOD.stored.files")



cleanEx()
nameEx("BIOMOD_EnsembleForecasting")
### * BIOMOD_EnsembleForecasting

flush(stderr()); flush(stdout())

### Name: BIOMOD_EnsembleForecasting
### Title: Ensemble projections of species over space and time
### Aliases: BIOMOD_EnsembleForecasting
### Keywords: models

### ** Examples

# 0. Load data & Selecting Data
# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species 
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd", 
                     package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio7.grd", 
                             package="biomod2"),  
                system.file( "external/bioclim/current/bio11.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio12.grd", 
                             package="biomod2"))

# 1. Formatting Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)
                                                                     
# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Running the models 
myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                       models = c('RF'), 
                                       models.options = myBiomodOption, 
                                       NbRunEval=2, 
                                       DataSplit=60, 
                                       Yweights=NULL, 
                                       VarImport=0, 
                                       models.eval.meth = c('TSS'),
                                       SaveObj = TRUE,
                                       rescal.all.models = FALSE,
                                       do.full.models = FALSE)
                                       
# 4. Creating the ensemble models 
myBiomodEM <- BIOMOD_EnsembleModeling( 
                 modeling.output = myBiomodModelOut,
                 chosen.models = grep('_RF', get_built_models(myBiomodModelOut), 
                                  value=TRUE),
                 em.by = 'algo',
                 eval.metric = c('TSS'),
                 eval.metric.quality.threshold = c(0.7),
                 prob.mean = TRUE,
                 prob.cv = FALSE,
                 prob.ci = FALSE,
                 prob.ci.alpha = 0.05,
                 prob.median = FALSE,
                 committee.averaging = FALSE,
                 prob.mean.weight = FALSE,
                 prob.mean.weight.decay = 'proportional' )
                                       
# 5. Individual models projections on current environmental conditions
myBiomodProjection <- BIOMOD_Projection(
                        modeling.output = myBiomodModelOut,
                        new.env = myExpl,
                        proj.name = 'current',
                        selected.models = grep('_RF', get_built_models(
                                                myBiomodModelOut), value=TRUE),
                        compress = FALSE,
                        build.clamping.mask = FALSE)
                          

# 4. Creating the ensemble projections
BIOMOD_EnsembleForecasting( projection.output = myBiomodProjection,
                            EM.output = myBiomodEM)



cleanEx()
nameEx("BIOMOD_EnsembleModeling")
### * BIOMOD_EnsembleModeling

flush(stderr()); flush(stdout())

### Name: BIOMOD_EnsembleModeling
### Title: Create and evaluate an ensemble set of models and predictions
### Aliases: BIOMOD_EnsembleModeling
### Keywords: models

### ** Examples

# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd",
                     package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio7.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio11.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio12.grd",
                             package="biomod2"))

# 1. Formatting Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation

myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,
                                       models = c('SRE','CTA','RF'),
                                       models.options = myBiomodOption,
                                       NbRunEval=1,
                                       DataSplit=80,
                                       Yweights=NULL,
                                       VarImport=3,
                                       models.eval.meth = c('TSS'),
                                       SaveObj = TRUE,
                                       rescal.all.models = FALSE,
                                       do.full.models = FALSE)

# 4. Doing Ensemble Modelling
myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
                           chosen.models = 'all',
                           em.by = 'all',
                           eval.metric = c('TSS'),
                           eval.metric.quality.threshold = c(0.7),
                           models.eval.meth = c('TSS','ROC'),
                           prob.mean = TRUE,
                           prob.cv = FALSE,
                           prob.ci = FALSE,
                           prob.ci.alpha = 0.05,
                           prob.median = FALSE,
                           committee.averaging = FALSE,
                           prob.mean.weight = TRUE,
                           prob.mean.weight.decay = 'proportional' )

# print summary
myBiomodEM

# get evaluation scores
get_evaluations(myBiomodEM)





cleanEx()
nameEx("BIOMOD_FormatingData")
### * BIOMOD_FormatingData

flush(stderr()); flush(stdout())

### Name: BIOMOD_FormatingData
### Title: Initialize the datasets for usage in 'biomod2'
### Aliases: BIOMOD_FormatingData
### Keywords: models datasets

### ** Examples


# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd",
                     package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio7.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio11.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio12.grd",
                             package="biomod2"))
# 1. Formatting Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

myBiomodData
plot(myBiomodData)




cleanEx()
nameEx("BIOMOD_LoadModels")
### * BIOMOD_LoadModels

flush(stderr()); flush(stdout())

### Name: BIOMOD_LoadModels
### Title: Load models built within BIOMOD_Modeling function
### Aliases: BIOMOD_LoadModels
### Keywords: models datasets

### ** Examples

# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd",
                     package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio7.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio11.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio12.grd",
                             package="biomod2"))

# 1. Formatting Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)


# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation

myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,
                                       models = c('RF'),
                                       models.options = myBiomodOption,
                                       NbRunEval=2,
                                       DataSplit=70,
                                       models.eval.meth = c('TSS'),
                                       SaveObj = TRUE,
                                       do.full.models = FALSE)

# 4. Loading some models built

myLoadedModels <- BIOMOD_LoadModels(myBiomodModelOut, models='RF')

myLoadedModels




cleanEx()
nameEx("BIOMOD_Modeling")
### * BIOMOD_Modeling

flush(stderr()); flush(stdout())

### Name: BIOMOD_Modeling
### Title: Run a range of species distribution models
### Aliases: BIOMOD_Modeling
### Keywords: models multivariate nonlinear nonparametric regression tree

### ** Examples

##' species occurrences
DataSpecies <- 
  read.csv(
    system.file(
      "external/species/mammals_table.csv",
      package="biomod2"
    )
  )
head(DataSpecies)

##' the name of studied species
myRespName <- 'GuloGulo'

##' the presence/absences data for our species
myResp <- as.numeric(DataSpecies[, myRespName])

##' the XY coordinates of species data
myRespXY <- DataSpecies[, c("X_WGS84", "Y_WGS84")]


##' Environmental variables extracted from BIOCLIM (bio_3, 
##' bio_4, bio_7, bio_11 & bio_12)
myExpl <- 
  raster::stack(
    system.file("external/bioclim/current/bio3.grd", package = "biomod2"),
    system.file("external/bioclim/current/bio4.grd", package = "biomod2"),
    system.file("external/bioclim/current/bio7.grd", package = "biomod2"),
    system.file("external/bioclim/current/bio11.grd", package = "biomod2"),
    system.file("external/bioclim/current/bio12.grd", package = "biomod2")
  )

##' 1. Formatting Data
myBiomodData <- 
  BIOMOD_FormatingData(
    resp.var = myResp,
    expl.var = myExpl,
    resp.xy = myRespXY,
    resp.name = myRespName
  )

##' 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

##' 3. Doing Modelisation
myBiomodModelOut <- 
  BIOMOD_Modeling(
    myBiomodData,
    models = c('SRE','RF'),
    models.options = myBiomodOption,
    NbRunEval = 2,
    DataSplit = 80,
    VarImport = 0,
    models.eval.meth = c('TSS','ROC'),
    do.full.models = FALSE,
    modeling.id = "test"
  )

##' print a summary of modeling stuff
myBiomodModelOut




cleanEx()
nameEx("BIOMOD_ModelingOptions")
### * BIOMOD_ModelingOptions

flush(stderr()); flush(stdout())

### Name: BIOMOD_ModelingOptions
### Title: Configure the modeling options for each selected model
### Aliases: BIOMOD_ModelingOptions
### Keywords: models options

### ** Examples

  ## default BIOMOD.model.option object
  myBiomodOptions <- BIOMOD_ModelingOptions()

  ## print the object
  myBiomodOptions

  ## you can copy a part of the print, change it and custom your options
  ## here we want to compute quadratic GLM and select best model with 'BIC' criterium
  myBiomodOptions <- BIOMOD_ModelingOptions(
    GLM = list( type = 'quadratic',
                interaction.level = 0,
                myFormula = NULL,
                test = 'BIC',
                family = 'binomial',
                control = glm.control(epsilon = 1e-08,
                                      maxit = 1000,
                                      trace = FALSE) ))

  ## check changes was done
  myBiomodOptions

  ##' you can prefer to establish your own GLM formula
  myBiomodOptions <- BIOMOD_ModelingOptions(
    GLM = list( myFormula = formula("Sp277 ~ bio3 +
                    log(bio10) + poly(bio16,2) + bio19 + bio3:bio19")))

  ## check changes was done
  myBiomodOptions

  ##' you also can directly print default parameters and then follow the same processus
  Print_Default_ModelingOptions()




cleanEx()
nameEx("BIOMOD_Projection")
### * BIOMOD_Projection

flush(stderr()); flush(stdout())

### Name: BIOMOD_Projection
### Title: Project the calibrated models within 'biomod2' into new space or
###   time
### Aliases: BIOMOD_Projection
### Keywords: models regression nonlinear multivariate nonparametric tree

### ** Examples

# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd",
                     package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio7.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio11.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio12.grd",
                             package="biomod2"))
# 1. Formatting Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation

myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,
                                       models = c('SRE','RF'),
                                       models.options = myBiomodOption,
                                       NbRunEval=1,
                                       DataSplit=70,
                                       models.eval.meth = c('TSS'),
                                       do.full.models = FALSE)


# 4.1 Projection on current environemental conditions

myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                          new.env = myExpl,
                                          proj.name = 'current',
                                          selected.models = 'all',
                                          binary.meth = 'TSS',
                                          compress = FALSE,
                                          build.clamping.mask = FALSE)

## Not run: 
##D # 4.2 Projection on future environemental conditions
##D myExplFuture = raster::stack(system.file("external/bioclim/future/bio3.grd",package="biomod2"),
##D                      system.file("external/bioclim/future/bio4.grd",package="biomod2"),
##D                      system.file("external/bioclim/future/bio7.grd",package="biomod2"),
##D                      system.file("external/bioclim/future/bio11.grd",package="biomod2"),
##D                      system.file("external/bioclim/future/bio12.grd",package="biomod2"))
##D 
##D myBiomodProjectionFuture <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
##D                                               new.env = myExplFuture,
##D                                               proj.name = 'future',
##D                                               selected.models = 'all',
##D                                               binary.meth = 'TSS',
##D                                               compress = FALSE,
##D                                               build.clamping.mask = TRUE)
##D 
##D # print summary and plot projections
##D myBiomodProjectionFuture
##D plot(myBiomodProjectionFuture)
## End(Not run)




cleanEx()
nameEx("BIOMOD_RangeSize")
### * BIOMOD_RangeSize

flush(stderr()); flush(stdout())

### Name: BIOMOD_RangeSize
### Title: Analysis of the range size changes
### Aliases: BIOMOD_RangeSize BIOMOD_RangeSize-methods
###   BIOMOD_RangeSize,data.frame,data.frame-method
###   BIOMOD_RangeSize,array,array-method
###   BIOMOD_RangeSize,RasterStack,RasterStack-method
###   BIOMOD_RangeSize,RasterLayer,RasterLayer-method
###   BIOMOD_RangeSize,RasterLayer,RasterStack-method
### Keywords: IO

### ** Examples

# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species 
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd", 
                             package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio7.grd", 
                             package="biomod2"),  
                system.file( "external/bioclim/current/bio11.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio12.grd", 
                             package="biomod2"))

# 1. Formatting Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)


# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation

myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                       models = c('CTA','RF'), 
                                       models.options = myBiomodOption, 
                                       models.eval.meth ='TSS',
                                       rescal.all.models=FALSE)


# 4.1 Projection on current environemental conditions

myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                          new.env = myExpl,
                                          proj.name = 'current',
                                          selected.models = 'all',
                                          binary.meth = 'TSS',
                                          compress = FALSE,
                                          build.clamping.mask = FALSE)

# 4.2 Projection on future environemental conditions

myExplFuture = raster::stack(system.file("external/bioclim/future/bio3.grd",package="biomod2"),
                     system.file("external/bioclim/future/bio4.grd",package="biomod2"),
                     system.file("external/bioclim/future/bio7.grd",package="biomod2"),
                     system.file("external/bioclim/future/bio11.grd",package="biomod2"),
                     system.file("external/bioclim/future/bio12.grd",package="biomod2"))

myBiomodProjectionFuture <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                              new.env = myExplFuture,
                                              proj.name = 'future',
                                              selected.models = 'all',
                                              binary.meth = 'TSS',
                                              compress = FALSE,
                                              build.clamping.mask = TRUE)

# 5. Detect where our species occurances state is forecasted to change

# load binary projections
# here is rasters objects ('.grd')
currentPred <- raster::stack("GuloGulo/proj_current/proj_current_GuloGulo_TSSbin.grd")
futurePred <- raster::stack("GuloGulo/proj_future/proj_future_GuloGulo_TSSbin.grd")


# call the Range size function
myBiomodRangeSize <- BIOMOD_RangeSize(
  CurrentPred=currentPred,
  FutureProj=futurePred)

# see the results
myBiomodRangeSize$Compt.By.Models
plot(myBiomodRangeSize$Diff.By.Pixel)




cleanEx()
nameEx("BIOMOD_cv")
### * BIOMOD_cv

flush(stderr()); flush(stdout())

### Name: BIOMOD_cv
### Title: Custom models cross-validation procedure
### Aliases: BIOMOD_cv

### ** Examples

## Not run: 
##D # species occurrences
##D DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##D                                     package="biomod2"))
##D head(DataSpecies)
##D 
##D the name of studied species
##D myRespName <- 'GuloGulo'
##D 
##D # the presence/absences data for our species 
##D myResp <- as.numeric(DataSpecies[,myRespName])
##D 
##D # the XY coordinates of species data
##D myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
##D 
##D 
##D # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##D myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
##D                              package="biomod2"),
##D                 system.file( "external/bioclim/current/bio4.grd", 
##D                              package="biomod2"), 
##D                 system.file( "external/bioclim/current/bio7.grd", 
##D                              package="biomod2"),  
##D                 system.file( "external/bioclim/current/bio11.grd", 
##D                              package="biomod2"), 
##D                 system.file( "external/bioclim/current/bio12.grd", 
##D                              package="biomod2"))
##D 
##D # 1. Formatting Data
##D myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##D                                      expl.var = myExpl,
##D                                      resp.xy = myRespXY,
##D                                      resp.name = myRespName)
##D 
##D # 2. Defining Models Options using default options.
##D myBiomodOption <- BIOMOD_ModelingOptions()
##D 
##D 
##D # 3. Creating DataSplitTable
##D 
##D DataSplitTable <- BIOMOD_cv(myBiomodData, k=5, rep=2, do.full.models=F)
##D DataSplitTable.y <- BIOMOD_cv(myBiomodData,stratified.cv=T, stratify="y", k=2)
##D colnames(DataSplitTable.y)[1:2] <- c("RUN11","RUN12")
##D DataSplitTable <- cbind(DataSplitTable,DataSplitTable.y)
##D head(DataSplitTable)
##D 
##D # 4. Doing Modelisation
##D 
##D myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
##D                                      models = c('RF'), 
##D                                      models.options = myBiomodOption, 
##D                                      DataSplitTable = DataSplitTable,
##D                                      VarImport=0, 
##D                                      models.eval.meth = c('ROC'),
##D                                      do.full.models=FALSE,
##D                                      modeling.id="test")
##D 
##D ## get cv evaluations
##D eval <- get_evaluations(myBiomodModelOut,as.data.frame=T)
##D 
##D eval$strat <- NA
##D eval$strat[grepl("13",eval$Model.name)] <- "Full"
##D eval$strat[!(grepl("11",eval$Model.name)|
##D              grepl("12",eval$Model.name)|
##D              grepl("13",eval$Model.name))] <- "Random"
##D eval$strat[grepl("11",eval$Model.name)|grepl("12",eval$Model.name)] <- "Strat"
##D 
##D boxplot(eval$Testing.data~ eval$strat, ylab="ROC AUC")
## End(Not run)



cleanEx()
nameEx("BIOMOD_presenceonly")
### * BIOMOD_presenceonly

flush(stderr()); flush(stdout())

### Name: BIOMOD_presenceonly
### Title: evaluate models with presences only metrics
### Aliases: BIOMOD_presenceonly

### ** Examples

## Not run: 
##D requireNamesapce(PresenceAbsence, 'PresenceAbsence', quietly = TRUE)
##D 
##D # species occurrences
##D DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##D                                     package="biomod2"), row.names = 1)
##D head(DataSpecies)
##D 
##D # the name of studied species
##D myRespName <- 'GuloGulo'
##D 
##D # the presence/absences data for our species 
##D myResp <- as.numeric(DataSpecies[,myRespName])
##D 
##D # the XY coordinates of species data
##D myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
##D 
##D 
##D # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##D myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
##D                              package="biomod2"),
##D                 system.file( "external/bioclim/current/bio4.grd", 
##D                              package="biomod2"), 
##D                 system.file( "external/bioclim/current/bio7.grd", 
##D                              package="biomod2"),  
##D                 system.file( "external/bioclim/current/bio11.grd", 
##D                              package="biomod2"), 
##D                 system.file( "external/bioclim/current/bio12.grd", 
##D                              package="biomod2"))
##D 
##D # 1. Formatting Data
##D myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##D                                      expl.var = myExpl,
##D                                      resp.xy = myRespXY,
##D                                      resp.name = myRespName)
##D 
##D # 2. Defining Models Options using default options.
##D myBiomodOption <- BIOMOD_ModelingOptions()
##D 
##D # 3. Doing Modelisation
##D 
##D myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
##D                                      models = c('SRE','CTA','RF'), 
##D                                      models.options = myBiomodOption, 
##D                                      NbRunEval=1, 
##D                                      DataSplit=80, 
##D                                      Yweights=NULL, 
##D                                      VarImport=3, 
##D                                      models.eval.meth = c('TSS','ROC'),
##D                                      SaveObj = TRUE,
##D                                      rescal.all.models = FALSE,
##D                                      do.full.models = FALSE)
##D 
##D # 4. Doing Ensemble Modelling
##D myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
##D                                        chosen.models = 'all',
##D                                        em.by = 'PA_dataset+repet',
##D                                        eval.metric = c('TSS'),
##D                                        eval.metric.quality.threshold = c(0.7),
##D                                        models.eval.meth = c('TSS','ROC'),
##D                                        prob.mean = TRUE,
##D                                        prob.cv = FALSE,
##D                                        prob.ci = FALSE,
##D                                        prob.ci.alpha = 0.05,
##D                                        prob.median = FALSE,
##D                                        committee.averaging = FALSE,
##D                                        prob.mean.weight = TRUE,
##D                                        prob.mean.weight.decay = 'proportional' )   
##D 
##D # evaluate Biomod models with the Boyce index and MPA
##D pres.only.eval <- BIOMOD_presenceonly(myBiomodModelOut, myBiomodEM)
##D pres.only.eval$eval
##D 
##D # evaluate Biomod models with the Boyce index and MPA using Background data
##D bg.Values <- getValues(myExpl)
##D 
##D pres.only.eval <- BIOMOD_presenceonly(myBiomodModelOut, myBiomodEM, bg.env = bg.Values)
##D pres.only.eval$eval
## End(Not run)



cleanEx()
nameEx("BIOMOD_tuning")
### * BIOMOD_tuning

flush(stderr()); flush(stdout())

### Name: BIOMOD_tuning
### Title: Tune models parameters
### Aliases: BIOMOD_tuning

### ** Examples

## Not run: 
##D # species occurrences
##D DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##D                                     package="biomod2"))
##D head(DataSpecies)
##D 
##D # the name of studied species
##D myRespName <- 'GuloGulo'
##D 
##D # the presence/absences data for our species 
##D myResp <- as.numeric(DataSpecies[,myRespName])
##D 
##D # the XY coordinates of species data
##D myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
##D 
##D # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##D myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
##D                              package="biomod2"),
##D                 system.file( "external/bioclim/current/bio4.grd", 
##D                              package="biomod2"), 
##D                 system.file( "external/bioclim/current/bio7.grd", 
##D                              package="biomod2"),  
##D                 system.file( "external/bioclim/current/bio11.grd", 
##D                              package="biomod2"), 
##D                 system.file( "external/bioclim/current/bio12.grd", 
##D                              package="biomod2"))
##D # 1. Formatting Data
##D myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##D                                      expl.var = myExpl,
##D                                      resp.xy = myRespXY,
##D                                      resp.name = myRespName)
##D 
##D # 2. Defining Models Options using default options.
##D ### Duration for turing all models sequential with default settings 
##D ### on 3.4 GHz processor: approx. 45 min tuning all models in parallel
##D ### (on 8 cores) using foreach loops runs much faster: approx. 14 min
##D 
##D #library(doParallel);cl<-makeCluster(8);doParallel::registerDoParallel(cl) 
##D 
##D 
##D time.seq<-system.time(Biomod.tuning <- BIOMOD_tuning(myBiomodData,
##D                                                              env.ME = myExpl,
##D                                                              n.bg.ME = ncell(myExpl)))
##D #stopCluster(cl)
##D 
##D myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
##D                                      models = c('RF','CTA'), 
##D                                      models.options = Biomod.tuning$models.options, 
##D                                      NbRunEval=1, 
##D                                      DataSplit=100, 
##D                                      VarImport=0, 
##D                                      models.eval.meth = c('ROC'),
##D                                      do.full.models=FALSE,
##D                                      modeling.id="test")
##D 
##D 
##D #  eval.plot(Biomod.tuning$tune.MAXENT.Phillips at results)
##D par(mfrow=c(1,3))
##D plot(Biomod.tuning$tune.CTA.rpart)
##D plot(Biomod.tuning$tune.CTA.rpart2)
##D plot(Biomod.tuning$tune.RF)
## End(Not run)



cleanEx()
nameEx("BinaryTransformation-methods")
### * BinaryTransformation-methods

flush(stderr()); flush(stdout())

### Name: BinaryTransformation
### Title: Convert species' probability of occurrence into binary
###   presence-absence data using a predefined threshold
### Aliases: BinaryTransformation BinaryTransformation,data.frame-method
###   BinaryTransformation, data.frame-method
###   BinaryTransformation,matrix-method matrix-method
###   BinaryTransformation,numeric-method numeric-method
###   BinaryTransformation,array-method array-method
###   BinaryTransformation,RasterLayer-method RasterLayer-method
###   BinaryTransformation,RasterStack-method RasterStack-method
###   BinaryTransformation,RasterBrick-method RasterBrick-method
### Keywords: models

### ** Examples

  xx <- rnorm(50,10)
  yy <- BinaryTransformation(xx, 10)

  cbind(xx,yy)




cleanEx()
nameEx("CustomIndexMaker")
### * CustomIndexMaker

flush(stderr()); flush(stdout())

### Name: CustomIndexMaker
### Title: Replace default package Index help file by a custom one.
### Aliases: CustomIndexMaker
### Keywords: models datasets

### ** Examples


## Automaticly done at buildinfg package state
# CustomIndexMaker()




cleanEx()
nameEx("DF_to_ARRAY")
### * DF_to_ARRAY

flush(stderr()); flush(stdout())

### Name: DF_to_ARRAY
### Title: Convert a biomod2 data.frame (or list) into array
### Aliases: DF_to_ARRAY LIST_to_ARRAY
### Keywords: models formula options

### ** Examples


# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species 
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd", 
                     package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio7.grd", 
                             package="biomod2"),  
                system.file( "external/bioclim/current/bio11.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio12.grd", 
                             package="biomod2"))
                             
# Keep only points where we have info                       
myExpl <- raster::extract(myExpl, myRespXY)

# 1. Formatting Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)
                                                                     
# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation

myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                       models = c('SRE','RF'), 
                                       models.options = myBiomodOption, 
                                       NbRunEval=1, 
                                       DataSplit=70, 
                                       Yweights=NULL, 
                                       VarImport=0, 
                                       models.eval.meth = c('ROC'),
                                       rescal.all.models = FALSE,
                                       do.full.models = FALSE)
                                       
                                       
# 4 Projection on current environemental conditions

myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                          new.env = data.frame(myExpl),
                                          proj.name = 'current',
                                          selected.models = 'all')
                                          


# 5. Get projection under data.frame format
myProjDF <- get_predictions(myBiomodProjection, as.data.frame=TRUE)
class(myProjDF)
dim(myProjDF)
dimnames(myProjDF)

# 6. Transform data.frame into array
myProjArray <- DF_to_ARRAY(myProjDF)
class(myProjArray)
dim(myProjArray)
dimnames(myProjArray)




cleanEx()
nameEx("FilteringTransformation")
### * FilteringTransformation

flush(stderr()); flush(stdout())

### Name: FilteringTransformation
### Title: Convert species' probability of occurrence into binary
###   presence-absence data using a predefined threshold
### Aliases: FilteringTransformation FilteringTransformation-methods
###   FilteringTransformation,data.frame-method
###   FilteringTransformation,matrix-method
###   FilteringTransformation,numeric-method
###   FilteringTransformation,array-method
###   FilteringTransformation,RasterBrick-method
###   FilteringTransformation,RasterLayer-method
###   FilteringTransformation,RasterStack-method

### ** Examples

xx <- rnorm(50,10)
yy <- FilteringTransformation(xx, 10)

cbind(xx,yy)



cleanEx()
nameEx("Find.Optim.Stat")
### * Find.Optim.Stat

flush(stderr()); flush(stdout())

### Name: Find.Optim.Stat
### Title: Calculate the best score according to a given evaluation method
### Aliases: Find.Optim.Stat
### Keywords: evaluation models options

### ** Examples

  a <- sample(c(0,1),100, replace=TRUE)

  ##' random drawing
  b <- runif(100,min=0,max=1000)
  Find.Optim.Stat(Stat='TSS',
                  Fit=b,
                  Obs=a)

  ##' biased drawing
  BiasedDrawing <- function(x, m1=300, sd1=200, m2=700, sd2=200){
    return(ifelse(x<0.5, rnorm(1,m1,sd1), rnorm(1,m2,sd2)))
  }

  c <- sapply(a,BiasedDrawing)

  Find.Optim.Stat(Stat='TSS',
                  Fit=c,
                  Obs=a,
                  Nb.thresh.test = 100)





cleanEx()
nameEx("Print_Default_ModelingOptions")
### * Print_Default_ModelingOptions

flush(stderr()); flush(stdout())

### Name: Print_Default_ModelingOptions
### Title: Get default values of BIOMOD inner models' options
### Aliases: Print_Default_ModelingOptions
### Keywords: models options

### ** Examples

# print default models options
Print_Default_ModelingOptions()



cleanEx()
nameEx("ProbDensFunc")
### * ProbDensFunc

flush(stderr()); flush(stdout())

### Name: ProbDensFunc
### Title: Probability Density Function
### Aliases: ProbDensFunc
### Keywords: distribution optimize

### ** Examples

## Not run: 
##D DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##D                                     package="biomod2"), row.names = 1)
##D head(DataSpecies)
##D 
##D ##' the name of studied species
##D myRespName <- 'GuloGulo'
##D 
##D ##' the presence/absences data for our species
##D myResp <- as.numeric(DataSpecies[,myRespName])
##D 
##D ##' remove all 0 from response vector to work with
##D ##' presence only data (Pseudo Absences selections)
##D rm_id <- which(myResp==0)
##D myResp <- myResp[-rm_id]
##D 
##D 
##D ##' the XY coordinates of species data
##D myRespXY <- DataSpecies[-rm_id,c("X_WGS84","Y_WGS84")]
##D 
##D 
##D ##' Environmental variables extracted from BIOCLIM
##D myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd",
##D                              package="biomod2"),
##D                 system.file( "external/bioclim/current/bio4.grd",
##D                              package="biomod2"),
##D                 system.file( "external/bioclim/current/bio7.grd",
##D                              package="biomod2"),
##D                 system.file( "external/bioclim/current/bio11.grd",
##D                              package="biomod2"),
##D                 system.file( "external/bioclim/current/bio12.grd",
##D                              package="biomod2"))
##D 
##D ##' 1. Formatting Data
##D myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##D                                      expl.var = myExpl,
##D                                      resp.xy = myRespXY,
##D                                      resp.name = myRespName,
##D                                      PA.nb.rep=3)
##D 
##D ##' 2. Defining Models Options using default options.
##D myBiomodOption <- BIOMOD_ModelingOptions()
##D 
##D ##' 3. Doing Modelisation
##D myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,
##D                                      models = c('CTA','RF','GLM','GAM','ANN','MARS'),
##D                                      models.options = myBiomodOption,
##D                                      NbRunEval=5,
##D                                      DataSplit=70,
##D                                      Prevalence=0.5,
##D                                      models.eval.meth = c('TSS'),
##D                                      do.full.models = FALSE,
##D                                      rescal.all.models=T,
##D                                      modeling.id='test')
##D 
##D ##' 4. Build ensemble-models that will be taken as reference
##D myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
##D                                        chosen.models = 'all',
##D                                        em.by = 'all',
##D                                        eval.metric = c('TSS'),
##D                                        eval.metric.quality.threshold = c(0.7),
##D                                        prob.mean = TRUE,
##D                                        prob.median = TRUE)
##D 
##D ##' 5. Projection on future environmental conditions
##D 
##D ###' load future environmental conditions from biomod2 package
##D myExpl_fut <- raster::stack( system.file( "external/bioclim/future/bio3.grd",
##D                                   package="biomod2"),
##D                      system.file( "external/bioclim/future/bio4.grd",
##D                                   package="biomod2"),
##D                      system.file( "external/bioclim/future/bio7.grd",
##D                                   package="biomod2"),
##D                      system.file( "external/bioclim/future/bio11.grd",
##D                                   package="biomod2"),
##D                      system.file( "external/bioclim/future/bio12.grd",
##D                                   package="biomod2"))
##D 
##D myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
##D                                         new.env = myExpl_fut,
##D                                         proj.name = 'future',
##D                                         selected.models = 'all',
##D                                         binary.meth = 'TSS',
##D                                         compress = FALSE,
##D                                         build.clamping.mask = TRUE)
##D 
##D BIOMOD_EnsembleForecasting(projection.output=myBiomodProjection,
##D                            EM.output=myBiomodEM,
##D                            binary.meth='TSS')
##D 
##D ##' 6. load binary projections
##D consensusBin <- raster::stack('GuloGulo/proj_future/proj_future_GuloGulo_ensemble_TSSbin.grd')
##D projectionsBin <- raster::stack('GuloGulo/proj_future/proj_future_GuloGulo_TSSbin.grd')
##D 
##D ##' 7. build a ref state based on ensemble-models
##D ref <- sampleRandom(subset(consensusBin, 1, drop=T), size=5000, sp=T, na.rm=T)
##D 
##D ##' 8. autoatic creation of groups matrix
##D find_groups <- function(diff_by_pix){
##D   data.set <- sapply(names(diff_by_pix),biomod2:::.extractModelNamesInfo,info='data.set')
##D   run.eval <- sapply(names(diff_by_pix),biomod2:::.extractModelNamesInfo,info='run.eval')
##D   models <- sapply(names(diff_by_pix),biomod2:::.extractModelNamesInfo,info='models')
##D   return(rbind(data.set,run.eval,models))
##D }
##D 
##D groups <- find_groups(projectionsBin)
##D 
##D ##' 9. plot ProbDensFunct graphs
##D ProbDensFunc(initial = ref,
##D              projections = projectionsBin,
##D              plothist=TRUE,
##D              cvsn=TRUE,
##D              groups=groups,
##D              resolution=2,
##D              filename=NULL,
##D              lim=c(0.5,0.8,0.95))
##D 
##D ###' 3 plots should be produced.. Should be convenient to save it within a device
##D ###' supporting multiple plots.
##D 
## End(Not run)




cleanEx()
nameEx("SampleMat2")
### * SampleMat2

flush(stderr()); flush(stdout())

### Name: SampleMat2
### Title: Sample binary vector
### Aliases: SampleMat2
### Keywords: formula models options

### ** Examples

a <- sample(c(0,1),100, replace=TRUE)
SampleMat2(ref=a, ratio=0.7)




cleanEx()
nameEx("biomod2_model-class")
### * biomod2_model-class

flush(stderr()); flush(stdout())

### Name: biomod2_model-class
### Title: biomod2 models objects class and functions
### Aliases: biomod2_model-class biomod2_model biomod2_ensemble_model-class
###   biomod2_ensemble_model show,biomod2_model-method get_formal_model
###   get_formal_model,biomod2_model-method get_scaling_model
###   get_scaling_model,biomod2_model-method check_data_range get_var_range
###   get_var_type ANN_biomod2_model-class ANN_biomod2_model
###   predict,ANN_biomod2_model-method CTA_biomod2_model-class
###   CTA_biomod2_model predict,CTA_biomod2_model-method
###   FDA_biomod2_model-class FDA_biomod2_model
###   predict,FDA_biomod2_model-method GAM_biomod2_model-class
###   GAM_biomod2_model predict,GAM_biomod2_model-method
###   GLM_biomod2_model-class GLM_biomod2_model
###   predict,GLM_biomod2_model-method GBM_biomod2_model-class
###   GBM_biomod2_model predict,GBM_biomod2_model-method
###   MARS_biomod2_model-class MARS_biomod2_model
###   predict,MARS_biomod2_model-method MAXENT.Phillips_biomod2_model-class
###   MAXENT.Phillips_biomod2_model
###   predict,MAXENT.Phillips_biomod2_model-method
###   MAXENT.Phillips.2_biomod2_model-class MAXENT.Phillips.2_biomod2_model
###   predict,MAXENT.Phillips.2_biomod2_model-method RF_biomod2_model-class
###   RF_biomod2_model predict,RF_biomod2_model-method
###   SRE_biomod2_model-class SRE_biomod2_model
###   predict,SRE_biomod2_model-method EMca_biomod2_model-class
###   EMca_biomod2_model predict,EMca_biomod2_model-method
###   EMci_biomod2_model-class EMci_biomod2_model
###   predict,EMci_biomod2_model-method EMcv_biomod2_model-class
###   EMcv_biomod2_model predict,EMcv_biomod2_model-method
###   EMmean_biomod2_model-class EMmean_biomod2_model
###   predict,EMmean_biomod2_model-method EMmedian_biomod2_model-class
###   EMmedian_biomod2_model predict,EMmedian_biomod2_model-method
###   EMwmean_biomod2_model-class EMwmean_biomod2_model
###   predict,EMwmean_biomod2_model-method
### Keywords: models predict

### ** Examples

showClass("ANN_biomod2_model")



cleanEx()
nameEx("calculate.stat")
### * calculate.stat

flush(stderr()); flush(stdout())

### Name: calculate.stat
### Title: Calculate evaluation metrics based on a misclassification table
### Aliases: calculate.stat
### Keywords: models formula options

### ** Examples

  a <- sample(c(0,1),100, replace=TRUE)
  b <- sample(c(0,1),100, replace=TRUE)
  
  miscTab_aa <- table(a,a)
  miscTab_ab <- table(a,b)
  
  # perfect score
  calculate.stat( miscTab_aa, stat='TSS')
  # random score
  calculate.stat( miscTab_ab, stat='TSS')
  



cleanEx()
nameEx("evaluate")
### * evaluate

flush(stderr()); flush(stdout())

### Name: evaluate
### Title: biomod2 modelling outputs evaluation
### Aliases: evaluate
### Keywords: evaluation models score

### ** Examples



# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species 
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd", 
                     package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio7.grd", 
                             package="biomod2"),  
                system.file( "external/bioclim/current/bio11.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio12.grd", 
                             package="biomod2"))

# 1. Formatting Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation

myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                     models = c('SRE','CTA','RF'), 
                                     models.options = myBiomodOption, 
                                     NbRunEval=1, 
                                     DataSplit=80, 
                                     Yweights=NULL, 
                                     VarImport=3, 
                                     models.eval.meth = c('TSS'),
                                     SaveObj = TRUE,
                                     rescal.all.models = FALSE,
                                     do.full.models = FALSE,
                                     modeling.id='test')
                                     
# 4. Evaluate model over another dataset (here the full one)

## creation of suitable dataset
data <- cbind(GuloGulo=get_formal_data(myBiomodModelOut,'resp.var'), 
              get_formal_data(myBiomodModelOut,'expl.var'))

## evaluation
evaluate(myBiomodModelOut, data=data, stat=c('ROC','TSS'))




cleanEx()
nameEx("full_shuffling")
### * full_shuffling

flush(stderr()); flush(stdout())

### Name: full_suffling
### Title: data set shuffling tool
### Aliases: full_suffling
### Keywords: shuffle random importance

### ** Examples

xx <- matrix(rep(1:10,3),10,3)
full_suffling(xx,c(1,2))



cleanEx()
nameEx("getStatOptimValue")
### * getStatOptimValue

flush(stderr()); flush(stdout())

### Name: getStatOptimValue
### Title: get the optimal score of evaluation statistical metrics
### Aliases: getStatOptimValue
### Keywords: models formula options

### ** Examples

  getStatOptimValue('TSS')
  getStatOptimValue('KAPPA')
  getStatOptimValue('POFD')



cleanEx()
nameEx("level.plot")
### * level.plot

flush(stderr()); flush(stdout())

### Name: level.plot
### Title: Plot 2-dimensional data for visualizing distribution of species
###   or environment
### Aliases: level.plot
### Keywords: plot

### ** Examples

## Not run: 
##D # species occurrences
##D DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##D                                     package="biomod2"), row.names = 1)
##D 
##D # the name of studied species
##D myRespName <- 'GuloGulo'
##D 
##D # the presence/absences data for our species 
##D myResp <- as.numeric(DataSpecies[,myRespName])
##D 
##D # the XY coordinates of species data
##D myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
##D 
##D 
##D level.plot(data.in=myResp, XY=myRespXY)
## End(Not run)



cleanEx()
nameEx("makeFormula")
### * makeFormula

flush(stderr()); flush(stdout())

### Name: makeFormula
### Title: Standardized formula maker
### Aliases: makeFormula
### Keywords: formula models options

### ** Examples

##' create simulated data
myResp <- sample(c(0, 1), 20, replace = TRUE)
myExpl <- 
  matrix(
    runif(60), 
    ncol = 3, 
    dimnames=list(NULL, c('var1', 'var2', 'var3'))
  )

##' create a formula
myFormula <- 
  makeFormula( 
    respName = 'myResp',
    explVar = head(myExpl),
    type = 'quadratic',
    interaction.level = 0
  )
  
##' show formula created
myFormula




cleanEx()
nameEx("models_scores_graph")
### * models_scores_graph

flush(stderr()); flush(stdout())

### Name: models_scores_graph
### Title: Produce models evaluation bi-dimensional graph
### Aliases: models_scores_graph
### Keywords: evaluation scores graph

### ** Examples


## this example is based on BIOMOD_Modeling function example
example(BIOMOD_Modeling)

## we will need ggplot2 package to produce our custom version of the graphs
require(ggplot2)

## plot evaluation models score graph

### by models
gg1 <- models_scores_graph( myBiomodModelOut,
                            by = 'models',
                            metrics = c('ROC','TSS') )
## we see a influence of model selected on models capabilities
## e.g. RF are much better than SRE

### by cross validation run
gg2 <- models_scores_graph( myBiomodModelOut,
                            by = 'cv_run',
                            metrics = c('ROC','TSS') )
## there is no difference in models quality if we focus on 
## cross validation sampling

### some graphical customisations
gg1_custom <- 
  gg1 + 
  ggtitle("Diff between RF and SRE evaluation scores") + ## add title
  scale_colour_manual(values=c("green", "blue")) ## change colors

gg1_custom




cleanEx()
nameEx("multiple.plot")
### * multiple.plot

flush(stderr()); flush(stdout())

### Name: multiple.plot
### Title: Plot and compare prediction maps within BIOMOD
### Aliases: multiple.plot
### Keywords: plot

### ** Examples

## Not run: 
##D # species occurrences
##D DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##D                                     package="biomod2"), row.names = 1)
##D 
##D # the name of studied species
##D myRespName <- c("ConnochaetesGnou", "GuloGulo", "PantheraOnca", 
##D                 "PteropusGiganteus", "TenrecEcaudatus", "VulpesVulpes")
##D 
##D # the presence/absences data for our species 
##D myResp <- DataSpecies[,myRespName]
##D 
##D # the XY coordinates of species data
##D myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
##D 
##D multiple.plot(Data = myResp,
##D               coor = myRespCoord )
## End(Not run)



cleanEx()
nameEx("randomise_data")
### * randomise_data

flush(stderr()); flush(stdout())

### Name: randomise_data
### Title: data set shuffling tool
### Aliases: randomise_data
### Keywords: suffle random importance

### ** Examples

xx <- data.frame(a=1:10,b=11:20,c=21:30)
randomise_data(data=xx, variable='b', method='full_rand')



cleanEx()
nameEx("response.plot2")
### * response.plot2

flush(stderr()); flush(stdout())

### Name: response.plot2
### Title: Function for for plotting predicted responses from species
###   distribution models in 2 or 3 dimensions
### Aliases: response.plot2
### Keywords: models multivariate nonlinear nonparametric plot regression
###   tree

### ** Examples

## Not run: 
##D ##' species occurrences
##D DataSpecies <- 
##D   read.csv(
##D     system.file("external/species/mammals_table.csv", package="biomod2"), 
##D     row.names = 1
##D   )
##D head(DataSpecies)
##D ##' the name of studied species
##D myRespName <- 'VulpesVulpes'
##D     
##D ##' the presence/absences data for our species 
##D myResp <- as.numeric(DataSpecies[, myRespName])
##D     
##D ##' the XY coordinates of species data
##D myRespXY <- DataSpecies[, c("X_WGS84", "Y_WGS84")]
##D 
##D myExpl <- 
##D   raster::stack(
##D     system.file("external/bioclim/current/bio3.grd", package = "biomod2"),
##D     system.file("external/bioclim/current/bio4.grd", package = "biomod2"),
##D     system.file("external/bioclim/current/bio7.grd", package = "biomod2"),
##D     system.file("external/bioclim/current/bio11.grd", package = "biomod2"),
##D     system.file("external/bioclim/current/bio12.grd", package = "biomod2")
##D   )
##D 
##D ##' 1. Formatting Data
##D myBiomodData <- 
##D   BIOMOD_FormatingData(
##D     resp.var = myResp,
##D     expl.var = myExpl,
##D     resp.xy = myRespXY,
##D     resp.name = myRespName
##D   )
##D 
##D ##' 2. Defining Models Options using default options.
##D myBiomodOption <- BIOMOD_ModelingOptions()
##D 
##D ##' 3. Doing Modelisation
##D myBiomodModelOut <- 
##D   BIOMOD_Modeling(
##D     myBiomodData,
##D     models = c('GLM','RF'),
##D     models.options = myBiomodOption,
##D     NbRunEval = 2,
##D     DataSplit = 80,
##D     VarImport = 0,
##D     models.eval.meth = c('TSS','ROC'),
##D     do.full.models = FALSE,
##D     modeling.id = "test"
##D   )
##D ##' 4. Plot response curves
##D ##' 4.1 Load the models for which we want to extract the predicted
##D ##' response curves
##D myGLMs <- BIOMOD_LoadModels(myBiomodModelOut, models = 'GLM')
##D     
##D ##' 4.2 plot 2D response plots
##D myRespPlot2D <- 
##D   response.plot2(
##D     models = myGLMs,
##D     Data = get_formal_data(myBiomodModelOut, 'expl.var'),
##D     show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
##D     do.bivariate = FALSE,
##D     fixed.var.metric = 'median',
##D     col = c("blue", "red"),
##D     legend = TRUE,
##D     data_species = get_formal_data(myBiomodModelOut, 'resp.var')
##D   )
##D     
##D ##' 4.2 plot 3D response plots
##D ###' here only for a lone model (i.e "VulpesVulpes_PA1_AllData_GLM")
##D myRespPlot3D <- 
##D   response.plot2(
##D   models = myGLMs[1],
##D   Data = get_formal_data(myBiomodModelOut, 'expl.var'), 
##D   show.variables = get_formal_data(myBiomodModelOut, 'expl.var.names'),
##D   do.bivariate = TRUE,
##D   fixed.var.metric = 'median',
##D   data_species = get_formal_data(myBiomodModelOut, 'resp.var'),
##D   display_title = FALSE
##D )
##D     
##D ##' all the values used to produce this plot are stored into the
##D ##' returned object you can redo plots by yourself and customised 
##D ##' them
##D dim(myRespPlot2D)
##D dimnames(myRespPlot2D)
##D     
##D dim(myRespPlot3D)
##D dimnames(myRespPlot3D)
## End(Not run)




cleanEx()
nameEx("sample.factor.levels")
### * sample.factor.levels

flush(stderr()); flush(stdout())

### Name: sample.factor.levels
### Title: Tool to ensure the sampling of all levels of a factorial
###   variable
### Aliases: sample.factor.levels

### ** Examples

## example with raster* object ---------- 
library(raster)
## create a factorial raster
r1 <- raster()
r1[] <- 1; r1[1] <- 2; r1[2:3] <- 3
r1 <- as.factor(r1)
## create a continuous raster
r2 <- raster()
r2[] <- rnorm(ncell(r2))
## pull the raster into a RasterStack
stk <- stack(r1, r2)
is.factor(stk)

## define a mask for already sampled points
mask.out <- r1
mask.out[] <- NA; mask.out[2:3] <- 1

## define a list of mask where we want to sample in priority
mask.in.1 <- mask.in.2 <- r1
mask.in.1[1:10] <- NA ## only level 1 should be sampled in this mask
mask.in.2[1] <- NA ## only levels 1 and 3 should be sampled in this mask
mask.in <- list(mask.in.1 = mask.in.1, 
                mask.in.2 = mask.in.2)

## test different version of the function
sample.factor.levels(stk, mask.out = mask.out)
sample.factor.levels(stk, mask.in = mask.in)
sample.factor.levels(stk, mask.out = mask.out, mask.in = mask.in)




cleanEx()
nameEx("sre")
### * sre

flush(stderr()); flush(stdout())

### Name: sre
### Title: Surface Range Envelope
### Aliases: sre
### Keywords: models multivariate

### ** Examples

require(raster)
##' species occurrences
DataSpecies <- 
  read.csv(
    system.file("external/species/mammals_table.csv", package = "biomod2"), 
    row.names = 1
  )
head(DataSpecies)

##' the name of studied species
myRespName <- 'GuloGulo'

##' the presence/absences data for our species 
myResp <- as.numeric(DataSpecies[,myRespName])

##' the XY coordinates of species data
myRespXY <- DataSpecies[which(myResp==1),c("X_WGS84","Y_WGS84")]

##' Environmental variables extracted from BIOCLIM (bio_3, 
##' bio_4, bio_7, bio_11 & bio_12)
myExpl <- 
  raster::stack(
    system.file("external/bioclim/current/bio3.grd", package = "biomod2"),
    system.file("external/bioclim/current/bio4.grd", package = "biomod2"),
    system.file("external/bioclim/current/bio7.grd", package = "biomod2"),
    system.file("external/bioclim/current/bio11.grd", package = "biomod2"),
    system.file("external/bioclim/current/bio12.grd", package = "biomod2")
  )
myResp <- 
  raster::reclassify(
    subset(myExpl, 1, drop = TRUE), c(-Inf, Inf, 0)
  )
myResp[cellFromXY(myResp,myRespXY)] <- 1

##' Compute some SRE for several quantile values
sre.100 <- 
  sre(
    Response = myResp, 
    Explanatory = myExpl, 
    NewData=myExpl, 
    Quant = 0
  )
  
sre.095 <- 
  sre(
    Response = myResp, 
    Explanatory = myExpl, 
    NewData=myExpl, 
    Quant = 0.025
  )

sre.090 <- 
  sre(
    Response = myResp, 
    Explanatory = myExpl, 
    NewData=myExpl, 
    Quant = 0.05
  )
  
##' visualise results
par(mfrow=c(2,2),mar=c(6, 5, 5, 3))
plot(myResp, main = paste(myRespName, "original distrib."))
plot(sre.100, main="full data calibration")
plot(sre.095, main="95 %")
plot(sre.090, main="90 %")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("variables_importance")
### * variables_importance

flush(stderr()); flush(stdout())

### Name: variables_importance
### Title: Variables importance calculation
### Aliases: variables_importance
### Keywords: importance random suffle

### ** Examples

xx <- 
  data.frame( 
    a = sample(c(0, 1), 100, replace = TRUE),
    b = rnorm(100),
    c = 1:100
  )
  
mod <- glm(a ~ b + c, data = xx)

variables_importance(
  model = mod, 
  data = xx[, c('b', 'c')], 
  method = "full_rand", 
  nb_rand = 3
)




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
