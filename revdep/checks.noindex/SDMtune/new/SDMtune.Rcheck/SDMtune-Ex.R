pkgname <- "SDMtune"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SDMtune')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("SDMmodel2MaxEnt")
### * SDMmodel2MaxEnt

flush(stderr()); flush(stdout())

### Name: SDMmodel2MaxEnt
### Title: SDMmodel2MaxEnt
### Aliases: SDMmodel2MaxEnt

### ** Examples





cleanEx()
nameEx("addSamplesToBg")
### * addSamplesToBg

flush(stderr()); flush(stdout())

### Name: addSamplesToBg
### Title: Add Samples to Background
### Aliases: addSamplesToBg

### ** Examples

# Acquire environmental variables
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)

# Prepare presence and background locations
p_coords <- virtualSp$presence
bg_coords <- virtualSp$background

# Create SWD object
data <- prepareSWD(species = "Virtual species", p = p_coords, a = bg_coords,
                   env = predictors, categorical = "biome")

# Add presence locations with values not included in the background to the
# background locations
new_data <- addSamplesToBg(data)
new_data

# Add all the presence locations to the background locations, even if they
# have values already included in the background
new_data <- addSamplesToBg(data, all = TRUE)
new_data



cleanEx()
nameEx("aicc")
### * aicc

flush(stderr()); flush(stdout())

### Name: aicc
### Title: AICc
### Aliases: aicc

### ** Examples

# Acquire environmental variables
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)

# Prepare presence and background locations
p_coords <- virtualSp$presence
bg_coords <- virtualSp$background

# Create SWD object
data <- prepareSWD(species = "Virtual species", p = p_coords, a = bg_coords,
                   env = predictors, categorical = "biome")

# Train a model
model <- train(method = "Maxnet", data = data, fc = "l")

# Compute the AICc
aicc(model, predictors)




cleanEx()
nameEx("auc")
### * auc

flush(stderr()); flush(stdout())

### Name: auc
### Title: AUC
### Aliases: auc

### ** Examples

# Acquire environmental variables
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)

# Prepare presence and background locations
p_coords <- virtualSp$presence
bg_coords <- virtualSp$background

# Create SWD object
data <- prepareSWD(species = "Virtual species", p = p_coords, a = bg_coords,
                   env = predictors, categorical = "biome")

# Split presence locations in training (80%) and testing (20%) datasets
datasets <- trainValTest(data, test = 0.2, only_presence = TRUE)
train <- datasets[[1]]
test <- datasets[[2]]

# Train a model
model <- train(method = "Maxnet", data = train, fc = "l")

# Compute the training AUC
auc(model)

# Compute the testing AUC
auc(model, test = test)




cleanEx()
nameEx("confMatrix")
### * confMatrix

flush(stderr()); flush(stdout())

### Name: confMatrix
### Title: Confusion Matrix
### Aliases: confMatrix

### ** Examples

# Acquire environmental variables
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)

# Prepare presence and background locations
p_coords <- virtualSp$presence
bg_coords <- virtualSp$background

# Create SWD object
data <- prepareSWD(species = "Virtual species", p = p_coords, a = bg_coords,
                   env = predictors, categorical = "biome")

# Train a model
model <- train(method = "Maxnet", data = data, fc = "l")

# Get the confusion matrix for thresholds ranging from 0 to 1
cm <- confMatrix(model, type = "cloglog")
head(cm)
tail(cm)

# Get the confusion matrix for a specific threshold
confMatrix(model, type = "logistic", th = 0.6)



cleanEx()
nameEx("corVar")
### * corVar

flush(stderr()); flush(stdout())

### Name: corVar
### Title: Print Correlated Variables
### Aliases: corVar

### ** Examples

# Acquire environmental variables
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)

# Prepare background locations
bg_coords <- dismo::randomPoints(predictors, 10000)

# Create SWD object
bg <- prepareSWD(species = "Virtual species", a = bg_coords,
                 env = predictors, categorical = "biome")

# Get the correlation among all the environmental variables
corVar(bg, method = "spearman")

# Get the environmental variables that have a correlation greater or equal to
# the given threshold
corVar(bg, method = "pearson", cor_th = 0.8)



cleanEx()
nameEx("doJk")
### * doJk

flush(stderr()); flush(stdout())

### Name: doJk
### Title: Jackknife Test
### Aliases: doJk

### ** Examples




cleanEx()
nameEx("get_tunable_args")
### * get_tunable_args

flush(stderr()); flush(stdout())

### Name: get_tunable_args
### Title: Get Tunable Arguments
### Aliases: get_tunable_args

### ** Examples

# Acquire environmental variables
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)

# Prepare presence and background locations
p_coords <- virtualSp$presence
bg_coords <- virtualSp$background

# Create SWD object
data <- prepareSWD(species = "Virtual species", p = p_coords, a = bg_coords,
                   env = predictors, categorical = "biome")

# Train a Maxent model and get tunable hyperparameters
model <- train(method = "Maxnet", data = data, fc = "l")
get_tunable_args(model)




cleanEx()
nameEx("gridSearch")
### * gridSearch

flush(stderr()); flush(stdout())

### Name: gridSearch
### Title: Grid Search
### Aliases: gridSearch

### ** Examples




cleanEx()
nameEx("maxentTh")
### * maxentTh

flush(stderr()); flush(stdout())

### Name: maxentTh
### Title: MaxEnt Thresholds
### Aliases: maxentTh

### ** Examples




cleanEx()
nameEx("maxentVarImp")
### * maxentVarImp

flush(stderr()); flush(stdout())

### Name: maxentVarImp
### Title: Maxent Variable Importance
### Aliases: maxentVarImp

### ** Examples




cleanEx()
nameEx("mergeSWD")
### * mergeSWD

flush(stderr()); flush(stdout())

### Name: mergeSWD
### Title: Merge SWD Objects
### Aliases: mergeSWD

### ** Examples

# Acquire environmental variables
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)

# Prepare presence and background locations
p_coords <- virtualSp$presence
bg_coords <- virtualSp$background

# Create SWD object
data <- prepareSWD(species = "Virtual species", p = p_coords, a = bg_coords,
                   env = predictors, categorical = "biome")

# Split only presence locations in training (80%) and testing (20%) datasets
datasets <- trainValTest(data, test = 0.2, only_presence = TRUE)
train <- datasets[[1]]
test <- datasets[[2]]

# Merge the training and the testing datasets together
merged <- mergeSWD(train, test, only_presence = TRUE)

# Split presence and absence locations in training (80%) and testing (20%)
datasets
datasets <- trainValTest(data, test = 0.2)
train <- datasets[[1]]
test <- datasets[[2]]

# Merge the training and the testing datasets together
merged <- mergeSWD(train, test)



cleanEx()
nameEx("modelReport")
### * modelReport

flush(stderr()); flush(stdout())

### Name: modelReport
### Title: Model Report
### Aliases: modelReport

### ** Examples




cleanEx()
nameEx("optimizeModel")
### * optimizeModel

flush(stderr()); flush(stdout())

### Name: optimizeModel
### Title: Optimize Model
### Aliases: optimizeModel

### ** Examples




cleanEx()
nameEx("plot-methods")
### * plot-methods

flush(stderr()); flush(stdout())

### Name: plot,SDMtune,missing-method
### Title: Plot SDMtune object
### Aliases: plot,SDMtune,missing-method

### ** Examples




cleanEx()
nameEx("plotCor")
### * plotCor

flush(stderr()); flush(stdout())

### Name: plotCor
### Title: Plot Correlation
### Aliases: plotCor

### ** Examples




cleanEx()
nameEx("plotJk")
### * plotJk

flush(stderr()); flush(stdout())

### Name: plotJk
### Title: Plot Jackknife Test
### Aliases: plotJk

### ** Examples




cleanEx()
nameEx("plotPA")
### * plotPA

flush(stderr()); flush(stdout())

### Name: plotPA
### Title: Plot Presence Absence Map
### Aliases: plotPA

### ** Examples




cleanEx()
nameEx("plotPred")
### * plotPred

flush(stderr()); flush(stdout())

### Name: plotPred
### Title: Plot Prediction
### Aliases: plotPred

### ** Examples




cleanEx()
nameEx("plotROC")
### * plotROC

flush(stderr()); flush(stdout())

### Name: plotROC
### Title: Plot ROC curve
### Aliases: plotROC

### ** Examples





cleanEx()
nameEx("plotResponse")
### * plotResponse

flush(stderr()); flush(stdout())

### Name: plotResponse
### Title: Plot Response Curve
### Aliases: plotResponse

### ** Examples




cleanEx()
nameEx("plotVarImp")
### * plotVarImp

flush(stderr()); flush(stdout())

### Name: plotVarImp
### Title: Plot Variable Importance
### Aliases: plotVarImp

### ** Examples




cleanEx()
nameEx("predict-SDMmodel-method")
### * predict-SDMmodel-method

flush(stderr()); flush(stdout())

### Name: predict,SDMmodel-method
### Title: Predict
### Aliases: predict,SDMmodel-method

### ** Examples

# Acquire environmental variables
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)

# Prepare presence and background locations
p_coords <- virtualSp$presence
bg_coords <- virtualSp$background

# Create SWD object
data <- prepareSWD(species = "Virtual species", p = p_coords, a = bg_coords,
                   env = predictors, categorical = "biome")

# Split presence locations in training (80%) and testing (20%) datasets
datasets <- trainValTest(data, test = 0.2, only_presence = TRUE)
train <- datasets[[1]]
test <- datasets[[2]]

# Train a model
model <- train(method = "Maxnet", data = train, fc = "l")

# Make cloglog prediction for the test dataset
predict(model, data = test, type = "cloglog")

# Make logistic prediction for the all study area
predict(model, data = predictors, type = "logistic")




cleanEx()
nameEx("predict-SDMmodelCV-method")
### * predict-SDMmodelCV-method

flush(stderr()); flush(stdout())

### Name: predict,SDMmodelCV-method
### Title: Predict for Cross Validation
### Aliases: predict,SDMmodelCV-method

### ** Examples




cleanEx()
nameEx("prepareSWD")
### * prepareSWD

flush(stderr()); flush(stdout())

### Name: prepareSWD
### Title: Prepare an SWD object
### Aliases: prepareSWD

### ** Examples

# Acquire environmental variables
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)

# Prepare presence and background locations
p_coords <- virtualSp$presence
bg_coords <- virtualSp$background

# Create the SWD object
data <- prepareSWD(species = "Virtual species", p = p_coords, a = bg_coords,
                   env = predictors, categorical = "biome")
data



cleanEx()
nameEx("randomFolds")
### * randomFolds

flush(stderr()); flush(stdout())

### Name: randomFolds
### Title: Create Random Folds
### Aliases: randomFolds

### ** Examples

# Acquire environmental variables
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)

# Prepare presence and background locations
p_coords <- virtualSp$presence
bg_coords <- virtualSp$background

data <- prepareSWD(species = "Virtual species", p = p_coords, a = bg_coords,
                   env = predictors, categorical = "biome")

# Create 4 random folds splitting presence and absence locations
folds <- randomFolds(data, k = 4)

# Create 4 random folds splitting only the presence locations
folds <- randomFolds(data, k = 4, only_presence = TRUE)



cleanEx()
nameEx("randomSearch")
### * randomSearch

flush(stderr()); flush(stdout())

### Name: randomSearch
### Title: Random Search
### Aliases: randomSearch

### ** Examples




cleanEx()
nameEx("reduceVar")
### * reduceVar

flush(stderr()); flush(stdout())

### Name: reduceVar
### Title: Reduce Variables
### Aliases: reduceVar

### ** Examples




cleanEx()
nameEx("swd2csv")
### * swd2csv

flush(stderr()); flush(stdout())

### Name: swd2csv
### Title: SWD to csv
### Aliases: swd2csv

### ** Examples




cleanEx()
nameEx("thinData")
### * thinData

flush(stderr()); flush(stdout())

### Name: thinData
### Title: Thin Data
### Aliases: thinData

### ** Examples




cleanEx()
nameEx("thresholds")
### * thresholds

flush(stderr()); flush(stdout())

### Name: thresholds
### Title: Thresholds
### Aliases: thresholds

### ** Examples

# Acquire environmental variables
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)

# Prepare presence and background locations
p_coords <- virtualSp$presence
bg_coords <- virtualSp$background

# Create SWD object
data <- prepareSWD(species = "Virtual species", p = p_coords, a = bg_coords,
                   env = predictors, categorical = "biome")

# Split presence locations in training (80%) and testing (20%) datasets
datasets <- trainValTest(data, test = 0.2, only_presence = TRUE)
train <- datasets[[1]]
test <- datasets[[2]]

# Train a model
model <- train(method = "Maxnet", data = train, fc = "l")

# Get the cloglog thresholds
thresholds(model, type = "cloglog")

# Get the logistic thresholds passing the test dataset
thresholds(model, type = "logistic", test = test)



cleanEx()
nameEx("train")
### * train

flush(stderr()); flush(stdout())

### Name: train
### Title: Train
### Aliases: train

### ** Examples




cleanEx()
nameEx("trainValTest")
### * trainValTest

flush(stderr()); flush(stdout())

### Name: trainValTest
### Title: Train, Validation and Test datasets
### Aliases: trainValTest

### ** Examples

# Acquire environmental variables
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)

# Prepare presence and background locations
p_coords <- virtualSp$presence
bg_coords <- virtualSp$background

# Create SWD object
data <- prepareSWD(species = "Virtual species", p = p_coords, a = bg_coords,
                   env = predictors, categorical = "biome")

# Split presence locations in training (80%) and testing (20%) datasets
# and splitting only the presence locations
datasets <- trainValTest(data, test = 0.2, only_presence = TRUE)
train <- datasets[[1]]
test <- datasets[[2]]

# Split presence locations in training (60%), validation (20%) and testing
# (20%) datasets and splitting the presence and the absence locations
datasets <- trainValTest(data, val = 0.2, test = 0.2)
train <- datasets[[1]]
val <- datasets[[2]]
test <- datasets[[3]]



cleanEx()
nameEx("tss")
### * tss

flush(stderr()); flush(stdout())

### Name: tss
### Title: True Skill Statistics
### Aliases: tss

### ** Examples

# Acquire environmental variables
files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
                    pattern = "grd", full.names = TRUE)
predictors <- raster::stack(files)

# Prepare presence and background locations
p_coords <- virtualSp$presence
bg_coords <- virtualSp$background

# Create SWD object
data <- prepareSWD(species = "Virtual species", p = p_coords, a = bg_coords,
                   env = predictors, categorical = "biome")

# Split presence locations in training (80%) and testing (20%) datasets
datasets <- trainValTest(data, test = 0.2, only_presence = TRUE)
train <- datasets[[1]]
test <- datasets[[2]]

# Train a model
model <- train(method = "Maxnet", data = train, fc = "l")

# Compute the training TSS
tss(model)

# Compute the testing TSS
tss(model, test)




cleanEx()
nameEx("varImp")
### * varImp

flush(stderr()); flush(stdout())

### Name: varImp
### Title: Variable Importance
### Aliases: varImp

### ** Examples




cleanEx()
nameEx("varSel")
### * varSel

flush(stderr()); flush(stdout())

### Name: varSel
### Title: Variable Selection
### Aliases: varSel

### ** Examples




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
