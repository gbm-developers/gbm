pkgname <- "mlr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('mlr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("FailureModel")
### * FailureModel

flush(stderr()); flush(stdout())

### Name: FailureModel
### Title: Failure model.
### Aliases: FailureModel

### ** Examples

configureMlr(on.learner.error = "warn")
data = iris
data$newfeat = 1 # will make LDA crash
task = makeClassifTask(data = data, target = "Species")
m = train("classif.lda", task) # LDA crashed, but mlr catches this
print(m)
print(m$learner.model) # the error message
p = predict(m, task) # this will predict NAs
print(p)
print(performance(p))
configureMlr(on.learner.error = "stop")



cleanEx()
nameEx("Task")
### * Task

flush(stderr()); flush(stdout())

### Name: Task
### Title: Create a classification, regression, survival, cluster,
###   cost-sensitive classification or multilabel task.
### Aliases: Task

### ** Examples

if (requireNamespace("mlbench")) {
  library(mlbench)
  data(BostonHousing)
  data(Ionosphere)

  makeClassifTask(data = iris, target = "Species")
  makeRegrTask(data = BostonHousing, target = "medv")
  # an example of a classification task with more than those standard arguments:
  blocking = factor(c(rep(1, 51), rep(2, 300)))
  makeClassifTask(id = "myIonosphere", data = Ionosphere, target = "Class",
    positive = "good", blocking = blocking)
  makeClusterTask(data = iris[, -5L])
}



cleanEx()
nameEx("benchmark")
### * benchmark

flush(stderr()); flush(stdout())

### Name: benchmark
### Title: Benchmark experiment for multiple learners and tasks.
### Aliases: benchmark

### ** Examples

lrns = list(makeLearner("classif.lda"), makeLearner("classif.rpart"))
tasks = list(iris.task, sonar.task)
rdesc = makeResampleDesc("CV", iters = 2L)
meas = list(acc, ber)
bmr = benchmark(lrns, tasks, rdesc, measures = meas)
rmat = convertBMRToRankMatrix(bmr)
print(rmat)
plotBMRSummary(bmr)
plotBMRBoxplots(bmr, ber, style = "violin")
plotBMRRanksAsBarChart(bmr, pos = "stack")
friedmanTestBMR(bmr)
friedmanPostHocTestBMR(bmr, p.value = 0.05)



cleanEx()
nameEx("calculateConfusionMatrix")
### * calculateConfusionMatrix

flush(stderr()); flush(stdout())

### Name: calculateConfusionMatrix
### Title: Confusion matrix.
### Aliases: calculateConfusionMatrix print.ConfusionMatrix

### ** Examples

# get confusion matrix after simple manual prediction
allinds = 1:150
train = sample(allinds, 75)
test = setdiff(allinds, train)
mod = train("classif.lda", iris.task, subset = train)
pred = predict(mod, iris.task, subset = test)
print(calculateConfusionMatrix(pred))
print(calculateConfusionMatrix(pred, sums = TRUE))
print(calculateConfusionMatrix(pred, relative = TRUE))

# now after cross-validation
r = crossval("classif.lda", iris.task, iters = 2L)
print(calculateConfusionMatrix(r$pred))



cleanEx()
nameEx("calculateROCMeasures")
### * calculateROCMeasures

flush(stderr()); flush(stdout())

### Name: calculateROCMeasures
### Title: Calculate receiver operator measures.
### Aliases: calculateROCMeasures print.ROCMeasures

### ** Examples

lrn = makeLearner("classif.rpart", predict.type = "prob")
fit = train(lrn, sonar.task)
pred = predict(fit, task = sonar.task)
calculateROCMeasures(pred)



cleanEx()
nameEx("capLargeValues")
### * capLargeValues

flush(stderr()); flush(stdout())

### Name: capLargeValues
### Title: Convert large/infinite numeric values in a data.frame or task.
### Aliases: capLargeValues

### ** Examples

capLargeValues(iris, threshold = 5, impute = 5)



cleanEx()
nameEx("convertBMRToRankMatrix")
### * convertBMRToRankMatrix

flush(stderr()); flush(stdout())

### Name: convertBMRToRankMatrix
### Title: Convert BenchmarkResult to a rank-matrix.
### Aliases: convertBMRToRankMatrix

### ** Examples

# see benchmark



cleanEx()
nameEx("convertMLBenchObjToTask")
### * convertMLBenchObjToTask

flush(stderr()); flush(stdout())

### Name: convertMLBenchObjToTask
### Title: Convert a machine learning benchmark / demo object from package
###   mlbench to a task.
### Aliases: convertMLBenchObjToTask

### ** Examples

print(convertMLBenchObjToTask("Ionosphere"))
print(convertMLBenchObjToTask("mlbench.spirals", n = 100, sd = 0.1))



cleanEx()
nameEx("createSpatialResamplingPlots")
### * createSpatialResamplingPlots

flush(stderr()); flush(stdout())

### Name: createSpatialResamplingPlots
### Title: Create (spatial) resampling plot objects.
### Aliases: createSpatialResamplingPlots

### ** Examples




cleanEx()
nameEx("estimateRelativeOverfitting")
### * estimateRelativeOverfitting

flush(stderr()); flush(stdout())

### Name: estimateRelativeOverfitting
### Title: Estimate relative overfitting.
### Aliases: estimateRelativeOverfitting

### ** Examples

task = makeClassifTask(data = iris, target = "Species")
rdesc = makeResampleDesc("CV", iters = 2)
estimateRelativeOverfitting(rdesc, acc, task, makeLearner("classif.knn"))
estimateRelativeOverfitting(rdesc, acc, task, makeLearner("classif.lda"))
rpred = resample("classif.knn", task, rdesc)$pred
estimateRelativeOverfitting(rpred, acc, task)



cleanEx()
nameEx("extractFDAFeatures")
### * extractFDAFeatures

flush(stderr()); flush(stdout())

### Name: extractFDAFeatures
### Title: Extract features from functional data.
### Aliases: extractFDAFeatures

### ** Examples

df = data.frame(x = matrix(rnorm(24), ncol = 8), y = factor(c("a", "a", "b")))
fdf = makeFunctionalData(df, fd.features = list(x1 = 1:4, x2 = 5:8), exclude.cols = "y")
task = makeClassifTask(data = fdf, target = "y")
extracted = extractFDAFeatures(task,
  feat.methods = list("x1" = extractFDAFourier(), "x2" = extractFDAWavelets(filter = "haar")))
print(extracted$task)
reextractFDAFeatures(task, extracted$desc)



cleanEx()
nameEx("filterFeatures")
### * filterFeatures

flush(stderr()); flush(stdout())

### Name: filterFeatures
### Title: Filter features by thresholding filter values.
### Aliases: filterFeatures

### ** Examples

# simple filter
filterFeatures(iris.task, method = "FSelectorRcpp_gain.ratio", abs = 2)
# ensemble filter
filterFeatures(iris.task, method = "E-min",
  base.methods = c("FSelectorRcpp_gain.ratio",
    "FSelectorRcpp_information.gain"), abs = 2)



cleanEx()
nameEx("friedmanPostHocTestBMR")
### * friedmanPostHocTestBMR

flush(stderr()); flush(stdout())

### Name: friedmanPostHocTestBMR
### Title: Perform a posthoc Friedman-Nemenyi test.
### Aliases: friedmanPostHocTestBMR

### ** Examples

# see benchmark



cleanEx()
nameEx("friedmanTestBMR")
### * friedmanTestBMR

flush(stderr()); flush(stdout())

### Name: friedmanTestBMR
### Title: Perform overall Friedman test for a BenchmarkResult.
### Aliases: friedmanTestBMR

### ** Examples

# see benchmark



cleanEx()
nameEx("generateFeatureImportanceData")
### * generateFeatureImportanceData

flush(stderr()); flush(stdout())

### Name: generateFeatureImportanceData
### Title: Generate feature importance.
### Aliases: generateFeatureImportanceData FeatureImportanceData

### ** Examples


lrn = makeLearner("classif.rpart", predict.type = "prob")
fit = train(lrn, iris.task)
imp = generateFeatureImportanceData(iris.task, "permutation.importance",
  lrn, "Petal.Width", nmc = 10L, local = TRUE)



cleanEx()
nameEx("generateFilterValuesData")
### * generateFilterValuesData

flush(stderr()); flush(stdout())

### Name: generateFilterValuesData
### Title: Calculates feature filter values.
### Aliases: generateFilterValuesData FilterValues

### ** Examples

# two simple filter methods
fval = generateFilterValuesData(iris.task,
  method = c("FSelectorRcpp_gain.ratio", "FSelectorRcpp_information.gain"))
# using ensemble method "E-mean"
fval = generateFilterValuesData(iris.task,
  method = list("E-mean", c("FSelectorRcpp_gain.ratio",
    "FSelectorRcpp_information.gain")))



cleanEx()
nameEx("generateHyperParsEffectData")
### * generateHyperParsEffectData

flush(stderr()); flush(stdout())

### Name: generateHyperParsEffectData
### Title: Generate hyperparameter effect data.
### Aliases: generateHyperParsEffectData

### ** Examples

## Not run: 
##D # 3-fold cross validation
##D ps = makeParamSet(makeDiscreteParam("C", values = 2^(-4:4)))
##D ctrl = makeTuneControlGrid()
##D rdesc = makeResampleDesc("CV", iters = 3L)
##D res = tuneParams("classif.ksvm", task = pid.task, resampling = rdesc,
##D   par.set = ps, control = ctrl)
##D data = generateHyperParsEffectData(res)
##D plt = plotHyperParsEffect(data, x = "C", y = "mmce.test.mean")
##D plt + ylab("Misclassification Error")
##D 
##D # nested cross validation
##D ps = makeParamSet(makeDiscreteParam("C", values = 2^(-4:4)))
##D ctrl = makeTuneControlGrid()
##D rdesc = makeResampleDesc("CV", iters = 3L)
##D lrn = makeTuneWrapper("classif.ksvm", control = ctrl,
##D   resampling = rdesc, par.set = ps)
##D res = resample(lrn, task = pid.task, resampling = cv2,
##D   extract = getTuneResult)
##D data = generateHyperParsEffectData(res)
##D plotHyperParsEffect(data, x = "C", y = "mmce.test.mean", plot.type = "line")
## End(Not run)



cleanEx()
nameEx("generateLearningCurveData")
### * generateLearningCurveData

flush(stderr()); flush(stdout())

### Name: generateLearningCurveData
### Title: Generates a learning curve.
### Aliases: generateLearningCurveData LearningCurveData

### ** Examples

r = generateLearningCurveData(list("classif.rpart", "classif.knn"),
  task = sonar.task, percs = seq(0.2, 1, by = 0.2),
  measures = list(tp, fp, tn, fn),
  resampling = makeResampleDesc(method = "Subsample", iters = 5),
  show.info = FALSE)
plotLearningCurve(r)



cleanEx()
nameEx("generatePartialDependenceData")
### * generatePartialDependenceData

flush(stderr()); flush(stdout())

### Name: generatePartialDependenceData
### Title: Generate partial dependence.
### Aliases: generatePartialDependenceData PartialDependenceData

### ** Examples

lrn = makeLearner("regr.svm")
fit = train(lrn, bh.task)
pd = generatePartialDependenceData(fit, bh.task, "lstat")
plotPartialDependence(pd, data = getTaskData(bh.task))

lrn = makeLearner("classif.rpart", predict.type = "prob")
fit = train(lrn, iris.task)
pd = generatePartialDependenceData(fit, iris.task, "Petal.Width")
plotPartialDependence(pd, data = getTaskData(iris.task))



cleanEx()
nameEx("getCaretParamSet")
### * getCaretParamSet

flush(stderr()); flush(stdout())

### Name: getCaretParamSet
### Title: Get tuning parameters from a learner of the caret R-package.
### Aliases: getCaretParamSet

### ** Examples

if (requireNamespace("caret") && requireNamespace("mlbench")) {
  library(caret)
  classifTask = makeClassifTask(data = iris, target = "Species")

  # (1) classification (random forest) with discretized parameters
  getCaretParamSet("rf", length = 9L, task = classifTask, discretize = TRUE)

  # (2) regression (gradient boosting machine) without discretized parameters
  library(mlbench)
  data(BostonHousing)
  regrTask = makeRegrTask(data = BostonHousing, target = "medv")
  getCaretParamSet("gbm", length = 9L, task = regrTask, discretize = FALSE)
}



cleanEx()
nameEx("getHyperPars")
### * getHyperPars

flush(stderr()); flush(stdout())

### Name: getHyperPars
### Title: Get current parameter settings for a learner.
### Aliases: getHyperPars

### ** Examples

getHyperPars(makeLearner("classif.ranger"))

## set learner hyperparameter `mtry` manually
getHyperPars(makeLearner("classif.ranger", mtry = 100))



cleanEx()
nameEx("getMultilabelBinaryPerformances")
### * getMultilabelBinaryPerformances

flush(stderr()); flush(stdout())

### Name: getMultilabelBinaryPerformances
### Title: Retrieve binary classification measures for multilabel
###   classification predictions.
### Aliases: getMultilabelBinaryPerformances

### ** Examples

# see makeMultilabelBinaryRelevanceWrapper



cleanEx()
nameEx("getNestedTuneResultsOptPathDf")
### * getNestedTuneResultsOptPathDf

flush(stderr()); flush(stdout())

### Name: getNestedTuneResultsOptPathDf
### Title: Get the 'opt.path's from each tuning step from the outer
###   resampling.
### Aliases: getNestedTuneResultsOptPathDf

### ** Examples

# see example of makeTuneWrapper



cleanEx()
nameEx("getNestedTuneResultsX")
### * getNestedTuneResultsX

flush(stderr()); flush(stdout())

### Name: getNestedTuneResultsX
### Title: Get the tuned hyperparameter settings from a nested tuning.
### Aliases: getNestedTuneResultsX

### ** Examples

# see example of makeTuneWrapper



cleanEx()
nameEx("getOOBPreds")
### * getOOBPreds

flush(stderr()); flush(stdout())

### Name: getOOBPreds
### Title: Extracts out-of-bag predictions from trained models.
### Aliases: getOOBPreds

### ** Examples

training.set = sample(1:150, 50)
lrn = makeLearner("classif.ranger", predict.type = "prob", predict.threshold = 0.6)
mod = train(lrn, sonar.task, subset = training.set)
oob = getOOBPreds(mod, sonar.task)
oob
performance(oob, measures = list(auc, mmce))



cleanEx()
nameEx("getPredictionProbabilities")
### * getPredictionProbabilities

flush(stderr()); flush(stdout())

### Name: getPredictionProbabilities
### Title: Get probabilities for some classes.
### Aliases: getPredictionProbabilities

### ** Examples

task = makeClassifTask(data = iris, target = "Species")
lrn = makeLearner("classif.lda", predict.type = "prob")
mod = train(lrn, task)
# predict probabilities
pred = predict(mod, newdata = iris)

# Get probabilities for all classes
head(getPredictionProbabilities(pred))

# Get probabilities for a subset of classes
head(getPredictionProbabilities(pred, c("setosa", "virginica")))



cleanEx()
nameEx("getResamplingIndices")
### * getResamplingIndices

flush(stderr()); flush(stdout())

### Name: getResamplingIndices
### Title: Get the resampling indices from a tuning or feature selection
###   wrapper..
### Aliases: getResamplingIndices

### ** Examples

task = makeClassifTask(data = iris, target = "Species")
lrn = makeLearner("classif.rpart")
# stupid mini grid
ps = makeParamSet(
  makeDiscreteParam("cp", values = c(0.05, 0.1)),
  makeDiscreteParam("minsplit", values = c(10, 20))
)
ctrl = makeTuneControlGrid()
inner = makeResampleDesc("Holdout")
outer = makeResampleDesc("CV", iters = 2)
lrn = makeTuneWrapper(lrn, resampling = inner, par.set = ps, control = ctrl)
# nested resampling for evaluation
# we also extract tuned hyper pars in each iteration and by that the resampling indices
r = resample(lrn, task, outer, extract = getTuneResult)
# get tuning indices
getResamplingIndices(r, inner = TRUE)



cleanEx()
nameEx("getTaskData")
### * getTaskData

flush(stderr()); flush(stdout())

### Name: getTaskData
### Title: Extract data in task.
### Aliases: getTaskData

### ** Examples

library("mlbench")
data(BreastCancer)

df = BreastCancer
df$Id = NULL
task = makeClassifTask(id = "BreastCancer", data = df, target = "Class", positive = "malignant")
head(getTaskData)
head(getTaskData(task, features = c("Cell.size", "Cell.shape"), recode.target = "-1+1"))
head(getTaskData(task, subset = 1:100, recode.target = "01"))



cleanEx()
nameEx("getTaskTargets")
### * getTaskTargets

flush(stderr()); flush(stdout())

### Name: getTaskTargets
### Title: Get target data of task.
### Aliases: getTaskTargets

### ** Examples

task = makeClassifTask(data = iris, target = "Species")
getTaskTargets(task)



cleanEx()
nameEx("impute")
### * impute

flush(stderr()); flush(stdout())

### Name: impute
### Title: Impute and re-impute data
### Aliases: impute

### ** Examples

df = data.frame(x = c(1, 1, NA), y = factor(c("a", "a", "b")), z = 1:3)
imputed = impute(df, target = character(0), cols = list(x = 99, y = imputeMode()))
print(imputed$data)
reimpute(data.frame(x = NA_real_), imputed$desc)



cleanEx()
nameEx("joinClassLevels")
### * joinClassLevels

flush(stderr()); flush(stdout())

### Name: joinClassLevels
### Title: Join some class existing levels to new, larger class levels for
###   classification problems.
### Aliases: joinClassLevels

### ** Examples

joinClassLevels(iris.task, new.levels = list(foo = c("setosa", "virginica")))



cleanEx()
nameEx("listLearners")
### * listLearners

flush(stderr()); flush(stdout())

### Name: listLearners
### Title: Find matching learning algorithms.
### Aliases: listLearners listLearners.default listLearners.character
###   listLearners.Task

### ** Examples

## Not run: 
##D listLearners("classif", properties = c("multiclass", "prob"))
##D data = iris
##D task = makeClassifTask(data = data, target = "Species")
##D listLearners(task)
## End(Not run)



cleanEx()
nameEx("makeAggregation")
### * makeAggregation

flush(stderr()); flush(stdout())

### Name: makeAggregation
### Title: Specify your own aggregation of measures.
### Aliases: makeAggregation

### ** Examples

# computes the interquartile range on all performance values
test.iqr = makeAggregation(
  id = "test.iqr", name = "Test set interquartile range",
  properties = "req.test",
  fun = function(task, perf.test, perf.train, measure, group, pred) IQR(perf.test)
)



cleanEx()
nameEx("makeClassificationViaRegressionWrapper")
### * makeClassificationViaRegressionWrapper

flush(stderr()); flush(stdout())

### Name: makeClassificationViaRegressionWrapper
### Title: Classification via regression wrapper.
### Aliases: makeClassificationViaRegressionWrapper

### ** Examples

lrn = makeLearner("regr.rpart")
lrn = makeClassificationViaRegressionWrapper(lrn)
mod = train(lrn, sonar.task, subset = 1:140)
predictions = predict(mod, newdata = getTaskData(sonar.task)[141:208, 1:60])



cleanEx()
nameEx("makeFeatSelWrapper")
### * makeFeatSelWrapper

flush(stderr()); flush(stdout())

### Name: makeFeatSelWrapper
### Title: Fuse learner with feature selection.
### Aliases: makeFeatSelWrapper

### ** Examples

# nested resampling with feature selection (with a nonsense algorithm for selection)
outer = makeResampleDesc("CV", iters = 2L)
inner = makeResampleDesc("Holdout")
ctrl = makeFeatSelControlRandom(maxit = 1)
lrn = makeFeatSelWrapper("classif.ksvm", resampling = inner, control = ctrl)
# we also extract the selected features for all iteration here
r = resample(lrn, iris.task, outer, extract = getFeatSelResult)



cleanEx()
nameEx("makeFilterWrapper")
### * makeFilterWrapper

flush(stderr()); flush(stdout())

### Name: makeFilterWrapper
### Title: Fuse learner with a feature filter method.
### Aliases: makeFilterWrapper

### ** Examples




cleanEx()
nameEx("makeFunctionalData")
### * makeFunctionalData

flush(stderr()); flush(stdout())

### Name: makeFunctionalData
### Title: Create a data.frame containing functional features from a normal
###   data.frame.
### Aliases: makeFunctionalData

### ** Examples

# data.frame where columns 1:6 and 8:10 belong to a functional feature
d1 = data.frame(matrix(rnorm(100), nrow = 10), "target" = seq_len(10))
# Transform to functional data
d2 = makeFunctionalData(d1, fd.features = list("fd1" = 1:6, "fd2" = 8:10))
# Create a regression task
makeRegrTask(data = d2, target = "target")



cleanEx()
nameEx("makeLearner")
### * makeLearner

flush(stderr()); flush(stdout())

### Name: makeLearner
### Title: Create learner object.
### Aliases: makeLearner Learner

### ** Examples

makeLearner("classif.rpart")
makeLearner("classif.lda", predict.type = "prob")
lrn = makeLearner("classif.lda", method = "t", nu = 10)
getHyperPars(lrn)



cleanEx()
nameEx("makeLearners")
### * makeLearners

flush(stderr()); flush(stdout())

### Name: makeLearners
### Title: Create multiple learners at once.
### Aliases: makeLearners

### ** Examples

makeLearners(c("rpart", "lda"), type = "classif", predict.type = "prob")



cleanEx()
nameEx("makeMeasure")
### * makeMeasure

flush(stderr()); flush(stdout())

### Name: makeMeasure
### Title: Construct performance measure.
### Aliases: makeMeasure Measure

### ** Examples

f = function(task, model, pred, extra.args) {
  sum((pred$data$response - pred$data$truth)^2)
}
makeMeasure(id = "my.sse", minimize = TRUE,
  properties = c("regr", "response"), fun = f)



cleanEx()
nameEx("makeModelMultiplexer")
### * makeModelMultiplexer

flush(stderr()); flush(stdout())

### Name: makeModelMultiplexer
### Title: Create model multiplexer for model selection to tune over
###   multiple possible models.
### Aliases: makeModelMultiplexer ModelMultiplexer

### ** Examples

set.seed(123)



cleanEx()
nameEx("makeModelMultiplexerParamSet")
### * makeModelMultiplexerParamSet

flush(stderr()); flush(stdout())

### Name: makeModelMultiplexerParamSet
### Title: Creates a parameter set for model multiplexer tuning.
### Aliases: makeModelMultiplexerParamSet

### ** Examples

# See makeModelMultiplexer



cleanEx()
nameEx("makeMultilabelBinaryRelevanceWrapper")
### * makeMultilabelBinaryRelevanceWrapper

flush(stderr()); flush(stdout())

### Name: makeMultilabelBinaryRelevanceWrapper
### Title: Use binary relevance method to create a multilabel learner.
### Aliases: makeMultilabelBinaryRelevanceWrapper

### ** Examples

d = getTaskData(yeast.task)
# drop some labels so example runs faster
d = d[seq(1, nrow(d), by = 20), c(1:2, 15:17)]
task = makeMultilabelTask(data = d, target = c("label1", "label2"))
lrn = makeLearner("classif.rpart")
lrn = makeMultilabelBinaryRelevanceWrapper(lrn)
lrn = setPredictType(lrn, "prob")
# train, predict and evaluate
mod = train(lrn, task)
pred = predict(mod, task)
performance(pred, measure = list(multilabel.hamloss, multilabel.subset01, multilabel.f1))
# the next call basically has the same structure for any multilabel meta wrapper
getMultilabelBinaryPerformances(pred, measures = list(mmce, auc))
# above works also with predictions from resample!




cleanEx()
nameEx("makeMultilabelClassifierChainsWrapper")
### * makeMultilabelClassifierChainsWrapper

flush(stderr()); flush(stdout())

### Name: makeMultilabelClassifierChainsWrapper
### Title: Use classifier chains method (CC) to create a multilabel
###   learner.
### Aliases: makeMultilabelClassifierChainsWrapper

### ** Examples

d = getTaskData(yeast.task)
# drop some labels so example runs faster
d = d[seq(1, nrow(d), by = 20), c(1:2, 15:17)]
task = makeMultilabelTask(data = d, target = c("label1", "label2"))
lrn = makeLearner("classif.rpart")
lrn = makeMultilabelBinaryRelevanceWrapper(lrn)
lrn = setPredictType(lrn, "prob")
# train, predict and evaluate
mod = train(lrn, task)
pred = predict(mod, task)
performance(pred, measure = list(multilabel.hamloss, multilabel.subset01, multilabel.f1))
# the next call basically has the same structure for any multilabel meta wrapper
getMultilabelBinaryPerformances(pred, measures = list(mmce, auc))
# above works also with predictions from resample!




cleanEx()
nameEx("makeMultilabelDBRWrapper")
### * makeMultilabelDBRWrapper

flush(stderr()); flush(stdout())

### Name: makeMultilabelDBRWrapper
### Title: Use dependent binary relevance method (DBR) to create a
###   multilabel learner.
### Aliases: makeMultilabelDBRWrapper

### ** Examples

d = getTaskData(yeast.task)
# drop some labels so example runs faster
d = d[seq(1, nrow(d), by = 20), c(1:2, 15:17)]
task = makeMultilabelTask(data = d, target = c("label1", "label2"))
lrn = makeLearner("classif.rpart")
lrn = makeMultilabelBinaryRelevanceWrapper(lrn)
lrn = setPredictType(lrn, "prob")
# train, predict and evaluate
mod = train(lrn, task)
pred = predict(mod, task)
performance(pred, measure = list(multilabel.hamloss, multilabel.subset01, multilabel.f1))
# the next call basically has the same structure for any multilabel meta wrapper
getMultilabelBinaryPerformances(pred, measures = list(mmce, auc))
# above works also with predictions from resample!




cleanEx()
nameEx("makeMultilabelNestedStackingWrapper")
### * makeMultilabelNestedStackingWrapper

flush(stderr()); flush(stdout())

### Name: makeMultilabelNestedStackingWrapper
### Title: Use nested stacking method to create a multilabel learner.
### Aliases: makeMultilabelNestedStackingWrapper

### ** Examples

d = getTaskData(yeast.task)
# drop some labels so example runs faster
d = d[seq(1, nrow(d), by = 20), c(1:2, 15:17)]
task = makeMultilabelTask(data = d, target = c("label1", "label2"))
lrn = makeLearner("classif.rpart")
lrn = makeMultilabelBinaryRelevanceWrapper(lrn)
lrn = setPredictType(lrn, "prob")
# train, predict and evaluate
mod = train(lrn, task)
pred = predict(mod, task)
performance(pred, measure = list(multilabel.hamloss, multilabel.subset01, multilabel.f1))
# the next call basically has the same structure for any multilabel meta wrapper
getMultilabelBinaryPerformances(pred, measures = list(mmce, auc))
# above works also with predictions from resample!




cleanEx()
nameEx("makeMultilabelStackingWrapper")
### * makeMultilabelStackingWrapper

flush(stderr()); flush(stdout())

### Name: makeMultilabelStackingWrapper
### Title: Use stacking method (stacked generalization) to create a
###   multilabel learner.
### Aliases: makeMultilabelStackingWrapper

### ** Examples

d = getTaskData(yeast.task)
# drop some labels so example runs faster
d = d[seq(1, nrow(d), by = 20), c(1:2, 15:17)]
task = makeMultilabelTask(data = d, target = c("label1", "label2"))
lrn = makeLearner("classif.rpart")
lrn = makeMultilabelBinaryRelevanceWrapper(lrn)
lrn = setPredictType(lrn, "prob")
# train, predict and evaluate
mod = train(lrn, task)
pred = predict(mod, task)
performance(pred, measure = list(multilabel.hamloss, multilabel.subset01, multilabel.f1))
# the next call basically has the same structure for any multilabel meta wrapper
getMultilabelBinaryPerformances(pred, measures = list(mmce, auc))
# above works also with predictions from resample!




cleanEx()
nameEx("makeResampleDesc")
### * makeResampleDesc

flush(stderr()); flush(stdout())

### Name: makeResampleDesc
### Title: Create a description object for a resampling strategy.
### Aliases: makeResampleDesc ResampleDesc hout cv2 cv3 cv5 cv10

### ** Examples

# Bootstraping
makeResampleDesc("Bootstrap", iters = 10)
makeResampleDesc("Bootstrap", iters = 10, predict = "both")

# Subsampling
makeResampleDesc("Subsample", iters = 10, split = 3 / 4)
makeResampleDesc("Subsample", iters = 10)

# Holdout a.k.a. test sample estimation
makeResampleDesc("Holdout")



cleanEx()
nameEx("makeResampleInstance")
### * makeResampleInstance

flush(stderr()); flush(stdout())

### Name: makeResampleInstance
### Title: Instantiates a resampling strategy object.
### Aliases: makeResampleInstance ResampleInstance

### ** Examples

rdesc = makeResampleDesc("Bootstrap", iters = 10)
rin = makeResampleInstance(rdesc, task = iris.task)

rdesc = makeResampleDesc("CV", iters = 50)
rin = makeResampleInstance(rdesc, size = nrow(iris))

rin = makeResampleInstance("CV", iters = 10, task = iris.task)



cleanEx()
nameEx("makeStackedLearner")
### * makeStackedLearner

flush(stderr()); flush(stdout())

### Name: makeStackedLearner
### Title: Create a stacked learner object.
### Aliases: makeStackedLearner

### ** Examples

# Classification
data(iris)
tsk = makeClassifTask(data = iris, target = "Species")
base = c("classif.rpart", "classif.lda", "classif.svm")
lrns = lapply(base, makeLearner)
lrns = lapply(lrns, setPredictType, "prob")
m = makeStackedLearner(base.learners = lrns,
  predict.type = "prob", method = "hill.climb")
tmp = train(m, tsk)
res = predict(tmp, tsk)

# Regression
data(BostonHousing, package = "mlbench")
tsk = makeRegrTask(data = BostonHousing, target = "medv")
base = c("regr.rpart", "regr.svm")
lrns = lapply(base, makeLearner)
m = makeStackedLearner(base.learners = lrns,
  predict.type = "response", method = "compress")
tmp = train(m, tsk)
res = predict(tmp, tsk)



cleanEx()
nameEx("makeTuneWrapper")
### * makeTuneWrapper

flush(stderr()); flush(stdout())

### Name: makeTuneWrapper
### Title: Fuse learner with tuning.
### Aliases: makeTuneWrapper

### ** Examples




cleanEx()
nameEx("makeWeightedClassesWrapper")
### * makeWeightedClassesWrapper

flush(stderr()); flush(stdout())

### Name: makeWeightedClassesWrapper
### Title: Wraps a classifier for weighted fitting where each class
###   receives a weight.
### Aliases: makeWeightedClassesWrapper

### ** Examples




cleanEx()
nameEx("performance")
### * performance

flush(stderr()); flush(stdout())

### Name: performance
### Title: Measure performance of prediction.
### Aliases: performance

### ** Examples

training.set = seq(1, nrow(iris), by = 2)
test.set = seq(2, nrow(iris), by = 2)

task = makeClassifTask(data = iris, target = "Species")
lrn = makeLearner("classif.lda")
mod = train(lrn, task, subset = training.set)
pred = predict(mod, newdata = iris[test.set, ])
performance(pred, measures = mmce)

# Compute multiple performance measures at once
ms = list("mmce" = mmce, "acc" = acc, "timetrain" = timetrain)
performance(pred, measures = ms, task, mod)



cleanEx()
nameEx("plotBMRBoxplots")
### * plotBMRBoxplots

flush(stderr()); flush(stdout())

### Name: plotBMRBoxplots
### Title: Create box or violin plots for a BenchmarkResult.
### Aliases: plotBMRBoxplots

### ** Examples

# see benchmark



cleanEx()
nameEx("plotBMRRanksAsBarChart")
### * plotBMRRanksAsBarChart

flush(stderr()); flush(stdout())

### Name: plotBMRRanksAsBarChart
### Title: Create a bar chart for ranks in a BenchmarkResult.
### Aliases: plotBMRRanksAsBarChart

### ** Examples

# see benchmark



cleanEx()
nameEx("plotBMRSummary")
### * plotBMRSummary

flush(stderr()); flush(stdout())

### Name: plotBMRSummary
### Title: Plot a benchmark summary.
### Aliases: plotBMRSummary

### ** Examples

# see benchmark



cleanEx()
nameEx("plotCalibration")
### * plotCalibration

flush(stderr()); flush(stdout())

### Name: plotCalibration
### Title: Plot calibration data using ggplot2.
### Aliases: plotCalibration

### ** Examples

## Not run: 
##D lrns = list(makeLearner("classif.rpart", predict.type = "prob"),
##D   makeLearner("classif.nnet", predict.type = "prob"))
##D fit = lapply(lrns, train, task = iris.task)
##D pred = lapply(fit, predict, task = iris.task)
##D names(pred) = c("rpart", "nnet")
##D out = generateCalibrationData(pred, groups = 3)
##D plotCalibration(out)
##D 
##D fit = lapply(lrns, train, task = sonar.task)
##D pred = lapply(fit, predict, task = sonar.task)
##D names(pred) = c("rpart", "lda")
##D out = generateCalibrationData(pred)
##D plotCalibration(out)
## End(Not run)



cleanEx()
nameEx("plotCritDifferences")
### * plotCritDifferences

flush(stderr()); flush(stdout())

### Name: plotCritDifferences
### Title: Plot critical differences for a selected measure.
### Aliases: plotCritDifferences

### ** Examples

# see benchmark



cleanEx()
nameEx("plotFilterValues")
### * plotFilterValues

flush(stderr()); flush(stdout())

### Name: plotFilterValues
### Title: Plot filter values using ggplot2.
### Aliases: plotFilterValues

### ** Examples

fv = generateFilterValuesData(iris.task, method = "variance")
plotFilterValues(fv)



cleanEx()
nameEx("plotHyperParsEffect")
### * plotHyperParsEffect

flush(stderr()); flush(stdout())

### Name: plotHyperParsEffect
### Title: Plot the hyperparameter effects data
### Aliases: plotHyperParsEffect

### ** Examples

# see generateHyperParsEffectData



cleanEx()
nameEx("plotROCCurves")
### * plotROCCurves

flush(stderr()); flush(stdout())

### Name: plotROCCurves
### Title: Plots a ROC curve using ggplot2.
### Aliases: plotROCCurves

### ** Examples




cleanEx()
nameEx("plotThreshVsPerf")
### * plotThreshVsPerf

flush(stderr()); flush(stdout())

### Name: plotThreshVsPerf
### Title: Plot threshold vs. performance(s) for 2-class classification
###   using ggplot2.
### Aliases: plotThreshVsPerf

### ** Examples

lrn = makeLearner("classif.rpart", predict.type = "prob")
mod = train(lrn, sonar.task)
pred = predict(mod, sonar.task)
pvs = generateThreshVsPerfData(pred, list(acc, setAggregation(acc, train.mean)))
plotThreshVsPerf(pvs)



cleanEx()
nameEx("plotTuneMultiCritResult")
### * plotTuneMultiCritResult

flush(stderr()); flush(stdout())

### Name: plotTuneMultiCritResult
### Title: Plots multi-criteria results after tuning using ggplot2.
### Aliases: plotTuneMultiCritResult

### ** Examples

# see tuneParamsMultiCrit



cleanEx()
nameEx("predict.WrappedModel")
### * predict.WrappedModel

flush(stderr()); flush(stdout())

### Name: predict.WrappedModel
### Title: Predict new data.
### Aliases: predict.WrappedModel

### ** Examples

# train and predict
train.set = seq(1, 150, 2)
test.set = seq(2, 150, 2)
model = train("classif.lda", iris.task, subset = train.set)
p = predict(model, newdata = iris, subset = test.set)
print(p)
predict(model, task = iris.task, subset = test.set)

# predict now probabiliies instead of class labels
lrn = makeLearner("classif.lda", predict.type = "prob")
model = train(lrn, iris.task, subset = train.set)
p = predict(model, task = iris.task, subset = test.set)
print(p)
getPredictionProbabilities(p)



cleanEx()
nameEx("resample")
### * resample

flush(stderr()); flush(stdout())

### Name: resample
### Title: Fit models according to a resampling strategy.
### Aliases: resample crossval repcv holdout subsample bootstrapOOB
###   bootstrapB632 bootstrapB632plus growingcv fixedcv

### ** Examples

task = makeClassifTask(data = iris, target = "Species")
rdesc = makeResampleDesc("CV", iters = 2)
r = resample(makeLearner("classif.qda"), task, rdesc)
print(r$aggr)
print(r$measures.test)
print(r$pred)

# include the training set performance as well
rdesc = makeResampleDesc("CV", iters = 2, predict = "both")
r = resample(makeLearner("classif.qda"), task, rdesc,
  measures = list(mmce, setAggregation(mmce, train.mean)))
print(r$aggr)



cleanEx()
nameEx("selectFeatures")
### * selectFeatures

flush(stderr()); flush(stdout())

### Name: selectFeatures
### Title: Feature selection by wrapper approach.
### Aliases: selectFeatures

### ** Examples




cleanEx()
nameEx("setHyperPars")
### * setHyperPars

flush(stderr()); flush(stdout())

### Name: setHyperPars
### Title: Set the hyperparameters of a learner object.
### Aliases: setHyperPars

### ** Examples

cl1 = makeLearner("classif.ksvm", sigma = 1)
cl2 = setHyperPars(cl1, sigma = 10, par.vals = list(C = 2))
print(cl1)
# note the now set and altered hyperparameters:
print(cl2)



cleanEx()
nameEx("setThreshold")
### * setThreshold

flush(stderr()); flush(stdout())

### Name: setThreshold
### Title: Set threshold of prediction object.
### Aliases: setThreshold

### ** Examples

# create task and train learner (LDA)
task = makeClassifTask(data = iris, target = "Species")
lrn = makeLearner("classif.lda", predict.type = "prob")
mod = train(lrn, task)

# predict probabilities and compute performance
pred = predict(mod, newdata = iris)
performance(pred, measures = mmce)
head(as.data.frame(pred))

# adjust threshold and predict probabilities again
threshold = c(setosa = 0.4, versicolor = 0.3, virginica = 0.3)
pred = setThreshold(pred, threshold = threshold)
performance(pred, measures = mmce)
head(as.data.frame(pred))



cleanEx()
nameEx("subsetTask")
### * subsetTask

flush(stderr()); flush(stdout())

### Name: subsetTask
### Title: Subset data in task.
### Aliases: subsetTask

### ** Examples

task = makeClassifTask(data = iris, target = "Species")
subsetTask(task, subset = 1:100)



cleanEx()
nameEx("summarizeColumns")
### * summarizeColumns

flush(stderr()); flush(stdout())

### Name: summarizeColumns
### Title: Summarize columns of data.frame or task.
### Aliases: summarizeColumns

### ** Examples

summarizeColumns(iris)



cleanEx()
nameEx("summarizeLevels")
### * summarizeLevels

flush(stderr()); flush(stdout())

### Name: summarizeLevels
### Title: Summarizes factors of a data.frame by tabling them.
### Aliases: summarizeLevels

### ** Examples

summarizeLevels(iris)



cleanEx()
nameEx("train")
### * train

flush(stderr()); flush(stdout())

### Name: train
### Title: Train a learning algorithm.
### Aliases: train

### ** Examples

training.set = sample(seq_len(nrow(iris)), nrow(iris) / 2)

## use linear discriminant analysis to classify iris data
task = makeClassifTask(data = iris, target = "Species")
learner = makeLearner("classif.lda", method = "mle")
mod = train(learner, task, subset = training.set)
print(mod)

## use random forest to classify iris data
task = makeClassifTask(data = iris, target = "Species")
learner = makeLearner("classif.rpart", minsplit = 7, predict.type = "prob")
mod = train(learner, task, subset = training.set)
print(mod)



cleanEx()
nameEx("tuneParams")
### * tuneParams

flush(stderr()); flush(stdout())

### Name: tuneParams
### Title: Hyperparameter tuning.
### Aliases: tuneParams

### ** Examples

set.seed(123)
# a grid search for an SVM (with a tiny number of points...)
# note how easily we can optimize on a log-scale
ps = makeParamSet(
  makeNumericParam("C", lower = -12, upper = 12, trafo = function(x) 2^x),
  makeNumericParam("sigma", lower = -12, upper = 12, trafo = function(x) 2^x)
)
ctrl = makeTuneControlGrid(resolution = 2L)
rdesc = makeResampleDesc("CV", iters = 2L)
res = tuneParams("classif.ksvm", iris.task, rdesc, par.set = ps, control = ctrl)
print(res)
# access data for all evaluated points
df = as.data.frame(res$opt.path)
df1 = as.data.frame(res$opt.path, trafo = TRUE)
print(head(df[, -ncol(df)]))
print(head(df1[, -ncol(df)]))
# access data for all evaluated points - alternative
df2 = generateHyperParsEffectData(res)
df3 = generateHyperParsEffectData(res, trafo = TRUE)
print(head(df2$data[, -ncol(df2$data)]))
print(head(df3$data[, -ncol(df3$data)]))
## Not run: 
##D # we optimize the SVM over 3 kernels simultanously
##D # note how we use dependent params (requires = ...) and iterated F-racing here
##D ps = makeParamSet(
##D   makeNumericParam("C", lower = -12, upper = 12, trafo = function(x) 2^x),
##D   makeDiscreteParam("kernel", values = c("vanilladot", "polydot", "rbfdot")),
##D   makeNumericParam("sigma", lower = -12, upper = 12, trafo = function(x) 2^x,
##D     requires = quote(kernel == "rbfdot")),
##D   makeIntegerParam("degree", lower = 2L, upper = 5L,
##D     requires = quote(kernel == "polydot"))
##D )
##D print(ps)
##D ctrl = makeTuneControlIrace(maxExperiments = 5, nbIterations = 1, minNbSurvival = 1)
##D rdesc = makeResampleDesc("Holdout")
##D res = tuneParams("classif.ksvm", iris.task, rdesc, par.set = ps, control = ctrl)
##D print(res)
##D df = as.data.frame(res$opt.path)
##D print(head(df[, -ncol(df)]))
##D 
##D # include the training set performance as well
##D rdesc = makeResampleDesc("Holdout", predict = "both")
##D res = tuneParams("classif.ksvm", iris.task, rdesc, par.set = ps,
##D   control = ctrl, measures = list(mmce, setAggregation(mmce, train.mean)))
##D print(res)
##D df2 = as.data.frame(res$opt.path)
##D print(head(df2[, -ncol(df2)]))
## End(Not run)



cleanEx()
nameEx("tuneParamsMultiCrit")
### * tuneParamsMultiCrit

flush(stderr()); flush(stdout())

### Name: tuneParamsMultiCrit
### Title: Hyperparameter tuning for multiple measures at once.
### Aliases: tuneParamsMultiCrit

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
