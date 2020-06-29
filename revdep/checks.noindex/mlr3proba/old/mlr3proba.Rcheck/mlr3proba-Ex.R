pkgname <- "mlr3proba"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('mlr3proba')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("LearnerDens")
### * LearnerDens

flush(stderr()); flush(stdout())

### Name: LearnerDens
### Title: Density Learner
### Aliases: LearnerDens

### ** Examples

library(mlr3)
# get all density learners from mlr_learners:
lrns = mlr_learners$mget(mlr_learners$keys("^dens"))
names(lrns)

# get a specific learner from mlr_learners:
mlr_learners$get("dens.hist")
lrn("dens.hist")



cleanEx()
nameEx("LearnerSurv")
### * LearnerSurv

flush(stderr()); flush(stdout())

### Name: LearnerSurv
### Title: Survival Learner
### Aliases: LearnerSurv

### ** Examples

library(mlr3)
# get all survival learners from mlr_learners:
lrns = mlr_learners$mget(mlr_learners$keys("^surv"))
names(lrns)

# get a specific learner from mlr_learners:
mlr_learners$get("surv.coxph")
lrn("surv.coxph")



cleanEx()
nameEx("LearnerSurvBlackboost")
### * LearnerSurvBlackboost

flush(stderr()); flush(stdout())

### Name: LearnerSurvBlackboost
### Title: Gradient Boosting with Regression Trees Survival Learner
### Aliases: LearnerSurvBlackboost mlr_learners_surv.blackboost

### ** Examples

library(mlr3)
task = tgen("simsurv")$generate(20)
learner = lrn("surv.blackboost")
resampling = rsmp("cv", folds = 2)
resample(task, learner, resampling)



cleanEx()
nameEx("LearnerSurvGamboost")
### * LearnerSurvGamboost

flush(stderr()); flush(stdout())

### Name: LearnerSurvGamboost
### Title: Gradient Boosting for Additive Models Survival Learner
### Aliases: LearnerSurvGamboost mlr_learners_surv.gamboost

### ** Examples

library(mlr3)
task = tgen("simsurv")$generate(20)
learner = lrn("surv.gamboost")
learner$param_set$values = mlr3misc::insert_named(
  learner$param_set$values,
  list(dfbase = 3, center = TRUE, baselearner = "bols"))
resampling = rsmp("cv", folds = 2)
resample(task, learner, resampling)



cleanEx()
nameEx("LearnerSurvGlmboost")
### * LearnerSurvGlmboost

flush(stderr()); flush(stdout())

### Name: LearnerSurvGlmboost
### Title: Gradient Boosting with Component-wise Linear Models Survival
###   Learner
### Aliases: LearnerSurvGlmboost mlr_learners_surv.glmboost

### ** Examples

library(mlr3)
task = tgen("simsurv")$generate(20)
learner = lrn("surv.glmboost")
resampling = rsmp("cv", folds = 3)
resample(task, learner, resampling)



cleanEx()
nameEx("LearnerSurvMboost")
### * LearnerSurvMboost

flush(stderr()); flush(stdout())

### Name: LearnerSurvMboost
### Title: Gradient Boosting for Generalized Additive Models Survival
###   Learner
### Aliases: LearnerSurvMboost mlr_learners_surv.mboost

### ** Examples

library(mlr3)
task = tgen("simsurv")$generate(20)
learner = lrn("surv.mboost")
learner$param_set$values = mlr3misc::insert_named(
  learner$param_set$values,
  list(center = TRUE, baselearner = "bols"))
resampling = rsmp("cv", folds = 2)
resample(task, learner, resampling)



cleanEx()
nameEx("PipeOpCrankCompositor")
### * PipeOpCrankCompositor

flush(stderr()); flush(stdout())

### Name: PipeOpCrankCompositor
### Title: PipeOpCrankCompositor
### Aliases: PipeOpCrankCompositor mlr_pipeops_crankcompose

### ** Examples

library(mlr3)
library(mlr3pipelines)
set.seed(1)

# Three methods to predict a `crank` from `surv.rpart`
task = tgen("simsurv")$generate(30)

# Method 1 - Train and predict separately then compose
learn = lrn("surv.coxph")$train(task)$predict(task)
poc = po("crankcompose", param_vals = list(method = "mean"))
poc$predict(list(learn))

# Examples not run to save run-time.
## Not run: 
##D # Method 2 - Create a graph manually
##D gr = Graph$new()$
##D   add_pipeop(po("learner", lrn("surv.ranger")))$
##D   add_pipeop(po("crankcompose"))$
##D   add_edge("surv.ranger", "crankcompose")
##D gr$train(task)
##D gr$predict(task)
##D 
##D # Method 3 - Syntactic sugar: Wrap the learner in a graph
##D ranger.crank = crankcompositor(
##D   learner = lrn("surv.ranger"),
##D   method = "median")
##D resample(task, ranger.crank, rsmp("cv", folds = 2))$predictions()
## End(Not run)



cleanEx()
nameEx("PipeOpDistrCompositor")
### * PipeOpDistrCompositor

flush(stderr()); flush(stdout())

### Name: PipeOpDistrCompositor
### Title: PipeOpDistrCompositor
### Aliases: PipeOpDistrCompositor mlr_pipeops_distrcompose

### ** Examples

library(mlr3)
library(mlr3pipelines)
set.seed(42)

# Three methods to transform the cox ph predicted `distr` to an
#  accelerated failure time model
task = tgen("simsurv")$generate(30)

# Method 1 - Train and predict separately then compose
base = lrn("surv.kaplan")$train(task)$predict(task)
pred = lrn("surv.coxph")$train(task)$predict(task)
pod = po("distrcompose", param_vals = list(form = "aft", overwrite = TRUE))
pod$predict(list(base = base, pred = pred))

# Examples not run to save run-time.
## Not run: 
##D # Method 2 - Create a graph manually
##D gr = Graph$new()$
##D   add_pipeop(po("learner", lrn("surv.kaplan")))$
##D   add_pipeop(po("learner", lrn("surv.glmnet")))$
##D   add_pipeop(po("distrcompose"))$
##D   add_edge("surv.kaplan", "distrcompose", dst_channel = "base")$
##D   add_edge("surv.glmnet", "distrcompose", dst_channel = "pred")
##D gr$train(task)$gr$predict(task)
##D 
##D # Method 3 - Syntactic sugar: Wrap the learner in a graph.
##D cvglm.distr = distrcompositor(
##D   learner = lrn("surv.cvglmnet"),
##D   estimator = "kaplan",
##D   form = "aft")
##D cvglm.distr$fit(task)$predict(task)
## End(Not run)



cleanEx()
nameEx("PredictionDens")
### * PredictionDens

flush(stderr()); flush(stdout())

### Name: PredictionDens
### Title: Prediction Object for Density
### Aliases: PredictionDens

### ** Examples

library(mlr3)
task = mlr_tasks$get("precip")
learner = mlr_learners$get("dens.hist")
p = learner$train(task)$predict(task)
head(as.data.table(p))



cleanEx()
nameEx("PredictionSurv")
### * PredictionSurv

flush(stderr()); flush(stdout())

### Name: PredictionSurv
### Title: Prediction Object for Survival
### Aliases: PredictionSurv

### ** Examples

library(mlr3)
task = tgen("simsurv")$generate(20)
learner = mlr_learners$get("surv.rpart")
p = learner$train(task)$predict(task)
head(as.data.table(p))



cleanEx()
nameEx("TaskDens")
### * TaskDens

flush(stderr()); flush(stdout())

### Name: TaskDens
### Title: Density Task
### Aliases: TaskDens

### ** Examples

task = TaskDens$new("precip", backend = data.frame(target = precip), target = "target")
task$task_type
task$truth()



cleanEx()
nameEx("TaskSurv")
### * TaskSurv

flush(stderr()); flush(stdout())

### Name: TaskSurv
### Title: Survival Task
### Aliases: TaskSurv

### ** Examples

library(mlr3)
lung = mlr3misc::load_dataset("lung", package = "survival")
lung$status = (lung$status == 2L)
b = as_data_backend(lung)
task = TaskSurv$new("lung",
  backend = b, time = "time",
  event = "status")

task$target_names
task$feature_names
task$formula()
task$truth()



cleanEx()
nameEx("crankcompositor")
### * crankcompositor

flush(stderr()); flush(stdout())

### Name: crankcompositor
### Title: Compose a Crank Predict Type for Survival Learners
### Aliases: crankcompositor

### ** Examples

## Not run: 
##D library(mlr3)
##D library(mlr3pipelines)
##D 
##D task = tgen("simsurv")$generate(20)
##D ranger.crank = crankcompositor(
##D   learner = lrn("surv.coxph"),
##D   method = "median")
##D resample(task, ranger.crank, rsmp("cv", folds = 2))$predictions()
## End(Not run)



cleanEx()
nameEx("distrcompositor")
### * distrcompositor

flush(stderr()); flush(stdout())

### Name: distrcompositor
### Title: Compose a Distr Predict Type for Survival Learners
### Aliases: distrcompositor

### ** Examples

## Not run: 
##D library("mlr3")
##D library("mlr3pipelines")
##D 
##D task = tgen("simsurv")$generate(20)
##D cvglm.distr = distrcompositor(
##D   learner = lrn("surv.cvglmnet"),
##D   estimator = "kaplan",
##D   form = "aft")
##D 
##D resample(task, cvglm.distr, rsmp("cv", folds = 2))$predictions()
## End(Not run)



cleanEx()
nameEx("mlr_task_generators_simdens")
### * mlr_task_generators_simdens

flush(stderr()); flush(stdout())

### Name: mlr_task_generators_simdens
### Title: Density Task Generator for Package 'distr6'
### Aliases: mlr_task_generators_simdens TaskGeneratorSimdens

### ** Examples

generator = mlr3::mlr_task_generators$get("simdens")
task = generator$generate(20)
task$head()



cleanEx()
nameEx("mlr_task_generators_simsurv")
### * mlr_task_generators_simsurv

flush(stderr()); flush(stdout())

### Name: mlr_task_generators_simsurv
### Title: Survival Task Generator for Package 'simsurv'
### Aliases: mlr_task_generators_simsurv TaskGeneratorSimsurv

### ** Examples

generator = mlr3::mlr_task_generators$get("simsurv")
task = generator$generate(20)
task$head()



cleanEx()
nameEx("pecs")
### * pecs

flush(stderr()); flush(stdout())

### Name: pecs
### Title: Prediction Error Curves for PredictionSurv and LearnerSurv
### Aliases: pecs pecs.list pecs.PredictionSurv

### ** Examples

## Not run: 
##D library(mlr3)
##D task = tsk("rats")
##D 
##D # Prediction Error Curves for prediction object
##D learn = lrn("surv.coxph")
##D p = learn$train(task)$predict(task)
##D pecs(p)
##D pecs(p, measure = "logloss", times = c(20, 40, 60, 80)) +
##D   ggplot2::geom_point() +
##D   ggplot2::ggtitle("Logloss Prediction Error Curve for Cox PH")
##D 
##D # Access underlying data
##D x = pecs(p)
##D x$data
##D 
##D # Prediction Error Curves for fitted learners
##D learns = lrns(c("surv.kaplan", "surv.coxph", "surv.ranger"))
##D lapply(learns, function(x) x$train(task))
##D pecs(learns, task = task, measure = "logloss", times = c(20, 90), n = 10)
## End(Not run)




cleanEx()
nameEx("plot.LearnerSurv")
### * plot.LearnerSurv

flush(stderr()); flush(stdout())

### Name: plot.LearnerSurv
### Title: Visualization of fitted 'LearnerSurv' objects
### Aliases: plot.LearnerSurv

### ** Examples

library(mlr3)
task = tsk("rats")

# Prediction Error Curves for prediction object
learn = lrn("surv.coxph")
learn$train(task)

plot(learn, task, "survival", ind = 10)
plot(learn, task, "survival", row_ids = 1:5)
plot(learn, task, "survival", newdata = task$data()[1:5, ])
plot(learn, task, "survival", newdata = task$data()[1:5, ], ylim = c(0, 1))



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
