pkgname <- "DALEXtra"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('DALEXtra')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("aspect_importance")
### * aspect_importance

flush(stderr()); flush(stdout())

### Name: aspect_importance
### Title: Calculates the feature groups importance (called aspects
###   importance) for a selected observation
### Aliases: aspect_importance aspect_importance.explainer
###   aspect_importance.default lime

### ** Examples

library("DALEX")

model_titanic_glm <- glm(survived == 1 ~
                         class+gender+age+sibsp+parch+fare+embarked,
                         data = titanic_imputed,
                         family = "binomial")

explain_titanic_glm <- explain(model_titanic_glm,
                               data = titanic_imputed[,-8],
                               y = titanic_imputed$survived == 1,
                               verbose = FALSE)

aspects <- list(wealth = c("class", "fare"),
                family = c("sibsp", "parch"),
                personal = c("gender", "age"),
                embarked = "embarked")

aspect_importance(explain_titanic_glm,
                  new_observation = titanic_imputed[1,],
                  aspects = aspects)





cleanEx()
nameEx("aspect_importance_single")
### * aspect_importance_single

flush(stderr()); flush(stdout())

### Name: aspect_importance_single
### Title: Aspects importance for single aspects
### Aliases: aspect_importance_single aspect_importance_single.explainer
###   aspect_importance_single.default

### ** Examples

library("DALEX")

model_titanic_glm <- glm(survived == 1 ~ class + gender + age +
                         sibsp + parch + fare + embarked,
                         data = titanic_imputed,
                         family = "binomial")

aspect_importance_single(model_titanic_glm, data = titanic_imputed[,-8],
                         new_observation = titanic_imputed[1,-8])




cleanEx()
nameEx("champion_challenger")
### * champion_challenger

flush(stderr()); flush(stdout())

### Name: champion_challenger
### Title: Compare machine learning models
### Aliases: champion_challenger

### ** Examples





cleanEx()
nameEx("create_env")
### * create_env

flush(stderr()); flush(stdout())

### Name: create_env
### Title: Create your conda virtual env with DALEX
### Aliases: create_env

### ** Examples




cleanEx()
nameEx("explain_h2o")
### * explain_h2o

flush(stderr()); flush(stdout())

### Name: explain_h2o
### Title: Create explainer from your h2o model
### Aliases: explain_h2o

### ** Examples




cleanEx()
nameEx("explain_keras")
### * explain_keras

flush(stderr()); flush(stdout())

### Name: explain_keras
### Title: Wrapper for Python Keras Models
### Aliases: explain_keras

### ** Examples

library("DALEXtra")




cleanEx()
nameEx("explain_mlr")
### * explain_mlr

flush(stderr()); flush(stdout())

### Name: explain_mlr
### Title: Create explainer from your mlr model
### Aliases: explain_mlr

### ** Examples

library("DALEXtra")
titanic_test <- read.csv(system.file("extdata", "titanic_test.csv", package = "DALEXtra"))
titanic_train <- read.csv(system.file("extdata", "titanic_train.csv", package = "DALEXtra"))
library("mlr")
task <- mlr::makeClassifTask(
id = "R",
data = titanic_train,
target = "survived"
)
learner <- mlr::makeLearner(
  "classif.gbm",
  par.vals = list(
    distribution = "bernoulli",
    n.trees = 500,
    interaction.depth = 4,
    n.minobsinnode = 12,
    shrinkage = 0.001,
    bag.fraction = 0.5,
    train.fraction = 1
  ),
  predict.type = "prob"
)
gbm <- mlr::train(learner, task)
explain_mlr(gbm, titanic_test[,1:17], titanic_test[,18])




cleanEx()
nameEx("explain_mlr3")
### * explain_mlr3

flush(stderr()); flush(stdout())

### Name: explain_mlr3
### Title: Create explainer from your mlr model
### Aliases: explain_mlr3

### ** Examples

library("DALEXtra")
library(mlr3)
titanic_imputed$survived <- as.factor(titanic_imputed$survived)
task_classif <- TaskClassif$new(id = "1", backend = titanic_imputed, target = "survived")
learner_classif <- lrn("classif.rpart", predict_type = "prob")
learner_classif$train(task_classif)
explain_mlr3(learner_classif, data = titanic_imputed,
             y = as.numeric(as.character(titanic_imputed$survived)))


task_regr <- TaskRegr$new(id = "2", backend = apartments, target = "m2.price")
learner_regr <- lrn("regr.rpart")
learner_regr$train(task_regr)
explain_mlr3(learner_regr, data = apartments, apartments$m2.price)




cleanEx()
nameEx("explain_scikitlearn")
### * explain_scikitlearn

flush(stderr()); flush(stdout())

### Name: explain_scikitlearn
### Title: Wrapper for Python Scikit-Learn Models
### Aliases: explain_scikitlearn

### ** Examples





cleanEx()
nameEx("funnel_measure")
### * funnel_measure

flush(stderr()); flush(stdout())

### Name: funnel_measure
### Title: Caluculate difference in performance in models across different
###   categories
### Aliases: funnel_measure

### ** Examples




cleanEx()
nameEx("get_sample")
### * get_sample

flush(stderr()); flush(stdout())

### Name: get_sample
### Title: Function for getting binary matrix
### Aliases: get_sample

### ** Examples

 ## Not run: 
##D  get_sample(100,6,"binom",3)
## End(Not run)



cleanEx()
nameEx("group_variables")
### * group_variables

flush(stderr()); flush(stdout())

### Name: group_variables
### Title: Groups numeric features into aspects
### Aliases: group_variables

### ** Examples

library("DALEX")
dragons_data <- dragons[,c(2,3,4,7,8)]
group_variables(dragons_data, p = 0.7, clust_method = "complete")




cleanEx()
nameEx("overall_comparison")
### * overall_comparison

flush(stderr()); flush(stdout())

### Name: overall_comparison
### Title: Compare champion with challengers globally
### Aliases: overall_comparison

### ** Examples




cleanEx()
nameEx("plot.aspect_importance")
### * plot.aspect_importance

flush(stderr()); flush(stdout())

### Name: plot.aspect_importance
### Title: Function for plotting aspect_importance results
### Aliases: plot.aspect_importance

### ** Examples

library("DALEX")

model_titanic_glm <- glm(survived == 1 ~
                         class+gender+age+sibsp+parch+fare+embarked,
                         data = titanic_imputed,
                         family = "binomial")

explain_titanic_glm <- explain(model_titanic_glm,
                               data = titanic_imputed[,-8],
                               y = titanic_imputed$survived == 1,
                               verbose = FALSE)

aspects <- list(wealth = c("class", "fare"),
                family = c("sibsp", "parch"),
                personal = c("gender", "age"),
                embarked = "embarked")

plot(aspect_importance(explain_titanic_glm,
                  new_observation = titanic_imputed[1,],
                  aspects = aspects))




cleanEx()
nameEx("plot.funnel_measure")
### * plot.funnel_measure

flush(stderr()); flush(stdout())

### Name: plot.funnel_measure
### Title: Funnel plot for difference in measures
### Aliases: plot.funnel_measure

### ** Examples




cleanEx()
nameEx("plot.overall_comparison")
### * plot.overall_comparison

flush(stderr()); flush(stdout())

### Name: plot.overall_comparison
### Title: Plot function for overall_comparison
### Aliases: plot.overall_comparison

### ** Examples




cleanEx()
nameEx("plot.trainig_test_comparison")
### * plot.trainig_test_comparison

flush(stderr()); flush(stdout())

### Name: plot.training_test_comparison
### Title: Plot and compare performance of model between training and test
###   set
### Aliases: plot.training_test_comparison

### ** Examples




cleanEx()
nameEx("plot_aspects_importance_grouping")
### * plot_aspects_importance_grouping

flush(stderr()); flush(stdout())

### Name: plot_aspects_importance_grouping
### Title: Function plots tree with aspect importance values
### Aliases: plot_aspects_importance_grouping

### ** Examples

library(DALEX)
apartments_num <- apartments[,unlist(lapply(apartments, is.numeric))]
apartments_num_lm_model <- lm(m2.price ~ ., data = apartments_num)
apartments_num_new_observation <- apartments_num[2,-1]
apartments_num_mod <- apartments_num[,-1]
plot_aspects_importance_grouping(x = apartments_num_lm_model,
data = apartments_num_mod, new_observation = apartments_num_new_observation)





cleanEx()
nameEx("plot_group_variables")
### * plot_group_variables

flush(stderr()); flush(stdout())

### Name: plot_group_variables
### Title: Plots tree with correlation values
### Aliases: plot_group_variables

### ** Examples

library("DALEX")
dragons_data <- dragons[,c(2,3,4,7,8)]
group_variables(dragons_data, p = 0.7, clust_method = "complete",
                draw_tree = TRUE)




cleanEx()
nameEx("print.funnel_measure")
### * print.funnel_measure

flush(stderr()); flush(stdout())

### Name: print.funnel_measure
### Title: Print funnel_measure object
### Aliases: print.funnel_measure

### ** Examples




cleanEx()
nameEx("print.overall_comparison")
### * print.overall_comparison

flush(stderr()); flush(stdout())

### Name: print.overall_comparison
### Title: Print overall_comparison object
### Aliases: print.overall_comparison

### ** Examples




cleanEx()
nameEx("print.training_test_comparison")
### * print.training_test_comparison

flush(stderr()); flush(stdout())

### Name: print.training_test_comparison
### Title: Print funnel_measure object
### Aliases: print.training_test_comparison

### ** Examples




cleanEx()
nameEx("trainig_test_comparison")
### * trainig_test_comparison

flush(stderr()); flush(stdout())

### Name: training_test_comparison
### Title: Compare performance of model between training and test set
### Aliases: training_test_comparison

### ** Examples




cleanEx()
nameEx("triplot")
### * triplot

flush(stderr()); flush(stdout())

### Name: triplot
### Title: Three plots that sum up automatic aspect importance grouping
### Aliases: triplot triplot.explainer triplot.default

### ** Examples

library(DALEX)
apartments_num <- apartments[,unlist(lapply(apartments, is.numeric))]
apartments_num_lm_model <- lm(m2.price ~ ., data = apartments_num)
apartments_num_new_observation <- apartments_num[30,-1]
apartments_num_mod <- apartments_num[,-1]
triplot(x = apartments_num_lm_model,
  data = apartments_num_mod,
  new_observation = apartments_num_new_observation,
  add_importance_labels = FALSE)





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
