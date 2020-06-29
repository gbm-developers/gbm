pkgname <- "triplot"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('triplot')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("aspect_importance")
### * aspect_importance

flush(stderr()); flush(stdout())

### Name: aspect_importance
### Title: Calculates importance of variable groups (called aspects) for a
###   selected observation
### Aliases: aspect_importance aspect_importance.explainer
###   aspect_importance.default lime predict_aspects

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

predict_aspects(explain_titanic_glm,
                  new_observation = titanic_imputed[1,],
                  variable_groups = aspects)





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

explainer_titanic <- explain(model_titanic_glm,
                             data = titanic_imputed[,-8],
                             verbose = FALSE)
aspect_importance_single(explainer_titanic,
                         new_observation = titanic_imputed[1,-8])




cleanEx()
nameEx("calculate_triplot")
### * calculate_triplot

flush(stderr()); flush(stdout())

### Name: calculate_triplot
### Title: Calculate triplot that sums up automatic aspect/feature
###   importance grouping
### Aliases: calculate_triplot calculate_triplot.explainer
###   calculate_triplot.default print.triplot model_triplot predict_triplot

### ** Examples


library(DALEX)
set.seed(123)
apartments_num <- apartments[,unlist(lapply(apartments, is.numeric))]
apartments_num_lm_model <- lm(m2.price ~ ., data = apartments_num)
apartments_num_new_observation <- apartments_num[30, ]
explainer_apartments <- explain(model = apartments_num_lm_model,
                                data = apartments_num[,-1],
                                y = apartments_num[, 1],
                                verbose = FALSE)
apartments_tri <- calculate_triplot(x = explainer_apartments,
                                    new_observation =
                                      apartments_num_new_observation[-1])
apartments_tri
                                   




cleanEx()
nameEx("cluster_variables")
### * cluster_variables

flush(stderr()); flush(stdout())

### Name: cluster_variables
### Title: Creates a cluster tree from numeric features
### Aliases: cluster_variables cluster_variables.default

### ** Examples

library("DALEX")
dragons_data <- dragons[,c(2,3,4,7,8)]
cluster_variables(dragons_data, clust_method = "complete")




cleanEx()
nameEx("get_sample")
### * get_sample

flush(stderr()); flush(stdout())

### Name: get_sample
### Title: Function for getting binary matrix
### Aliases: get_sample

### ** Examples

 get_sample(100,6,"binom",3)




cleanEx()
nameEx("group_variables")
### * group_variables

flush(stderr()); flush(stdout())

### Name: group_variables
### Title: Helper function that combines clustering variables and creating
###   aspect list
### Aliases: group_variables

### ** Examples

library("DALEX")
dragons_data <- dragons[,c(2,3,4,7,8)]
group_variables(dragons_data, h = 0.5, clust_method = "complete")




cleanEx()
nameEx("hierarchical_importance")
### * hierarchical_importance

flush(stderr()); flush(stdout())

### Name: hierarchical_importance
### Title: Calculates importance of hierarchically grouped aspects
### Aliases: hierarchical_importance plot.hierarchical_importance

### ** Examples

library(DALEX)
apartments_num <- apartments[,unlist(lapply(apartments, is.numeric))]
apartments_num_lm_model <- lm(m2.price ~ ., data = apartments_num)
hi <- hierarchical_importance(x = apartments_num_lm_model,
 data = apartments_num[,-1],
 y = apartments_num[,1],
 type = "model")
plot(hi, add_last_group = TRUE, absolute_value = TRUE)




cleanEx()
nameEx("list_variables")
### * list_variables

flush(stderr()); flush(stdout())

### Name: list_variables
### Title: Cuts tree at custom height and returns a list
### Aliases: list_variables

### ** Examples

library("DALEX")
dragons_data <- dragons[,c(2,3,4,7,8)]
cv <- cluster_variables(dragons_data, clust_method = "complete")
list_variables(cv, h = 0.5)




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

titanic_ai <- predict_aspects(explain_titanic_glm,
                  new_observation = titanic_imputed[1,],
                  variable_groups = aspects)
plot(titanic_ai)




cleanEx()
nameEx("plot.cluster_variables")
### * plot.cluster_variables

flush(stderr()); flush(stdout())

### Name: plot.cluster_variables
### Title: Plots tree with correlation values
### Aliases: plot.cluster_variables

### ** Examples

library("DALEX")
dragons_data <- dragons[,c(2,3,4,7,8)]
cv <- cluster_variables(dragons_data, clust_method = "complete")
plot(cv, p = 0.7)




cleanEx()
nameEx("plot.triplot")
### * plot.triplot

flush(stderr()); flush(stdout())

### Name: plot.triplot
### Title: Plots triplot
### Aliases: plot.triplot

### ** Examples

library(DALEX)
set.seed(123)
apartments_num <- apartments[,unlist(lapply(apartments, is.numeric))]
apartments_num_lm_model <- lm(m2.price ~ ., data = apartments_num)
apartments_num_new_observation <- apartments_num[30, ]
explainer_apartments <- explain(model = apartments_num_lm_model,
                                data = apartments_num[,-1],
                                y = apartments_num[, 1],
                                verbose = FALSE)
apartments_tri <- calculate_triplot(x = explainer_apartments,
 new_observation = apartments_num_new_observation[-1])
plot(apartments_tri)




cleanEx()
nameEx("print.aspect_importance")
### * print.aspect_importance

flush(stderr()); flush(stdout())

### Name: print.aspect_importance
### Title: Function for printing aspect_importance results
### Aliases: print.aspect_importance

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

titanic_ai <- predict_aspects(explain_titanic_glm,
                  new_observation = titanic_imputed[1,],
                  variable_groups = aspects)
print(titanic_ai)




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
