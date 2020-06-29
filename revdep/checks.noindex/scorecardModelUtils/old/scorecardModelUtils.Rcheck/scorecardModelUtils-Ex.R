pkgname <- "scorecardModelUtils"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('scorecardModelUtils')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("cat_new_class")
### * cat_new_class

flush(stderr()); flush(stdout())

### Name: cat_new_class
### Title: Clubbing class of categorical variables with low population
###   percentage with another class of similar event rate
### Aliases: cat_new_class

### ** Examples

data <- iris[1:110,]
data$Species <- as.character(data$Species)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
data_newclass <- cat_new_class(base = data,target = "Y",cat_var_name = "Species",threshold = 0.1)



cleanEx()
nameEx("categorical_iv")
### * categorical_iv

flush(stderr()); flush(stdout())

### Name: categorical_iv
### Title: IV table for individual categorical variable
### Aliases: categorical_iv

### ** Examples

data <- iris
data$Species <- as.character(data$Species)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
cat_iv <- categorical_iv(base = data,target = "Y",variable = "Species",event = 1)



cleanEx()
nameEx("club_cat_class")
### * club_cat_class

flush(stderr()); flush(stdout())

### Name: club_cat_class
### Title: Clubbing class of a categorical variable with low population
###   percentage with another class of similar event rate
### Aliases: club_cat_class

### ** Examples

data <- iris[1:110,]
data$Species <- as.character(data$Species)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
data_clubclass <- club_cat_class(base = data,target = "Y",variable = "Species",threshold = 0.2)



cleanEx()
nameEx("cv_filter")
### * cv_filter

flush(stderr()); flush(stdout())

### Name: cv_filter
### Title: Variable reduction based on Cramer's V filter
### Aliases: cv_filter

### ** Examples

data <- iris
suppressWarnings(RNGversion('3.5.0'))
set.seed(11)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
cv_tab_list <- cv_table(data, c("Species", "Sepal.Length"))
cv_tab <- cv_tab_list$cv_val_tab
x <- c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")
iv_table_list <- iv_table(base = data,target = "Y",num_var_name = x,cat_var_name = "Species")
iv_tab <- iv_table_list$iv_table
cv_filter_list <- cv_filter(cv_table = cv_tab,iv_table = iv_tab,threshold = 0.5)
cv_filter_list$retain_var_list
cv_filter_list$dropped_var_list
cv_filter_list$dropped_var_tab
cv_filter_list$threshold



cleanEx()
nameEx("cv_table")
### * cv_table

flush(stderr()); flush(stdout())

### Name: cv_table
### Title: Pairwise Cramer's V among a list of categorical variables
### Aliases: cv_table

### ** Examples

data <- iris
data$Species <- as.character(data$Species)
data$Sepal.Length <- as.character(floor(data$Sepal.Length))
cv_tab_list <- cv_table(data, c("Species", "Sepal.Length"))
cv_tab_list$cv_val_tab
cv_tab_list$single_class_var_index



cleanEx()
nameEx("cv_test")
### * cv_test

flush(stderr()); flush(stdout())

### Name: cv_test
### Title: Cramer's V value between two categorical variables
### Aliases: cv_test

### ** Examples

data <- iris
data$Species <- as.character(data$Species)
data$Sepal.Length <- as.character(floor(data$Sepal.Length))
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
cv_result <- cv_test(base = data,var_1 = "Species",var_2 = "Sepal.Length")



cleanEx()
nameEx("dtree_split_val")
### * dtree_split_val

flush(stderr()); flush(stdout())

### Name: dtree_split_val
### Title: Getting the split value for terminal nodes from decision tree
### Aliases: dtree_split_val

### ** Examples

data <- iris
data$Y <- ifelse(data$Species=="setosa",1,0)



cleanEx()
nameEx("dtree_trend_iv")
### * dtree_trend_iv

flush(stderr()); flush(stdout())

### Name: dtree_trend_iv
### Title: Recursive Decision Tree partitioning with monotonic event rate
###   along with IV table for individual numerical variable
### Aliases: dtree_trend_iv

### ** Examples

data <- iris
data$Y <- ifelse(data$Species=="setosa",1,0)
dtree_trend_tab <- dtree_trend_iv(base = data,target = "Y",variable = "Sepal.Length",event = 1)



cleanEx()
nameEx("fn_conf_mat")
### * fn_conf_mat

flush(stderr()); flush(stdout())

### Name: fn_conf_mat
### Title: Creates confusion matrix and its related measures
### Aliases: fn_conf_mat

### ** Examples

data <- iris
data$Species <- as.character(data$Species)
suppressWarnings(RNGversion('3.5.0'))
set.seed(11)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
data$Y_pred <- sample(0:1,size=nrow(data),replace=TRUE)
fn_conf_mat_list <- fn_conf_mat(base = data,observed_col = "Y",predicted_col = "Y_pred",event = 1)
fn_conf_mat_list$confusion_mat
fn_conf_mat_list$accuracy
fn_conf_mat_list$precision
fn_conf_mat_list$recall
fn_conf_mat_list$sensitivity
fn_conf_mat_list$specificity
fn_conf_mat_list$f1_score



cleanEx()
nameEx("fn_cross_index")
### * fn_cross_index

flush(stderr()); flush(stdout())

### Name: fn_cross_index
### Title: Creates random index for k-fold cross validation
### Aliases: fn_cross_index

### ** Examples

data <- iris
data$Species <- as.character(data$Species)
suppressWarnings(RNGversion('3.5.0'))
set.seed(11)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
data$Y_pred <- sample(0:1,size=nrow(data),replace=TRUE)
data_k_list <- fn_cross_index(base = data,k = 5)
data_k_list$index1
data_k_list$index2
data_k_list$index3
data_k_list$index4
data_k_list$index5



cleanEx()
nameEx("fn_error")
### * fn_error

flush(stderr()); flush(stdout())

### Name: fn_error
### Title: Computes error measures between observed and predicted values
### Aliases: fn_error

### ** Examples

data <- iris
data$Species <- as.character(data$Species)
suppressWarnings(RNGversion('3.5.0'))
set.seed(11)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
data$Y_pred <- sample(0:1,size=nrow(data),replace=TRUE)
fn_error_list <- fn_error(base = data,observed_col = "Y",predicted_col = "Y_pred")
fn_error_list$mean_abs_error
fn_error_list$mean_sq_error
fn_error_list$root_mean_sq_error



cleanEx()
nameEx("fn_mode")
### * fn_mode

flush(stderr()); flush(stdout())

### Name: fn_mode
### Title: Calculating mode value of a vector
### Aliases: fn_mode

### ** Examples

fn_mode(c(1,2,3,1,4,1,7))



cleanEx()
nameEx("fn_target")
### * fn_target

flush(stderr()); flush(stdout())

### Name: fn_target
### Title: Redefines target value
### Aliases: fn_target

### ** Examples

data <- iris
data$Species <- as.character(data$Species)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)

data2 <- fn_target(base = data,target = "Y",event = 1)



cleanEx()
nameEx("gini_table")
### * gini_table

flush(stderr()); flush(stdout())

### Name: gini_table
### Title: Performance measure table with Gini coefficient, KS-statistics
###   and Gini lift curve
### Aliases: gini_table

### ** Examples

data <- iris
data$Species <- as.character(data$Species)
suppressWarnings(RNGversion('3.5.0'))
set.seed(11)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
suppressWarnings(RNGversion('3.5.0'))
set.seed(11)
data$Y_pred <- sample(300:900,size=nrow(data),replace=TRUE)
gini_tab_list <- gini_table(base = data,target = "Y",col_pred = "Y_pred",quantile_pt = 10)
gini_tab_list$prediction
gini_tab_list$gini_tab
gini_tab_list$gini_value
gini_tab_list$gini_plot
gini_tab_list$ks_value
gini_tab_list$breaks



cleanEx()
nameEx("gradient_boosting_parameters")
### * gradient_boosting_parameters

flush(stderr()); flush(stdout())

### Name: gradient_boosting_parameters
### Title: Hyperparameter optimisation or parameter tuning for Gradient
###   Boosting Regression Modelling by grid search
### Aliases: gradient_boosting_parameters

### ** Examples

data <- iris
suppressWarnings(RNGversion('3.5.0'))
set.seed(11)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
gbm_params_list <- gradient_boosting_parameters(base = data,target = "Y",ntree = 2,depth = 2,
                   shrinkage = 0.1,min_obs = 0.1,bag_fraction = 0.7)
gbm_params_list$error_tab_detailed
gbm_params_list$error_tab_summary
gbm_params_list$best_ntree
gbm_params_list$best_depth
gbm_params_list$best_shrinkage
gbm_params_list$best_min_obs
gbm_params_list$best_bag_fraction
gbm_params_list$runtime



cleanEx()
nameEx("iv_filter")
### * iv_filter

flush(stderr()); flush(stdout())

### Name: iv_filter
### Title: Variable reduction based on Information Value filter
### Aliases: iv_filter

### ** Examples

data <- iris
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
x <- c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")
iv_table_list <- iv_table(base = data,target = "Y",num_var_name = x,cat_var_name = "Species")
ivf_list <- iv_filter(base = data,iv_table = iv_table_list$iv_table,threshold = 0.02)
ivf_list$retain_var_tab
ivf_list$retain_var_name
ivf_list$dropped_var_tab
ivf_list$threshold



cleanEx()
nameEx("iv_table")
### * iv_table

flush(stderr()); flush(stdout())

### Name: iv_table
### Title: WOE and IV table for list of numerical and categorical variables
### Aliases: iv_table

### ** Examples

data <- iris
data$Species <- as.character(data$Species)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
x <- c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")
iv_table_list <- iv_table(base = data,target = "Y",num_var_name = x,cat_var_name = "Species")
iv_table_list$num_woe_table
iv_table_list$cat_woe_table
iv_table_list$woe_table
iv_table_list$iv_table



cleanEx()
nameEx("missing_val")
### * missing_val

flush(stderr()); flush(stdout())

### Name: missing_val
### Title: Missing value imputation
### Aliases: missing_val

### ** Examples

data <- iris
data$Species <- as.character(data$Species)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
data[sample(1:nrow(data),size=25),"Sepal.Length"] <- NA
data[sample(1:nrow(data),size=10),"Species"] <- NA

missing_list <- missing_val(base = data,target = "Y")
missing_list$base
missing_list$mapping_table



cleanEx()
nameEx("num_to_cat")
### * num_to_cat

flush(stderr()); flush(stdout())

### Name: num_to_cat
### Title: Binning numerical variables based on cuts from IV table
### Aliases: num_to_cat

### ** Examples

data <- iris
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
x <- c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")
iv_table_list <- iv_table(base = data,target = "Y",num_var_name = x,cat_var_name = "Species")
num_cat <- num_to_cat(base = data,num_woe_table = iv_table_list$num_woe_table)



cleanEx()
nameEx("others_class")
### * others_class

flush(stderr()); flush(stdout())

### Name: others_class
### Title: Clubbing of classes of categorical variable with low population
###   percentage into one class
### Aliases: others_class

### ** Examples

data <- iris[c(1:110),]
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
data$Species <- as.character(data$Species)
data_otherclass <- others_class(base = data,target = "Y",column_name = "Species",threshold = 0.15)



cleanEx()
nameEx("random_forest_parameters")
### * random_forest_parameters

flush(stderr()); flush(stdout())

### Name: random_forest_parameters
### Title: Hyperparameter optimisation or parameter tuning for Random
###   Forest by grid search
### Aliases: random_forest_parameters

### ** Examples

data <- iris
suppressWarnings(RNGversion('3.5.0'))
set.seed(11)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
rf_params_list <- random_forest_parameters(base = data,target = "Y",
                  model_type = "classification",ntree = 2,mtry = 1,nodesize = 3)
rf_params_list$error_tab_detailed
rf_params_list$error_tab_summary
rf_params_list$best_ntree
rf_params_list$best_mtry
rf_params_list$maxnodes
rf_params_list$best_nodesize
rf_params_list$runtime



cleanEx()
nameEx("sampling")
### * sampling

flush(stderr()); flush(stdout())

### Name: sampling
### Title: Random sampling of data into train and test
### Aliases: sampling

### ** Examples

data <- iris
sampling_list <- sampling(base = data,train_perc = 0.7,seed = 1234)
sampling_list$train
sampling_list$test
sampling_list$seed



cleanEx()
nameEx("scalling")
### * scalling

flush(stderr()); flush(stdout())

### Name: scalling
### Title: Converting coefficients of logistic regression into scores for
###   scorecard building
### Aliases: scalling

### ** Examples

data <- iris
suppressWarnings(RNGversion('3.5.0'))
set.seed(11)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
x <- c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")
iv_table_list <- iv_table(base = data,target = "Y",num_var_name = x,cat_var_name = "Species")
num_cat <- num_to_cat(base = data,num_woe_table = iv_table_list$num_woe_table)
log_model <- glm(Y ~ ., data = num_cat, family = "binomial")
scaling_tab <- scalling(base = num_cat,target = "Y",model = log_model)



cleanEx()
nameEx("scoring")
### * scoring

flush(stderr()); flush(stdout())

### Name: scoring
### Title: Scoring a dataset with class based on a scalling logic to arrive
###   at final score
### Aliases: scoring

### ** Examples

data <- iris
suppressWarnings(RNGversion('3.5.0'))
set.seed(11)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
x <- c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")
iv_table_list <- iv_table(base = data,target = "Y",num_var_name = x,cat_var_name = "Species")
num_cat <- num_to_cat(base = data,num_woe_table = iv_table_list$num_woe_table)
log_model <- glm(Y ~ ., data = num_cat, family = "binomial")
scaling_tab <- scalling(base = num_cat,target = "Y",model = log_model)
score_tab <- scoring(base = num_cat,target = "Y",scalling = scaling_tab)



cleanEx()
nameEx("support_vector_parameters")
### * support_vector_parameters

flush(stderr()); flush(stdout())

### Name: support_vector_parameters
### Title: Hyperparameter optimisation or parameter tuning for Suppert
###   Vector Machine by grid search
### Aliases: support_vector_parameters

### ** Examples

data <- iris
suppressWarnings(RNGversion('3.5.0'))
set.seed(11)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
svm_params_list <- support_vector_parameters(base = data,target = "Y",gamma = 0.1,
                   cost = 0.1,kernel = "radial")
svm_params_list$error_tab_detailed
svm_params_list$error_tab_summary
svm_params_list$best_kernel
svm_params_list$best_degree
svm_params_list$best_gamma
svm_params_list$best_cost
svm_params_list$runtime



cleanEx()
nameEx("univariate")
### * univariate

flush(stderr()); flush(stdout())

### Name: univariate
### Title: Univariate analysis of variables
### Aliases: univariate

### ** Examples

data <- iris
data$Species <- as.character(data$Species)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)

univariate_list <- univariate(base = data,target = "Y",threshold = 0.95)
univariate_list$univar_table
univariate_list$num_var_name
univariate_list$char_var_name
univariate_list$sparse_var_name



cleanEx()
nameEx("vif_filter")
### * vif_filter

flush(stderr()); flush(stdout())

### Name: vif_filter
### Title: Removing multicollinearity from a model using vif test
### Aliases: vif_filter

### ** Examples

data <- iris
suppressWarnings(RNGversion('3.5.0'))
set.seed(11)
data$Y <- sample(0:1,size=nrow(data),replace=TRUE)
vif_data_list <- vif_filter(base = data,target = "Y")
vif_data_list$vif_table
vif_data_list$model
vif_data_list$retain_var_list
vif_data_list$dropped_var_list
vif_data_list$threshold



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
