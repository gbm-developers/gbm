pkgname <- "regressoR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('regressoR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("as_string_c")
### * as_string_c

flush(stderr()); flush(stdout())

### Name: as_string_c
### Title: as_string_c
### Aliases: as_string_c

### ** Examples

as_string_c(c("A", "B", "C"))
as_string_c(c(5, 6, 7))
as_string_c(c(5, 6, 7), quote = FALSE)
as_string_c(iris$Species)




cleanEx()
nameEx("boosting_importance_plot")
### * boosting_importance_plot

flush(stderr()); flush(stdout())

### Name: boosting_importance_plot
### Title: boosting_importance_plot
### Aliases: boosting_importance_plot

### ** Examples

## Not run: 
##D library(gbm)
##D library(ggplot2)
##D library(forcats)
##D library(dplyr)
##D 
##D x <- boosting_model('iris', 'Petal.Length', "model_boosting")
##D exe(x)
##D 
##D x <- boosting_importance_plot('model_boosting', 'iris')
##D exe(x)
## End(Not run)



cleanEx()
nameEx("boosting_model")
### * boosting_model

flush(stderr()); flush(stdout())

### Name: boosting_model
### Title: boosting_model
### Aliases: boosting_model

### ** Examples

library(gbm)
library(dplyr)

x <- boosting_model('iris', 'Petal.Length')
exe(x)
print(modelo.boosting)




cleanEx()
nameEx("boosting_prediction")
### * boosting_prediction

flush(stderr()); flush(stdout())

### Name: boosting_prediction
### Title: boosting_prediction
### Aliases: boosting_prediction

### ** Examples

library(gbm)
library(dplyr)
x <- boosting_model('iris', 'Petal.Length', "model_boosting")
exe(x)
print(model_boosting)

x <- boosting_prediction('iris', 'Petal.Length', 'model_boosting', 'my_prediction')
exe(x)
print(my_prediction)




cleanEx()
nameEx("calibrate_boosting")
### * calibrate_boosting

flush(stderr()); flush(stdout())

### Name: calibrate_boosting
### Title: calibrate_boosting
### Aliases: calibrate_boosting

### ** Examples

calibrate_boosting(iris)




cleanEx()
nameEx("categorical_distribution")
### * categorical_distribution

flush(stderr()); flush(stdout())

### Name: categorical_distribution
### Title: categorical_distribution
### Aliases: categorical_distribution

### ** Examples

categorical_distribution(iris$Species)




cleanEx()
nameEx("categorical_summary")
### * categorical_summary

flush(stderr()); flush(stdout())

### Name: categorical_summary
### Title: categorical_summary
### Aliases: categorical_summary

### ** Examples

if(interactive()) {
  library(shiny)
  library(DT)
  shinyApp(ui = fluidPage(fluidRow(uiOutput("resumen"))),
           server = function(input, output) {
                         output$resumen = renderUI(categorical_summary(iris, "Species"))
           })
}



cleanEx()
nameEx("clean_report")
### * clean_report

flush(stderr()); flush(stdout())

### Name: clean_report
### Title: clean_report
### Aliases: clean_report

### ** Examples

new_report(iris, 'iris')
get_report()
clean_report()
get_report()




cleanEx()
nameEx("code_NA")
### * code_NA

flush(stderr()); flush(stdout())

### Name: code_NA
### Title: code_NA
### Aliases: code_NA

### ** Examples

iris2 <- iris
x <- code_NA(TRUE, 'iris2')
exe(x)
x <- code_NA(FALSE, 'iris2')
exe(x)




cleanEx()
nameEx("code_deactivate")
### * code_deactivate

flush(stderr()); flush(stdout())

### Name: code_deactivate
### Title: code_deactivate
### Aliases: code_deactivate

### ** Examples

iris2 <- iris
x <- code_deactivate('Species', 'iris2')
exe(x)
head(iris2)




cleanEx()
nameEx("code_load")
### * code_load

flush(stderr()); flush(stdout())

### Name: code_load
### Title: code_load
### Aliases: code_load

### ** Examples

code_load(TRUE, "MY/PATH/FILE.csv")




cleanEx()
nameEx("code_summary")
### * code_summary

flush(stderr()); flush(stdout())

### Name: code_summary
### Title: code_summary
### Aliases: code_summary

### ** Examples

x <- code_summary('iris')
exe(x)




cleanEx()
nameEx("code_transf")
### * code_transf

flush(stderr()); flush(stdout())

### Name: code_transf
### Title: code_transf
### Aliases: code_transf

### ** Examples

iris2 <- iris
x <-code_transf('Species', 'disyuntivo', 'iris', 'iris2')
exe(x)
head(iris2)




cleanEx()
nameEx("coef_lambda")
### * coef_lambda

flush(stderr()); flush(stdout())

### Name: coef_lambda
### Title: coef_lambda
### Aliases: coef_lambda

### ** Examples

library(glmnet)
x <- rlr_model('iris', 'Petal.Length')
exe(x)

x <- coef_lambda('iris','Petal.Length', 'modelo.rlr')
exe(x)




cleanEx()
nameEx("colnames_empty")
### * colnames_empty

flush(stderr()); flush(stdout())

### Name: colnames_empty
### Title: colnames_empty
### Aliases: colnames_empty

### ** Examples

colnames_empty(iris)
colnames_empty(NULL)




cleanEx()
nameEx("combine_names")
### * combine_names

flush(stderr()); flush(stdout())

### Name: combine_names
### Title: combine_names
### Aliases: combine_names

### ** Examples

x = c("A", "B", "C")
y = c("1", "2", "3")
combine_names(x, y)




cleanEx()
nameEx("comparative_table")
### * comparative_table

flush(stderr()); flush(stdout())

### Name: comparative_table
### Title: comparative_table
### Aliases: comparative_table

### ** Examples

models <- list('knnl-mode1' = list(0.11,0.22,0.33,0.44),
               'dtl-mode2'  = list(0.12,0.23,0.34,0.45),
               'rfl-mode1'  = list(0.51,0.42,0.13,0.24))
sel <- c("K Vecinos Más Cercanos-mode1", "Bosques Aleatorios-mode1")
comparative_table(sel, models)




cleanEx()
nameEx("cor_model")
### * cor_model

flush(stderr()); flush(stdout())

### Name: cor_model
### Title: cor_model
### Aliases: cor_model

### ** Examples

x <- cor_model('iris')
exe(x)
correlacion




cleanEx()
nameEx("correlations_plot")
### * correlations_plot

flush(stderr()); flush(stdout())

### Name: correlations_plot
### Title: correlations_plot
### Aliases: correlations_plot

### ** Examples

x <- cor_model('iris')
exe(x)
print(correlacion)

x <- correlations_plot()
exe(x)




cleanEx()
nameEx("def_code_cat")
### * def_code_cat

flush(stderr()); flush(stdout())

### Name: def_code_cat
### Title: def_code_cat
### Aliases: def_code_cat

### ** Examples

x <- def_code_cat('iris', 'Species')
exe(x)




cleanEx()
nameEx("def_code_num")
### * def_code_num

flush(stderr()); flush(stdout())

### Name: def_code_num
### Title: def_code_num
### Aliases: def_code_num

### ** Examples

x <- def_code_num('iris', 'Petal.Length')
exe(x)




cleanEx()
nameEx("default_calc_normal")
### * default_calc_normal

flush(stderr()); flush(stdout())

### Name: default_calc_normal
### Title: default_calc_normal
### Aliases: default_calc_normal

### ** Examples

x <- default_calc_normal('iris')
exe(x)




cleanEx()
nameEx("default_disp")
### * default_disp

flush(stderr()); flush(stdout())

### Name: default_disp
### Title: default_disp
### Aliases: default_disp

### ** Examples

library(scatterplot3d)

x <- default_disp('iris', c('Sepal.Length', 'Sepal.Width'))
exe(x)

x <- default_disp('iris', c('Sepal.Length', 'Sepal.Width', 'Petal.Length'))
exe(x)




cleanEx()
nameEx("disjunctive_data")
### * disjunctive_data

flush(stderr()); flush(stdout())

### Name: disjunctive_data
### Title: disjunctive_data
### Aliases: disjunctive_data

### ** Examples

disjunctive_data(iris, "Species")




cleanEx()
nameEx("disp_models")
### * disp_models

flush(stderr()); flush(stdout())

### Name: disp_models
### Title: disp_models
### Aliases: disp_models

### ** Examples

disp_models("prediction.knn", "KNN", "Species")




cleanEx()
nameEx("dt_model")
### * dt_model

flush(stderr()); flush(stdout())

### Name: dt_model
### Title: dt_model
### Aliases: dt_model

### ** Examples

library(rpart)

x <- dt_model('iris', 'Petal.Length')
exe(x)
print(modelo.dt)




cleanEx()
nameEx("dt_plot")
### * dt_plot

flush(stderr()); flush(stdout())

### Name: dt_plot
### Title: dt_plot
### Aliases: dt_plot

### ** Examples

## Not run: 
##D library(rpart)
##D 
##D x <- dt_model('iris', 'Petal.Length', model.var = 'model_dt')
##D exe(x)
##D print(model_dt)
##D 
##D x <- dt_plot('model_dt')
##D exe(x)
## End(Not run)



cleanEx()
nameEx("dt_prediction")
### * dt_prediction

flush(stderr()); flush(stdout())

### Name: dt_prediction
### Title: dt_prediction
### Aliases: dt_prediction

### ** Examples

library(rpart)

x <- dt_model('iris', 'Petal.Length', model.var = 'model_dt')
exe(x)
print(model_dt)

x <- dt_prediction('iris', 'model_dt', 'my_prediction')
exe(x)
print(my_prediction)




cleanEx()
nameEx("error_plot")
### * error_plot

flush(stderr()); flush(stdout())

### Name: error_plot
### Title: error_plot
### Aliases: error_plot

### ** Examples

error_plot("My Warning")




cleanEx()
nameEx("error_variables")
### * error_variables

flush(stderr()); flush(stdout())

### Name: error_variables
### Title: error_variables
### Aliases: error_variables

### ** Examples

error_variables(TRUE)
error_variables(FALSE)




cleanEx()
nameEx("exe")
### * exe

flush(stderr()); flush(stdout())

### Name: exe
### Title: exe
### Aliases: exe

### ** Examples

exe("5+5")
exe("5","+","5")
exe("plot(iris$Species)")




cleanEx()
nameEx("extract_code")
### * extract_code

flush(stderr()); flush(stdout())

### Name: extract_code
### Title: extract_code
### Aliases: extract_code

### ** Examples

extract_code("cat")
extract_code("plot")

parse(text = extract_code("plot"))




cleanEx()
nameEx("fisher_calc")
### * fisher_calc

flush(stderr()); flush(stdout())

### Name: fisher_calc
### Title: fisher_calc
### Aliases: fisher_calc

### ** Examples

fisher_calc(iris$Petal.Length)




cleanEx()
nameEx("general_indices")
### * general_indices

flush(stderr()); flush(stdout())

### Name: general_indices
### Title: general_indices
### Aliases: general_indices

### ** Examples

real <- rnorm(45)
prediction <- rnorm(45)
model <- "KNN"
general_indices(real, prediction)




cleanEx()
nameEx("get_env_report")
### * get_env_report

flush(stderr()); flush(stdout())

### Name: get_env_report
### Title: get_env_report
### Aliases: get_env_report

### ** Examples

e <- get_env_report()
e$codigo.reporte




cleanEx()
nameEx("get_report")
### * get_report

flush(stderr()); flush(stdout())

### Name: get_report
### Title: get_report
### Aliases: get_report

### ** Examples

get_report()




cleanEx()
nameEx("gg_color_hue")
### * gg_color_hue

flush(stderr()); flush(stdout())

### Name: gg_color_hue
### Title: gg_color_hue
### Aliases: gg_color_hue

### ** Examples

col <- gg_color_hue(3)
plot(iris$Species, col = col)




cleanEx()
nameEx("importance_plot_rf")
### * importance_plot_rf

flush(stderr()); flush(stdout())

### Name: importance_plot_rf
### Title: importance_plot_rf
### Aliases: importance_plot_rf

### ** Examples

library(randomForest)
x <- rf_model('iris', 'Petal.Length')
exe(x)
importance_plot_rf(modelo.rf, translate('impVarA'), translate('impVarRSS'))




cleanEx()
nameEx("init_regressor")
### * init_regressor

flush(stderr()); flush(stdout())

### Name: init_regressor
### Title: This function will start regressoR
### Aliases: init_regressor
### Keywords: regressoR

### ** Examples

 if(interactive()){
   init_regressor()
 }



cleanEx()
nameEx("insert_report")
### * insert_report

flush(stderr()); flush(stdout())

### Name: insert_report
### Title: insert_report
### Aliases: insert_report

### ** Examples

new_report(iris, "iris")
insert_report("1_part", 'Title 1', 'head(iris)\n', 'summary(iris)')
get_report()
clean_report()




cleanEx()
nameEx("kkn_model")
### * kkn_model

flush(stderr()); flush(stdout())

### Name: kkn_model
### Title: kkn_model
### Aliases: kkn_model

### ** Examples

library(kknn)
x <- kkn_model('iris', 'Petal.Length')
exe(x)
print(modelo.knn)




cleanEx()
nameEx("kkn_prediction")
### * kkn_prediction

flush(stderr()); flush(stdout())

### Name: kkn_prediction
### Title: kkn_prediction
### Aliases: kkn_prediction

### ** Examples

library(kknn)
library(dplyr)

x <- kkn_model('iris', 'Petal.Length', model.var = 'model_knn')
exe(x)
print(model_knn)

x <- kkn_prediction('iris', 'Petal.Length', 'model_knn', 'my_prediction')
exe(x)
print(my_prediction)




cleanEx()
nameEx("models_mode")
### * models_mode

flush(stderr()); flush(stdout())

### Name: models_mode
### Title: models_mode
### Aliases: models_mode

### ** Examples

x <- list('knnl-mode1' = 1, 'knnl-mode2' = 2, 'knnl-mode2' = 5)
models_mode(x)




cleanEx()
nameEx("new_col")
### * new_col

flush(stderr()); flush(stdout())

### Name: new_col
### Title: new_col
### Aliases: new_col

### ** Examples

new_col(iris)
new_col(iris, "var1", c(1,2,3))




cleanEx()
nameEx("new_report")
### * new_report

flush(stderr()); flush(stdout())

### Name: new_report
### Title: new_report
### Aliases: new_report

### ** Examples

new_report(iris, "iris")
get_report()
clean_report()




cleanEx()
nameEx("new_section_report")
### * new_section_report

flush(stderr()); flush(stdout())

### Name: new_section_report
### Title: new_section_report
### Aliases: new_section_report

### ** Examples

new_report(iris, 'iris')
insert_report('1_part', 'Title 1', 'head(iris)\n', 'summary(iris)')
get_report()

remove_report_elem('1_part')
get_report()

new_section_report()
insert_report('1_part', 'Title 1', 'head(iris)\n', 'summary(iris)')
get_report()

new_section_report()
insert_report('1_part', 'Title 1', 'head(iris)\n', 'summary(iris)')
get_report()

remove_report_elem('1_part')
get_report()

clean_report()




cleanEx()
nameEx("nn_model")
### * nn_model

flush(stderr()); flush(stdout())

### Name: nn_model
### Title: nn_model
### Aliases: nn_model

### ** Examples

## Not run: 
##D library(neuralnet)
##D library(dummies)
##D 
##D x <- nn_model('iris', 'Petal.Length','modelo.nn', 'mean.nn', 'sd.nn', 0.05, 2000, 3, 30, 50, 80)
##D exe(x)
##D 
##D print(modelo.nn)
##D print(mean.nn)
##D print(sd.nn)
## End(Not run)



cleanEx()
nameEx("nn_plot")
### * nn_plot

flush(stderr()); flush(stdout())

### Name: nn_plot
### Title: nn_plot
### Aliases: nn_plot

### ** Examples

## Not run: 
##D library(neuralnet)
##D library(dummies)
##D library(dplyr)
##D 
##D x <- nn_model('iris', 'Petal.Length','modelo.nn', 'mean.nn', 'sd.nn', 0.05, 2000, 3, 10, 10, 10)
##D exe(x)
##D 
##D x <- nn_plot('modelo.nn')
##D exe(x)
## End(Not run)



cleanEx()
nameEx("nn_prediction")
### * nn_prediction

flush(stderr()); flush(stdout())

### Name: nn_prediction
### Title: nn_prediction
### Aliases: nn_prediction

### ** Examples

## Not run: 
##D library(neuralnet)
##D library(dummies)
##D library(dplyr)
##D 
##D x <- nn_model('iris', 'Petal.Length','modelo.nn', 'mean.nn', 'sd.nn', 0.05, 2000, 3, 30, 50, 30)
##D exe(x)
##D 
##D x <- nn_prediction('iris', 'Petal.Length')
##D exe(x)
##D print(prediccion.nn)
## End(Not run)



cleanEx()
nameEx("normal_default")
### * normal_default

flush(stderr()); flush(stdout())

### Name: normal_default
### Title: normal_default
### Aliases: normal_default

### ** Examples

x <- normal_default('iris', 'Sepal.Length')
exe(x)




cleanEx()
nameEx("numerical_distribution")
### * numerical_distribution

flush(stderr()); flush(stdout())

### Name: numerical_distribution
### Title: numerical_distribution
### Aliases: numerical_distribution

### ** Examples

numerical_distribution(iris[,'Sepal.Length'], 'Sepal.Length', 'red')




cleanEx()
nameEx("numerical_summary")
### * numerical_summary

flush(stderr()); flush(stdout())

### Name: numerical_summary
### Title: numerical_summary
### Aliases: numerical_summary

### ** Examples

if(interactive()) {
  library(shiny)
  library(DT)
  shinyApp(ui = fluidPage(fluidRow(uiOutput("resumen"))),
           server = function(input, output) {
                         output$resumen = renderUI(numerical_summary(iris, "Sepal.Width"))
           })
}



cleanEx()
nameEx("options_regressor")
### * options_regressor

flush(stderr()); flush(stdout())

### Name: options_regressor
### Title: options_regressor
### Aliases: options_regressor

### ** Examples

options_regressor("language")
options_regressor(language = "en")
options_regressor("language")




cleanEx()
nameEx("pairs_power")
### * pairs_power

flush(stderr()); flush(stdout())

### Name: pairs_power
### Title: pairs_power
### Aliases: pairs_power

### ** Examples

## Not run: 
##D library(psych)
##D x <- pairs_power('iris')
##D exe(x)
## End(Not run)




cleanEx()
nameEx("partition_code")
### * partition_code

flush(stderr()); flush(stdout())

### Name: partition_code
### Title: partition_code
### Aliases: partition_code

### ** Examples

x <- partition_code('iris', 75, 'Species', 555, TRUE)
exe(x)
head(datos.aprendizaje)
head(datos.prueba)




cleanEx()
nameEx("plot_RMSE")
### * plot_RMSE

flush(stderr()); flush(stdout())

### Name: plot_RMSE
### Title: plot_RMSE
### Aliases: plot_RMSE

### ** Examples

library(pls)

x <- rd_model('iris', 'Petal.Length')
exe(x)

plot_RMSE(modelo.rd)




cleanEx()
nameEx("plot_coef_lambda")
### * plot_coef_lambda

flush(stderr()); flush(stdout())

### Name: plot_coef_lambda
### Title: plot_coef_lambda
### Aliases: plot_coef_lambda

### ** Examples

library(glmnet)
x <- rlr_model('iris', 'Petal.Length')
exe(x)

x <- plot_coef_lambda('modelo.rlr')
exe(x)




cleanEx()
nameEx("plot_pred_rd")
### * plot_pred_rd

flush(stderr()); flush(stdout())

### Name: plot_pred_rd
### Title: plot_pred_rd
### Aliases: plot_pred_rd

### ** Examples

library(pls)

x <- rd_model('iris', 'Petal.Length')
exe(x)

plot_pred_rd(modelo.rd)




cleanEx()
nameEx("plot_real_prediction")
### * plot_real_prediction

flush(stderr()); flush(stdout())

### Name: plot_real_prediction
### Title: plot_real_prediction
### Aliases: plot_real_prediction

### ** Examples

real <- rnorm(45)
prediction <- rnorm(45)
model <- "KNN"
plot_real_prediction(real, prediction, model)




cleanEx()
nameEx("plot_var_pred_rd")
### * plot_var_pred_rd

flush(stderr()); flush(stdout())

### Name: plot_var_pred_rd
### Title: plot_var_pred_rd
### Aliases: plot_var_pred_rd

### ** Examples

library(pls)

x <- rd_model('iris', 'Petal.Length')
exe(x)

plot_var_pred_rd(modelo.rd)




cleanEx()
nameEx("rd_model")
### * rd_model

flush(stderr()); flush(stdout())

### Name: rd_model
### Title: rd_model
### Aliases: rd_model

### ** Examples

library(pls)

x <- rd_model('iris', 'Petal.Length')
exe(x)
print(modelo.rd)




cleanEx()
nameEx("rd_prediction")
### * rd_prediction

flush(stderr()); flush(stdout())

### Name: rd_prediction
### Title: rd_prediction
### Aliases: rd_prediction

### ** Examples

library(pls)

x <- rd_model('iris', 'Petal.Length')
exe(x)
print(modelo.rd)

x <- rd_prediction('iris', 'modelo.rd', 'my_prediction')
exe(x)
print(my_prediction)




cleanEx()
nameEx("rd_type")
### * rd_type

flush(stderr()); flush(stdout())

### Name: rd_type
### Title: rd_type
### Aliases: rd_type

### ** Examples

rd_type(1)
rd_type(0)




cleanEx()
nameEx("remove_report_elem")
### * remove_report_elem

flush(stderr()); flush(stdout())

### Name: remove_report_elem
### Title: remove_report_elem
### Aliases: remove_report_elem

### ** Examples

new_report(iris, 'iris')
insert_report('1_part', 'Title 1', 'head(iris)\n', 'summary(iris)')
get_report()
remove_report_elem('1_part')
get_report()
clean_report()




cleanEx()
nameEx("render_index_table")
### * render_index_table

flush(stderr()); flush(stdout())

### Name: render_index_table
### Title: render_index_table
### Aliases: render_index_table

### ** Examples

if(interactive()) {
   library(shiny)
   shinyApp(
      ui = fluidPage(fluidRow(column(12, tableOutput('tbl')))),
      server = function(input, output) {
         output$tbl = render_index_table(iris)
      }
   )
}




cleanEx()
nameEx("render_table_data")
### * render_table_data

flush(stderr()); flush(stdout())

### Name: render_table_data
### Title: render_table_data
### Aliases: render_table_data

### ** Examples

if(interactive()) {
   library(shiny)
   library(DT)
   shinyApp(
      ui = fluidPage(fluidRow(column(12, DTOutput('tbl')))),
      server = function(input, output) {
         output$tbl = render_table_data(iris)
      }
   )
}




cleanEx()
nameEx("rf_model")
### * rf_model

flush(stderr()); flush(stdout())

### Name: rf_model
### Title: rf_model
### Aliases: rf_model

### ** Examples

library(randomForest)
x <- rf_model('iris', 'Petal.Length')
exe(x)
print(modelo.rf)




cleanEx()
nameEx("rf_prediction")
### * rf_prediction

flush(stderr()); flush(stdout())

### Name: rf_prediction
### Title: rf_prediction
### Aliases: rf_prediction

### ** Examples

library(randomForest)
library(dplyr)

x <- rf_model('iris', 'Petal.Length', model.var = 'model_rf')
exe(x)
print(model_rf)

x <- rf_prediction('iris', 'Petal.Length', 'model_rf', 'my_prediction')
exe(x)
print(my_prediction)




cleanEx()
nameEx("rl_coeff")
### * rl_coeff

flush(stderr()); flush(stdout())

### Name: rl_coeff
### Title: rl_coeff
### Aliases: rl_coeff

### ** Examples

x <- rl_model('iris', 'Petal.Length', 'model_rl')
exe(x)
print(model_rl)

x <- rl_coeff('model_rl')
exe(x)

print(df.rl)
print(r2)




cleanEx()
nameEx("rl_model")
### * rl_model

flush(stderr()); flush(stdout())

### Name: rl_model
### Title: rl_model
### Aliases: rl_model

### ** Examples

x <- rl_model('iris', 'Petal.Length')
exe(x)
print(modelo.rl)




cleanEx()
nameEx("rl_prediction")
### * rl_prediction

flush(stderr()); flush(stdout())

### Name: rl_prediction
### Title: rl_prediction
### Aliases: rl_prediction

### ** Examples

x <- rl_model('iris', 'Petal.Length', 'model_rl')
exe(x)
print(model_rl)

x <- rl_prediction('iris', 'model_rl', 'my_prediction')
exe(x)
print(my_prediction)




cleanEx()
nameEx("rlr_model")
### * rlr_model

flush(stderr()); flush(stdout())

### Name: rlr_model
### Title: rlr_model
### Aliases: rlr_model

### ** Examples

library(glmnet)
x <- rlr_model('iris', 'Petal.Length')
exe(x)
print(modelo.rlr)




cleanEx()
nameEx("rlr_prediction")
### * rlr_prediction

flush(stderr()); flush(stdout())

### Name: rlr_prediction
### Title: rlr_prediction
### Aliases: rlr_prediction

### ** Examples

library(glmnet)
x <- rlr_model('iris', 'Petal.Length')
exe(x)
print(modelo.rlr)

x <- rlr_prediction('iris', 'iris', 'Petal.Length', pred.var = 'my_prediction')
exe(x)
print(my_prediction)




cleanEx()
nameEx("rlr_type")
### * rlr_type

flush(stderr()); flush(stdout())

### Name: rlr_type
### Title: rlr_type
### Aliases: rlr_type

### ** Examples

rlr_type(1)
rlr_type(0)




cleanEx()
nameEx("summary_indices")
### * summary_indices

flush(stderr()); flush(stdout())

### Name: summary_indices
### Title: summary_indices
### Aliases: summary_indices

### ** Examples

summary_indices(iris$Sepal.Length)




cleanEx()
nameEx("svm_model")
### * svm_model

flush(stderr()); flush(stdout())

### Name: svm_model
### Title: svm_model
### Aliases: svm_model

### ** Examples

library(e1071)
x <- svm_model('iris', 'Petal.Length')
exe(x)
print(modelo.svm)




cleanEx()
nameEx("svm_prediction")
### * svm_prediction

flush(stderr()); flush(stdout())

### Name: svm_prediction
### Title: svm_prediction
### Aliases: svm_prediction

### ** Examples

library(e1071)
library(dplyr)

x <- svm_model('iris', 'Petal.Length', model.var = 'model_svm')
exe(x)
print(model_svm)

x <- svm_prediction('iris', 'Petal.Length', 'model_svm', 'my_prediction')
exe(x)
print(my_prediction)




cleanEx()
nameEx("tb_predic")
### * tb_predic

flush(stderr()); flush(stdout())

### Name: tb_predic
### Title: tb_predic
### Aliases: tb_predic

### ** Examples

if(interactive()) {
  library(shiny)
  library(DT)
  shinyApp( 
    ui = fluidPage(fluidRow(column(12, DTOutput('tbl')))),
   server = function(input, output) {
     real <- iris[,'Petal.Width',drop = F]
     pred <- sample(iris$Petal.Width, nrow(iris), replace =  T)
     output$tbl = DT::renderDT(tb_predic(real, pred))
   })
}




cleanEx()
nameEx("translate")
### * translate

flush(stderr()); flush(stdout())

### Name: translate
### Title: translate
### Aliases: translate

### ** Examples

translate("knnl")
translate("knnl", "en")




cleanEx()
nameEx("validate_pn_data")
### * validate_pn_data

flush(stderr()); flush(stdout())

### Name: validate_pn_data
### Title: validate_pn_data
### Aliases: validate_pn_data

### ** Examples

## Not run: 
##D validate_pn_data(iris, cars)
##D validate_pn_data(iris, iris)
##D x <- iris
##D x$Species <- as.numeric(x$Species)
##D validate_pn_data(iris, x)
## End(Not run)



cleanEx()
nameEx("var_categorical")
### * var_categorical

flush(stderr()); flush(stdout())

### Name: var_categorical
### Title: var_categorical
### Aliases: var_categorical

### ** Examples

var_categorical(iris)




cleanEx()
nameEx("var_numerical")
### * var_numerical

flush(stderr()); flush(stdout())

### Name: var_numerical
### Title: var_numerical
### Aliases: var_numerical

### ** Examples

var_numerical(iris)




cleanEx()
nameEx("word_report")
### * word_report

flush(stderr()); flush(stdout())

### Name: word_report
### Title: word_report
### Aliases: word_report

### ** Examples

new_report(iris, 'iris')

new_section_report()
insert_report('1_part', 'Title 1', 'head(iris)\n', 'summary(iris)')

new_section_report()
insert_report('1_part', 'Title 1', 'head(iris)\n', 'summary(iris)')

word_report(order_by_regressor = FALSE)




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
