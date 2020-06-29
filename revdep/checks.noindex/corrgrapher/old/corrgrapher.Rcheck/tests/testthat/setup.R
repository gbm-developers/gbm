library(ranger)
library(DALEX)

set.seed(2020)

data(dragons, package='DALEX')
model <- ranger::ranger(colour ~ ., data = dragons, num.trees = 100, probability = TRUE)
model_exp <- DALEX::explain(model, data = dragons[,-5], y = dragons$colour)
model_fi <- ingredients::feature_importance(model_exp, loss_function = loss_accuracy, type = 'raw')
model_pd <- ingredients::partial_dependence(model_exp, N=100, grid_points = 81)

cgr_exp <- corrgrapher(model_exp, 
                       feature_importance = model_fi,
                       partial_dependency = list(numerical = model_pd))

data("titanic_imputed", package='DALEX')
tit_model <- ranger(survived ~ ., data = titanic_imputed, num.trees = 100)
tit_model_exp <- DALEX::explain(tit_model, data = titanic_imputed[,-8], y = titanic_imputed[,8])
tit_model_fi <- ingredients::feature_importance(tit_model_exp)
tit_model_pd <- list(numerical = ingredients::partial_dependence(tit_model_exp, N=100, grid_points = 81),
                     categorical = ingredients::partial_dependence(tit_model_exp, N=100, grid_points = 81, variable_type = 'categorical'))

tit_cgr_exp <- corrgrapher(tit_model_exp, 
                       feature_importance = tit_model_fi,
                       partial_dependency = tit_model_pd)

data("freeny")
simple_model <- lm(y ~ ., data = freeny)
simple_model_exp <- DALEX::explain(simple_model, data = freeny[,-1], y = freeny$y)

data("Seatbelts")
df <- as.data.frame(Seatbelts)[,-8]
cgr_df <- corrgrapher(df)
cgr_df_mixed <- corrgrapher(titanic_imputed[,-8])

print('helper ended!')
