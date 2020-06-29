pkgname <- "MachineShop"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('MachineShop')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("AdaBagModel")
### * AdaBagModel

flush(stderr()); flush(stdout())

### Name: AdaBagModel
### Title: Bagging with Classification Trees
### Aliases: AdaBagModel

### ** Examples

fit(Species ~ ., data = iris, model = AdaBagModel(mfinal = 5))




cleanEx()
nameEx("AdaBoostModel")
### * AdaBoostModel

flush(stderr()); flush(stdout())

### Name: AdaBoostModel
### Title: Boosting with Classification Trees
### Aliases: AdaBoostModel

### ** Examples

fit(Species ~ ., data = iris, model = AdaBoostModel(mfinal = 5))




cleanEx()
nameEx("BARTMachineModel")
### * BARTMachineModel

flush(stderr()); flush(stdout())

### Name: BARTMachineModel
### Title: Bayesian Additive Regression Trees Model
### Aliases: BARTMachineModel

### ** Examples





cleanEx()
nameEx("BARTModel")
### * BARTModel

flush(stderr()); flush(stdout())

### Name: BARTModel
### Title: Bayesian Additive Regression Trees Model
### Aliases: BARTModel

### ** Examples





cleanEx()
nameEx("BlackBoostModel")
### * BlackBoostModel

flush(stderr()); flush(stdout())

### Name: BlackBoostModel
### Title: Gradient Boosting with Regression Trees
### Aliases: BlackBoostModel

### ** Examples

library(MASS)

fit(type ~ ., data = Pima.tr, model = BlackBoostModel)




cleanEx()
nameEx("C50Model")
### * C50Model

flush(stderr()); flush(stdout())

### Name: C50Model
### Title: C5.0 Decision Trees and Rule-Based Model
### Aliases: C50Model

### ** Examples

model_fit <- fit(Species ~ ., data = iris, model = C50Model)
varimp(model_fit, metric = "splits", scale = FALSE)




cleanEx()
nameEx("CForestModel")
### * CForestModel

flush(stderr()); flush(stdout())

### Name: CForestModel
### Title: Conditional Random Forest Model
### Aliases: CForestModel

### ** Examples

fit(sale_amount ~ ., data = ICHomes, model = CForestModel)




cleanEx()
nameEx("CoxModel")
### * CoxModel

flush(stderr()); flush(stdout())

### Name: CoxModel
### Title: Proportional Hazards Regression Model
### Aliases: CoxModel CoxStepAICModel

### ** Examples

library(survival)

fit(Surv(time, status) ~ ., data = veteran, model = CoxModel)




cleanEx()
nameEx("DiscreteVariate")
### * DiscreteVariate

flush(stderr()); flush(stdout())

### Name: DiscreteVariate
### Title: Discrete Variate Constructors
### Aliases: DiscreteVariate BinomialVariate NegBinomialVariate
###   PoissonVariate

### ** Examples

BinomialVariate(rbinom(25, 10, 0.5), size = 10)
PoissonVariate(rpois(25, 10))




cleanEx()
nameEx("EarthModel")
### * EarthModel

flush(stderr()); flush(stdout())

### Name: EarthModel
### Title: Multivariate Adaptive Regression Splines Model
### Aliases: EarthModel

### ** Examples

model_fit <- fit(Species ~ ., data = iris, model = EarthModel)
varimp(model_fit, metric = "nsubsets", scale = FALSE)




cleanEx()
nameEx("FDAModel")
### * FDAModel

flush(stderr()); flush(stdout())

### Name: FDAModel
### Title: Flexible and Penalized Discriminant Analysis Models
### Aliases: FDAModel PDAModel

### ** Examples

fit(Species ~ ., data = iris, model = FDAModel)

fit(Species ~ ., data = iris, model = PDAModel)




cleanEx()
nameEx("GAMBoostModel")
### * GAMBoostModel

flush(stderr()); flush(stdout())

### Name: GAMBoostModel
### Title: Gradient Boosting with Additive Models
### Aliases: GAMBoostModel

### ** Examples

library(MASS)

fit(type ~ ., data = Pima.tr, model = GAMBoostModel)




cleanEx()
nameEx("GBMModel")
### * GBMModel

flush(stderr()); flush(stdout())

### Name: GBMModel
### Title: Generalized Boosted Regression Model
### Aliases: GBMModel

### ** Examples

fit(Species ~ ., data = iris, model = GBMModel)




cleanEx()
nameEx("GLMBoostModel")
### * GLMBoostModel

flush(stderr()); flush(stdout())

### Name: GLMBoostModel
### Title: Gradient Boosting with Linear Models
### Aliases: GLMBoostModel

### ** Examples

library(MASS)

fit(type ~ ., data = Pima.tr, model = GLMBoostModel)




cleanEx()
nameEx("GLMModel")
### * GLMModel

flush(stderr()); flush(stdout())

### Name: GLMModel
### Title: Generalized Linear Model
### Aliases: GLMModel GLMStepAICModel

### ** Examples

fit(sale_amount ~ ., data = ICHomes, model = GLMModel)




cleanEx()
nameEx("GLMNetModel")
### * GLMNetModel

flush(stderr()); flush(stdout())

### Name: GLMNetModel
### Title: GLM Lasso or Elasticnet Model
### Aliases: GLMNetModel

### ** Examples

fit(sale_amount ~ ., data = ICHomes, model = GLMNetModel(lambda = 0.01))




cleanEx()
nameEx("Grid")
### * Grid

flush(stderr()); flush(stdout())

### Name: Grid
### Title: Tuning Grid Control
### Aliases: Grid

### ** Examples

TunedModel(GBMModel, grid = Grid(10, random = 5))




cleanEx()
nameEx("KNNModel")
### * KNNModel

flush(stderr()); flush(stdout())

### Name: KNNModel
### Title: Weighted k-Nearest Neighbor Model
### Aliases: KNNModel

### ** Examples

fit(Species ~ ., data = iris, model = KNNModel)




cleanEx()
nameEx("LARSModel")
### * LARSModel

flush(stderr()); flush(stdout())

### Name: LARSModel
### Title: Least Angle Regression, Lasso and Infinitesimal Forward
###   Stagewise Models
### Aliases: LARSModel

### ** Examples

fit(sale_amount ~ ., data = ICHomes, model = LARSModel)




cleanEx()
nameEx("LDAModel")
### * LDAModel

flush(stderr()); flush(stdout())

### Name: LDAModel
### Title: Linear Discriminant Analysis Model
### Aliases: LDAModel

### ** Examples

fit(Species ~ ., data = iris, model = LDAModel)




cleanEx()
nameEx("LMModel")
### * LMModel

flush(stderr()); flush(stdout())

### Name: LMModel
### Title: Linear Models
### Aliases: LMModel

### ** Examples

fit(sale_amount ~ ., data = ICHomes, model = LMModel)




cleanEx()
nameEx("MDAModel")
### * MDAModel

flush(stderr()); flush(stdout())

### Name: MDAModel
### Title: Mixture Discriminant Analysis Model
### Aliases: MDAModel

### ** Examples

fit(Species ~ ., data = iris, model = MDAModel)




cleanEx()
nameEx("MLControl")
### * MLControl

flush(stderr()); flush(stdout())

### Name: MLControl
### Title: Resampling Controls
### Aliases: MLControl controls BootControl BootOptimismControl CVControl
###   CVOptimismControl OOBControl SplitControl TrainControl

### ** Examples

## Bootstrapping with 100 samples
BootControl(samples = 100)

## Optimism-corrected bootstrapping with 100 samples
BootOptimismControl(samples = 100)

## Cross-validation with 5 repeats of 10 folds
CVControl(folds = 10, repeats = 5)

## Optimism-corrected cross-validation with 5 repeats of 10 folds
CVOptimismControl(folds = 10, repeats = 5)

## Out-of-bootstrap validation with 100 samples
OOBControl(samples = 100)

## Split sample validation with 2/3 training and 1/3 testing
SplitControl(prop = 2/3)

## Training set evaluation
TrainControl()




cleanEx()
nameEx("MLMetric")
### * MLMetric

flush(stderr()); flush(stdout())

### Name: MLMetric
### Title: MLMetric Class Constructor
### Aliases: MLMetric MLMetric<-

### ** Examples

f2_score <- function(observed, predicted, ...) {
  f_score(observed, predicted, beta = 2, ...)
}

MLMetric(f2_score) <- list(name = "f2_score",
                           label = "F Score (beta = 2)",
                           maximize = TRUE)




cleanEx()
nameEx("MLModel")
### * MLModel

flush(stderr()); flush(stdout())

### Name: MLModel
### Title: MLModel Class Constructor
### Aliases: MLModel

### ** Examples

## Logistic regression model
LogisticModel <- MLModel(
  name = "LogisticModel",
  response_types = "binary",
  fit = function(formula, data, weights, ...) {
    glm(formula, data = data, weights = weights, family = binomial, ...)
  },
  predict = function(object, newdata, ...) {
    predict(object, newdata = newdata, type = "response")
  },
  varimp = function(object, ...) {
    pchisq(coef(object)^2 / diag(vcov(object)), 1)
  }
)

library(MASS)
res <- resample(type ~ ., data = Pima.tr, model = LogisticModel)
summary(res)




cleanEx()
nameEx("ModelFrame-methods")
### * ModelFrame-methods

flush(stderr()); flush(stdout())

### Name: ModelFrame
### Title: ModelFrame Class
### Aliases: ModelFrame ModelFrame.formula ModelFrame.matrix

### ** Examples

mf <- ModelFrame(ncases / (ncases + ncontrols) ~ agegp + tobgp + alcgp,
                 data = esoph, weights = with(esoph, ncases + ncontrols))
gbm_fit <- fit(mf, model = GBMModel)
varimp(gbm_fit)




cleanEx()
nameEx("ModeledInput-methods")
### * ModeledInput-methods

flush(stderr()); flush(stdout())

### Name: ModeledInput
### Title: ModeledInput Classes
### Aliases: ModeledInput ModeledFrame ModeledRecipe ModeledInput.formula
###   ModeledInput.matrix ModeledInput.ModelFrame ModeledInput.recipe
###   ModeledInput.MLModel ModeledInput.MLModelFunction

### ** Examples

## Modeled model frame
mod_mf <- ModeledInput(sale_amount ~ ., data = ICHomes, model = GLMModel)
fit(mod_mf)

## Modeled recipe
library(recipes)

rec <- recipe(sale_amount ~ ., data = ICHomes)
mod_rec <- ModeledInput(rec, model = GLMModel)
fit(mod_rec)




cleanEx()
nameEx("NNetModel")
### * NNetModel

flush(stderr()); flush(stdout())

### Name: NNetModel
### Title: Neural Network Model
### Aliases: NNetModel

### ** Examples

fit(sale_amount ~ ., data = ICHomes, model = NNetModel)




cleanEx()
nameEx("NaiveBayesModel")
### * NaiveBayesModel

flush(stderr()); flush(stdout())

### Name: NaiveBayesModel
### Title: Naive Bayes Classifier Model
### Aliases: NaiveBayesModel

### ** Examples

fit(Species ~ ., data = iris, model = NaiveBayesModel)




cleanEx()
nameEx("PLSModel")
### * PLSModel

flush(stderr()); flush(stdout())

### Name: PLSModel
### Title: Partial Least Squares Model
### Aliases: PLSModel

### ** Examples

fit(sale_amount ~ ., data = ICHomes, model = PLSModel)




cleanEx()
nameEx("POLRModel")
### * POLRModel

flush(stderr()); flush(stdout())

### Name: POLRModel
### Title: Ordered Logistic or Probit Regression Model
### Aliases: POLRModel

### ** Examples

library(MASS)

df <- within(Boston,
             medv <- cut(medv,
                         breaks = c(0, 10, 15, 20, 25, 50),
                         ordered = TRUE))
fit(medv ~ ., data = df, model = POLRModel)




cleanEx()
nameEx("ParameterGrid")
### * ParameterGrid

flush(stderr()); flush(stdout())

### Name: ParameterGrid
### Title: Tuning Parameters Grid
### Aliases: ParameterGrid ParameterGrid.param ParameterGrid.list
###   ParameterGrid.parameters

### ** Examples

## GBMModel tuning parameters
library(dials)

grid <- ParameterGrid(
  n.trees = trees(),
  interaction.depth = tree_depth(),
  random = 5
)
TunedModel(GBMModel, grid = grid)




cleanEx()
nameEx("QDAModel")
### * QDAModel

flush(stderr()); flush(stdout())

### Name: QDAModel
### Title: Quadratic Discriminant Analysis Model
### Aliases: QDAModel

### ** Examples

fit(Species ~ ., data = iris, model = QDAModel)




cleanEx()
nameEx("RPartModel")
### * RPartModel

flush(stderr()); flush(stdout())

### Name: RPartModel
### Title: Recursive Partitioning and Regression Tree Models
### Aliases: RPartModel

### ** Examples

fit(Species ~ ., data = iris, model = RPartModel)




cleanEx()
nameEx("RandomForestModel")
### * RandomForestModel

flush(stderr()); flush(stdout())

### Name: RandomForestModel
### Title: Random Forest Model
### Aliases: RandomForestModel

### ** Examples

fit(sale_amount ~ ., data = ICHomes, model = RandomForestModel)




cleanEx()
nameEx("RangerModel")
### * RangerModel

flush(stderr()); flush(stdout())

### Name: RangerModel
### Title: Fast Random Forest Model
### Aliases: RangerModel

### ** Examples

fit(Species ~ ., data = iris, model = RangerModel)




cleanEx()
nameEx("SVMModel")
### * SVMModel

flush(stderr()); flush(stdout())

### Name: SVMModel
### Title: Support Vector Machine Models
### Aliases: SVMModel SVMANOVAModel SVMBesselModel SVMLaplaceModel
###   SVMLinearModel SVMPolyModel SVMRadialModel SVMSplineModel
###   SVMTanhModel

### ** Examples

fit(sale_amount ~ ., data = ICHomes, model = SVMRadialModel)




cleanEx()
nameEx("SelectedInput")
### * SelectedInput

flush(stderr()); flush(stdout())

### Name: SelectedInput
### Title: Selected Model Inputs
### Aliases: SelectedInput SelectedModelFrame SelectedModelRecipe
###   SelectedInput.formula SelectedInput.matrix SelectedInput.ModelFrame
###   SelectedInput.recipe SelectedInput.list

### ** Examples

## Selected model frame
sel_mf <- SelectedInput(
  sale_amount ~ sale_year + built + style + construction,
  sale_amount ~ sale_year + base_size + bedrooms + basement,
  data = ICHomes
)

fit(sel_mf, model = GLMModel)

## Selected recipe
library(recipes)
library(MASS)

rec1 <- recipe(medv ~ crim + zn + indus + chas + nox + rm, data = Boston)
rec2 <- recipe(medv ~ chas + nox + rm + age + dis + rad + tax, data = Boston)
sel_rec <- SelectedInput(rec1, rec2)

fit(sel_rec, model = GLMModel)




cleanEx()
nameEx("SelectedModel")
### * SelectedModel

flush(stderr()); flush(stdout())

### Name: SelectedModel
### Title: Selected Model
### Aliases: SelectedModel

### ** Examples

model_fit <- fit(sale_amount ~ ., data = ICHomes,
                 model = SelectedModel(GBMModel, GLMNetModel, SVMRadialModel))
(selected_model <- as.MLModel(model_fit))
summary(selected_model)




cleanEx()
nameEx("StackedModel")
### * StackedModel

flush(stderr()); flush(stdout())

### Name: StackedModel
### Title: Stacked Regression Model
### Aliases: StackedModel

### ** Examples

model <- StackedModel(GBMModel, SVMRadialModel, GLMNetModel(lambda = 0.01))
model_fit <- fit(sale_amount ~ ., data = ICHomes, model = model)
predict(model_fit, newdata = ICHomes)




cleanEx()
nameEx("SuperModel")
### * SuperModel

flush(stderr()); flush(stdout())

### Name: SuperModel
### Title: Super Learner Model
### Aliases: SuperModel

### ** Examples

model <- SuperModel(GBMModel, SVMRadialModel, GLMNetModel(lambda = 0.01))
model_fit <- fit(sale_amount ~ ., data = ICHomes, model = model)
predict(model_fit, newdata = ICHomes)




cleanEx()
nameEx("SurvRegModel")
### * SurvRegModel

flush(stderr()); flush(stdout())

### Name: SurvRegModel
### Title: Parametric Survival Model
### Aliases: SurvRegModel SurvRegStepAICModel

### ** Examples

library(survival)

fit(Surv(time, status) ~ ., data = veteran, model = SurvRegModel)




cleanEx()
nameEx("TreeModel")
### * TreeModel

flush(stderr()); flush(stdout())

### Name: TreeModel
### Title: Classification and Regression Tree Models
### Aliases: TreeModel

### ** Examples

fit(Species ~ ., data = iris, model = TreeModel)




cleanEx()
nameEx("TunedInput")
### * TunedInput

flush(stderr()); flush(stdout())

### Name: TunedInput
### Title: Tuned Model Inputs
### Aliases: TunedInput TunedModelRecipe TunedInput.recipe

### ** Examples

library(recipes)
library(MASS)

rec <- recipe(medv ~ ., data = Boston) %>%
  step_pca(all_numeric(), -all_outcomes(), id = "pca")

grid <- expand_steps(
  pca = list(num_comp = 1:2)
)

fit(TunedInput(rec, grid = grid), model = GLMModel)




cleanEx()
nameEx("TunedModel")
### * TunedModel

flush(stderr()); flush(stdout())

### Name: TunedModel
### Title: Tuned Model
### Aliases: TunedModel

### ** Examples





cleanEx()
nameEx("XGBModel")
### * XGBModel

flush(stderr()); flush(stdout())

### Name: XGBModel
### Title: Extreme Gradient Boosting Models
### Aliases: XGBModel XGBDARTModel XGBLinearModel XGBTreeModel

### ** Examples

model_fit <- fit(Species ~ ., data = iris, model = XGBTreeModel)
varimp(model_fit, metric = "Frequency", scale = FALSE)




cleanEx()
nameEx("calibration")
### * calibration

flush(stderr()); flush(stdout())

### Name: calibration
### Title: Model Calibration
### Aliases: calibration

### ** Examples

library(survival)

res <- resample(Surv(time, status) ~ ., data = veteran, model = GBMModel,
                control = CVControl(times = c(90, 180, 360)))
cal <- calibration(res)
plot(cal)




cleanEx()
nameEx("confusion")
### * confusion

flush(stderr()); flush(stdout())

### Name: confusion
### Title: Confusion Matrix
### Aliases: confusion ConfusionMatrix

### ** Examples

res <- resample(Species ~ ., data = iris, model = GBMModel)
(conf <- confusion(res))
plot(conf)




cleanEx()
nameEx("dependence")
### * dependence

flush(stderr()); flush(stdout())

### Name: dependence
### Title: Partial Dependence
### Aliases: dependence

### ** Examples

gbm_fit <- fit(Species ~ ., data = iris, model = GBMModel)
(pd <- dependence(gbm_fit, select = c(Petal.Length, Petal.Width)))
plot(pd)




cleanEx()
nameEx("diff-methods")
### * diff-methods

flush(stderr()); flush(stdout())

### Name: diff
### Title: Model Performance Differences
### Aliases: diff diff.MLModel diff.Performance diff.Resamples

### ** Examples

## Survival response example
library(survival)

fo <- Surv(time, status) ~ .
control <- CVControl()

gbm_res1 <- resample(fo, data = veteran, GBMModel(n.trees = 25), control)
gbm_res2 <- resample(fo, data = veteran, GBMModel(n.trees = 50), control)
gbm_res3 <- resample(fo, data = veteran, GBMModel(n.trees = 100), control)

res <- c(GBM1 = gbm_res1, GBM2 = gbm_res2, GBM3 = gbm_res3)
res_diff <- diff(res)
summary(res_diff)
plot(res_diff)




cleanEx()
nameEx("dot-")
### * dot-

flush(stderr()); flush(stdout())

### Name: .
### Title: Quote Operator
### Aliases: .

### ** Examples

## Stepwise variable selection with BIC
glm_fit <- fit(sale_amount ~ ., ICHomes, GLMStepAICModel(k = .(log(nobs))))
varimp(glm_fit)




cleanEx()
nameEx("expand_model")
### * expand_model

flush(stderr()); flush(stdout())

### Name: expand_model
### Title: Model Expansion Over Tuning Parameters
### Aliases: expand_model

### ** Examples

library(MASS)

models <- expand_model(GBMModel, n.trees = c(50, 100),
                                 interaction.depth = 1:2)

fit(medv ~ ., data = Boston, model = SelectedModel(models))




cleanEx()
nameEx("expand_params")
### * expand_params

flush(stderr()); flush(stdout())

### Name: expand_params
### Title: Model Parameters Expansion
### Aliases: expand_params

### ** Examples

library(MASS)

grid <- expand_params(
  n.trees = c(50, 100),
  interaction.depth = 1:2
)

fit(medv ~ ., data = Boston, model = TunedModel(GBMModel, grid = grid))




cleanEx()
nameEx("expand_steps")
### * expand_steps

flush(stderr()); flush(stdout())

### Name: expand_steps
### Title: Recipe Step Parameters Expansion
### Aliases: expand_steps

### ** Examples

library(recipes)
library(MASS)

rec <- recipe(medv ~ ., data = Boston) %>%
  step_corr(all_numeric(), -all_outcomes(), id = "corr") %>%
  step_pca(all_numeric(), -all_outcomes(), id = "pca")

expand_steps(
  corr = list(threshold = c(0.8, 0.9),
              method = c("pearson", "spearman")),
  pca = list(num_comp = 1:3)
)




cleanEx()
nameEx("fit-methods")
### * fit-methods

flush(stderr()); flush(stdout())

### Name: fit
### Title: Model Fitting
### Aliases: fit fit.formula fit.matrix fit.ModelFrame fit.recipe
###   fit.MLModel fit.MLModelFunction

### ** Examples

## Survival response example
library(survival)

gbm_fit <- fit(Surv(time, status) ~ ., data = veteran, model = GBMModel)
varimp(gbm_fit)




cleanEx()
nameEx("lift")
### * lift

flush(stderr()); flush(stdout())

### Name: lift
### Title: Model Lift Curves
### Aliases: lift

### ** Examples

library(MASS)

res <- resample(type ~ ., data = Pima.tr, model = GBMModel)
lf <- lift(res)
plot(lf)




cleanEx()
nameEx("metricinfo")
### * metricinfo

flush(stderr()); flush(stdout())

### Name: metricinfo
### Title: Display Performance Metric Information
### Aliases: metricinfo

### ** Examples

## All metrics
metricinfo()

## Metrics by observed and predicted response types
names(metricinfo(factor(0)))
names(metricinfo(factor(0), factor(0)))
names(metricinfo(factor(0), matrix(0)))
names(metricinfo(factor(0), numeric(0)))

## Metric-specific information
metricinfo(auc)




cleanEx()
nameEx("modelinfo")
### * modelinfo

flush(stderr()); flush(stdout())

### Name: modelinfo
### Title: Display Model Information
### Aliases: modelinfo

### ** Examples

## All models
modelinfo()

## Models by response types
names(modelinfo(factor(0)))
names(modelinfo(factor(0), numeric(0)))

## Model-specific information
modelinfo(GBMModel)




cleanEx()
nameEx("performance")
### * performance

flush(stderr()); flush(stdout())

### Name: performance
### Title: Model Performance Metrics
### Aliases: performance performance.BinomialVariate performance.factor
###   performance.matrix performance.numeric performance.Surv
###   performance.ConfusionList performance.ConfusionMatrix
###   performance.Resamples

### ** Examples

res <- resample(Species ~ ., data = iris, model = GBMModel)
(perf <- performance(res))
summary(perf)
plot(perf)

## Survival response example
library(survival)

gbm_fit <- fit(Surv(time, status) ~ ., data = veteran, model = GBMModel)

obs <- response(gbm_fit, newdata = veteran)
pred <- predict(gbm_fit, newdata = veteran, type = "prob")
performance(obs, pred)




cleanEx()
nameEx("performance_curve")
### * performance_curve

flush(stderr()); flush(stdout())

### Name: performance_curve
### Title: Model Performance Curves
### Aliases: performance_curve curves performance_curve.default
###   performance_curve.Resamples

### ** Examples

library(MASS)

res <- resample(type ~ ., data = Pima.tr, model = GBMModel)

## ROC curve
roc <- performance_curve(res)
plot(roc)
auc(roc)




cleanEx()
nameEx("plot-methods")
### * plot-methods

flush(stderr()); flush(stdout())

### Name: plot
### Title: Model Performance Plots
### Aliases: plot plot.Calibration plot.ConfusionList plot.ConfusionMatrix
###   plot.LiftCurve plot.MLModel plot.PartialDependence plot.Performance
###   plot.PerformanceCurve plot.Resamples plot.VarImp

### ** Examples

## Factor response example

fo <- Species ~ .
control <- CVControl()

gbm_fit <- fit(fo, data = iris, model = GBMModel, control = control)
plot(varimp(gbm_fit))

gbm_res1 <- resample(fo, iris, GBMModel(n.trees = 25), control)
gbm_res2 <- resample(fo, iris, GBMModel(n.trees = 50), control)
gbm_res3 <- resample(fo, iris, GBMModel(n.trees = 100), control)
plot(gbm_res3)

res <- c(GBM1 = gbm_res1, GBM2 = gbm_res2, GBM3 = gbm_res3)
plot(res)




cleanEx()
nameEx("predict")
### * predict

flush(stderr()); flush(stdout())

### Name: predict
### Title: Model Prediction
### Aliases: predict predict.MLModelFit

### ** Examples

## Survival response example
library(survival)

gbm_fit <- fit(Surv(time, status) ~ ., data = veteran, model = GBMModel)
predict(gbm_fit, newdata = veteran, times = c(90, 180, 360), type = "prob")




cleanEx()
nameEx("recipe_roles")
### * recipe_roles

flush(stderr()); flush(stdout())

### Name: recipe_roles
### Title: Set Recipe Roles
### Aliases: recipe_roles role_binom role_case role_pred role_surv

### ** Examples

library(survival)
library(recipes)

rec <- recipe(time + status ~ ., data = veteran) %>%
  role_surv(time = time, event = status) %>%
  role_case(stratum = status)

(res <- resample(rec, model = CoxModel))
summary(res)




cleanEx()
nameEx("resample-methods")
### * resample-methods

flush(stderr()); flush(stdout())

### Name: resample
### Title: Resample Estimation of Model Performance
### Aliases: resample resample.formula resample.matrix resample.ModelFrame
###   resample.recipe resample.MLModel resample.MLModelFunction

### ** Examples

## Factor response example

fo <- Species ~ .
control <- CVControl()

gbm_res1 <- resample(fo, iris, GBMModel(n.trees = 25), control)
gbm_res2 <- resample(fo, iris, GBMModel(n.trees = 50), control)
gbm_res3 <- resample(fo, iris, GBMModel(n.trees = 100), control)

summary(gbm_res1)
plot(gbm_res1)

res <- c(GBM1 = gbm_res1, GBM2 = gbm_res2, GBM3 = gbm_res3)
summary(res)
plot(res)




cleanEx()
nameEx("response-methods")
### * response-methods

flush(stderr()); flush(stdout())

### Name: response
### Title: Extract Response Variable
### Aliases: response response.MLModelFit response.ModelFrame
###   response.recipe

### ** Examples

## Survival response example
library(survival)

mf <- ModelFrame(Surv(time, status) ~ ., data = veteran)
response(mf)




cleanEx()
nameEx("settings")
### * settings

flush(stderr()); flush(stdout())

### Name: settings
### Title: MachineShop Settings
### Aliases: settings

### ** Examples

## View all current settings
settings()

## Change settings
presets <- settings(control = "BootControl", grid = 10)

## View one setting
settings("control")

## View multiple settings
settings("control", "grid")

## Restore the previous settings
settings(presets)




cleanEx()
nameEx("step_kmeans")
### * step_kmeans

flush(stderr()); flush(stdout())

### Name: step_kmeans
### Title: K-Means Clustering Variable Reduction
### Aliases: step_kmeans tidy.step_kmeans tunable.step_kmeans

### ** Examples

library(recipes)

rec <- recipe(rating ~ ., data = attitude)
kmeans_rec <- rec %>%
  step_kmeans(all_predictors(), k = 3)
kmeans_prep <- prep(kmeans_rec, training = attitude)
kmeans_data <- bake(kmeans_prep, attitude)

pairs(kmeans_data, lower.panel = NULL)

tidy(kmeans_rec, number = 1)
tidy(kmeans_prep, number = 1)




cleanEx()
nameEx("step_kmedoids")
### * step_kmedoids

flush(stderr()); flush(stdout())

### Name: step_kmedoids
### Title: K-Medoids Clustering Variable Selection
### Aliases: step_kmedoids tidy.step_kmedoids tunable.step_kmedoids

### ** Examples

library(recipes)

rec <- recipe(rating ~ ., data = attitude)
kmedoids_rec <- rec %>%
  step_kmedoids(all_predictors(), k = 3)
kmedoids_prep <- prep(kmedoids_rec, training = attitude)
kmedoids_data <- bake(kmedoids_prep, attitude)

pairs(kmedoids_data, lower.panel = NULL)

tidy(kmedoids_rec, number = 1)
tidy(kmedoids_prep, number = 1)




cleanEx()
nameEx("step_spca")
### * step_spca

flush(stderr()); flush(stdout())

### Name: step_spca
### Title: Sparse Principal Components Analysis Variable Reduction
### Aliases: step_spca tidy.step_spca tunable.step_spca

### ** Examples

library(recipes)

rec <- recipe(rating ~ ., data = attitude)
spca_rec <- rec %>%
  step_spca(all_predictors(), num_comp = 5, sparsity = 1)
spca_prep <- prep(spca_rec, training = attitude)
spca_data <- bake(spca_prep, attitude)

pairs(spca_data, lower.panel = NULL)

tidy(spca_rec, number = 1)
tidy(spca_prep, number = 1)




cleanEx()
nameEx("summary-methods")
### * summary-methods

flush(stderr()); flush(stdout())

### Name: summary
### Title: Model Performance Summaries
### Aliases: summary summary.ConfusionList summary.ConfusionMatrix
###   summary.MLModel summary.Performance summary.PerformanceCurve
###   summary.Resamples

### ** Examples

## Factor response example

fo <- Species ~ .
control <- CVControl()

gbm_res1 <- resample(fo, iris, GBMModel(n.trees = 25), control)
gbm_res2 <- resample(fo, iris, GBMModel(n.trees = 50), control)
gbm_res3 <- resample(fo, iris, GBMModel(n.trees = 100), control)
summary(gbm_res3)

res <- c(GBM1 = gbm_res1, GBM2 = gbm_res2, GBM3 = gbm_res3)
summary(res)




cleanEx()
nameEx("t.test")
### * t.test

flush(stderr()); flush(stdout())

### Name: t.test
### Title: Paired t-Tests for Model Comparisons
### Aliases: t.test t.test.PerformanceDiff

### ** Examples

## Numeric response example
fo <- sale_amount ~ .
control <- CVControl()

gbm_res1 <- resample(fo, ICHomes, GBMModel(n.trees = 25), control)
gbm_res2 <- resample(fo, ICHomes, GBMModel(n.trees = 50), control)
gbm_res3 <- resample(fo, ICHomes, GBMModel(n.trees = 100), control)

res <- c(GBM1 = gbm_res1, GBM2 = gbm_res2, GBM3 = gbm_res3)
res_diff <- diff(res)
t.test(res_diff)




cleanEx()
nameEx("varimp")
### * varimp

flush(stderr()); flush(stdout())

### Name: varimp
### Title: Variable Importance
### Aliases: varimp

### ** Examples

## Survival response example
library(survival)

gbm_fit <- fit(Surv(time, status) ~ ., data = veteran, model = GBMModel)
(vi <- varimp(gbm_fit))
plot(vi)




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
