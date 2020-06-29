pkgname <- "insight"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('insight')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("all_models_equal")
### * all_models_equal

flush(stderr()); flush(stdout())

### Name: all_models_equal
### Title: Checks if all objects are models of same class
### Aliases: all_models_equal all_models_same_class

### ** Examples

library(lme4)
data(mtcars)
data(sleepstudy)

m1 <- lm(mpg ~ wt + cyl + vs, data = mtcars)
m2 <- lm(mpg ~ wt + cyl, data = mtcars)
m3 <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
m4 <- glm(formula = vs ~ wt, family = binomial(), data = mtcars)

all_models_same_class(m1, m2)
all_models_same_class(m1, m2, m3)
all_models_same_class(m1, m4, m2, m3, verbose = TRUE)
all_models_same_class(m1, m4, mtcars, m2, m3, verbose = TRUE)



cleanEx()
nameEx("clean_names")
### * clean_names

flush(stderr()); flush(stdout())

### Name: clean_names
### Title: Get clean names of model terms
### Aliases: clean_names

### ** Examples

# example from ?stats::glm
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- as.numeric(gl(3, 1, 9))
treatment <- gl(3, 3)
m <- glm(counts ~ log(outcome) + as.factor(treatment), family = poisson())
clean_names(m)

# difference "clean_names()" and "find_variables()"
library(lme4)
m <- glmer(
  cbind(incidence, size - incidence) ~ period + (1 | herd),
  data = cbpp,
  family = binomial
)

clean_names(m)
find_variables(m)
find_variables(m, flatten = TRUE)



cleanEx()
nameEx("clean_parameters")
### * clean_parameters

flush(stderr()); flush(stdout())

### Name: clean_parameters
### Title: Get clean names of model parameters
### Aliases: clean_parameters

### ** Examples

## Not run: 
##D library(brms)
##D model <- download_model("brms_zi_2")
##D clean_parameters(model)
## End(Not run)



cleanEx()
nameEx("color_if")
### * color_if

flush(stderr()); flush(stdout())

### Name: color_if
### Title: Color-formatting for data columns based on condition
### Aliases: color_if colour_if

### ** Examples

# all values in Sepal.Length larger than 5 in green, all remaining in red
x <- color_if(iris[1:10, ], columns = "Sepal.Length", predicate = `>`, value = 5)
x
cat(x$Sepal.Length)

# all levels "setosa" in Species in green, all remaining in red
x <- color_if(iris, columns = "Species", predicate = `==`, value = "setosa")
cat(x$Species)

# own function, argument "value" not needed here
p <- function(x, y) {
  x >= 4.9 & x <= 5.1
}
# all values in Sepal.Length between 4.9 and 5.1 in green, all remaining in red
x <- color_if(iris[1:10, ], columns = "Sepal.Length", predicate = p)
cat(x$Sepal.Length)



cleanEx()
nameEx("find_algorithm")
### * find_algorithm

flush(stderr()); flush(stdout())

### Name: find_algorithm
### Title: Find sampling algorithm and optimizers
### Aliases: find_algorithm

### ** Examples

library(lme4)
data(sleepstudy)
m <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
find_algorithm(m)
## Not run: 
##D library(rstanarm)
##D m <- stan_lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
##D find_algorithm(m)
## End(Not run)




cleanEx()
nameEx("find_formula")
### * find_formula

flush(stderr()); flush(stdout())

### Name: find_formula
### Title: Find model formula
### Aliases: find_formula

### ** Examples

data(mtcars)
m <- lm(mpg ~ wt + cyl + vs, data = mtcars)
find_formula(m)



cleanEx()
nameEx("find_interactions")
### * find_interactions

flush(stderr()); flush(stdout())

### Name: find_interactions
### Title: Find interaction terms from models
### Aliases: find_interactions

### ** Examples

data(mtcars)

m <- lm(mpg ~ wt + cyl + vs, data = mtcars)
find_interactions(m)

m <- lm(mpg ~ wt * cyl + vs * hp * gear + carb, data = mtcars)
find_interactions(m)



cleanEx()
nameEx("find_parameters")
### * find_parameters

flush(stderr()); flush(stdout())

### Name: find_parameters
### Title: Find names of model parameters
### Aliases: find_parameters find_parameters.betamfx
###   find_parameters.logitmfx find_parameters.gam find_parameters.merMod
###   find_parameters.zeroinfl find_parameters.BFBayesFactor
###   find_parameters.brmsfit find_parameters.bayesx
###   find_parameters.stanreg find_parameters.sim.merMod
###   find_parameters.averaging

### ** Examples

data(mtcars)
m <- lm(mpg ~ wt + cyl + vs, data = mtcars)
find_parameters(m)



cleanEx()
nameEx("find_predictors")
### * find_predictors

flush(stderr()); flush(stdout())

### Name: find_predictors
### Title: Find names of model predictors
### Aliases: find_predictors

### ** Examples

data(mtcars)
m <- lm(mpg ~ wt + cyl + vs, data = mtcars)
find_predictors(m)



cleanEx()
nameEx("find_random")
### * find_random

flush(stderr()); flush(stdout())

### Name: find_random
### Title: Find names of random effects
### Aliases: find_random

### ** Examples

library(lme4)
data(sleepstudy)
sleepstudy$mygrp <- sample(1:5, size = 180, replace = TRUE)
sleepstudy$mysubgrp <- NA
for (i in 1:5) {
  filter_group <- sleepstudy$mygrp == i
  sleepstudy$mysubgrp[filter_group] <-
    sample(1:30, size = sum(filter_group), replace = TRUE)
}

m <- lmer(
  Reaction ~ Days + (1 | mygrp / mysubgrp) + (1 | Subject),
  data = sleepstudy
)

find_random(m)
find_random(m, split_nested = TRUE)



cleanEx()
nameEx("find_random_slopes")
### * find_random_slopes

flush(stderr()); flush(stdout())

### Name: find_random_slopes
### Title: Find names of random slopes
### Aliases: find_random_slopes

### ** Examples

library(lme4)
data(sleepstudy)

m <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)
find_random_slopes(m)



cleanEx()
nameEx("find_response")
### * find_response

flush(stderr()); flush(stdout())

### Name: find_response
### Title: Find name of the response variable
### Aliases: find_response

### ** Examples

library(lme4)
data(cbpp)
cbpp$trials <- cbpp$size - cbpp$incidence
m <- glm(cbind(incidence, trials) ~ period, data = cbpp, family = binomial)

find_response(m, combine = TRUE)
find_response(m, combine = FALSE)



cleanEx()
nameEx("find_statistic")
### * find_statistic

flush(stderr()); flush(stdout())

### Name: find_statistic
### Title: Find statistic for model
### Aliases: find_statistic

### ** Examples

# regression model object
data(mtcars)
m <- lm(mpg ~ wt + cyl + vs, data = mtcars)
find_statistic(m)



cleanEx()
nameEx("find_terms")
### * find_terms

flush(stderr()); flush(stdout())

### Name: find_terms
### Title: Find all model terms
### Aliases: find_terms

### ** Examples

library(lme4)
data(sleepstudy)
m <- lmer(
  log(Reaction) ~ Days + I(Days^2) + (1 + Days + exp(Days) | Subject),
  data = sleepstudy
)

find_terms(m)



cleanEx()
nameEx("find_variables")
### * find_variables

flush(stderr()); flush(stdout())

### Name: find_variables
### Title: Find names of all variables
### Aliases: find_variables

### ** Examples

library(lme4)
data(cbpp)
data(sleepstudy)
# some data preparation...
cbpp$trials <- cbpp$size - cbpp$incidence
sleepstudy$mygrp <- sample(1:5, size = 180, replace = TRUE)
sleepstudy$mysubgrp <- NA
for (i in 1:5) {
  filter_group <- sleepstudy$mygrp == i
  sleepstudy$mysubgrp[filter_group] <-
    sample(1:30, size = sum(filter_group), replace = TRUE)
}

m1 <- glmer(
  cbind(incidence, size - incidence) ~ period + (1 | herd),
  data = cbpp,
  family = binomial
)
find_variables(m1)

m2 <- lmer(
  Reaction ~ Days + (1 | mygrp / mysubgrp) + (1 | Subject),
  data = sleepstudy
)
find_variables(m2)
find_variables(m2, flatten = TRUE)



cleanEx()
nameEx("find_weights")
### * find_weights

flush(stderr()); flush(stdout())

### Name: find_weights
### Title: Find names of model weights
### Aliases: find_weights

### ** Examples

data(mtcars)
mtcars$weight <- rnorm(nrow(mtcars), 1, .3)
m <- lm(mpg ~ wt + cyl + vs, data = mtcars, weights = weight)
find_weights(m)



cleanEx()
nameEx("format_ci")
### * format_ci

flush(stderr()); flush(stdout())

### Name: format_ci
### Title: Confidence/Credible Interval (CI) Formatting
### Aliases: format_ci

### ** Examples

format_ci(1.20, 3.57, ci = 0.90)
format_ci(1.20, 3.57, ci = NULL)
format_ci(1.20, 3.57, ci = NULL, brackets = FALSE)
format_ci(c(1.205645, 23.4), c(3.57, -1.35), ci = 0.90)
format_ci(c(1.20, NA, NA), c(3.57, -1.35, NA), ci = 0.90)

# automatic alignment of width, useful for printing multiple CIs in columns
x <- format_ci(c(1.205, 23.4, 100.43), c(3.57, -13.35, 9.4))
cat(x, sep = "\n")

x <- format_ci(c(1.205, 23.4, 100.43), c(3.57, -13.35, 9.4), width = "auto")
cat(x, sep = "\n")



cleanEx()
nameEx("format_table")
### * format_table

flush(stderr()); flush(stdout())

### Name: format_table
### Title: Dataframe and Tables Pretty Formatting
### Aliases: format_table

### ** Examples

cat(format_table(iris))
cat(format_table(iris, sep = " ", header = "*", digits = 1))



cleanEx()
nameEx("format_value")
### * format_value

flush(stderr()); flush(stdout())

### Name: format_value
### Title: Numeric Values Formatting
### Aliases: format_value

### ** Examples

format_value(1.20)
format_value(1.2)
format_value(1.2012313)
format_value(c(0.0045, 234, -23))
format_value(c(0.0045, .12, .34))
format_value(c(0.0045, .12, .34), as_percent = TRUE)

format_value(as.factor(c("A", "B", "A")))
format_value(iris$Species)

format_value(3)
format_value(3, protect_integers = TRUE)

format_value(iris)



cleanEx()
nameEx("get_data")
### * get_data

flush(stderr()); flush(stdout())

### Name: get_data
### Title: Get the data that was used to fit the model
### Aliases: get_data get_data.gee get_data.rqss get_data.hurdle
###   get_data.zcpglm get_data.glmmTMB get_data.merMod get_data.glmmadmb
###   get_data.rlmerMod get_data.clmm get_data.mixed get_data.lme
###   get_data.MixMod get_data.brmsfit get_data.stanreg get_data.MCMCglmm

### ** Examples

data(cbpp, package = "lme4")
cbpp$trials <- cbpp$size - cbpp$incidence
m <- glm(cbind(incidence, trials) ~ period, data = cbpp, family = binomial)
head(get_data(m))



cleanEx()
nameEx("get_parameters")
### * get_parameters

flush(stderr()); flush(stdout())

### Name: get_parameters
### Title: Get model parameters
### Aliases: get_parameters get_parameters.betamfx get_parameters.logitmfx
###   get_parameters.averaging get_parameters.betareg
###   get_parameters.DirichletRegModel get_parameters.clm2
###   get_parameters.coxme get_parameters.merMod get_parameters.mixed
###   get_parameters.glmmTMB get_parameters.BBmm get_parameters.glimML
###   get_parameters.gam get_parameters.zeroinfl get_parameters.zcpglm
###   get_parameters.MCMCglmm get_parameters.BFBayesFactor
###   get_parameters.stanmvreg get_parameters.brmsfit
###   get_parameters.stanreg get_parameters.sim.merMod

### ** Examples

data(mtcars)
m <- lm(mpg ~ wt + cyl + vs, data = mtcars)
get_parameters(m)



cleanEx()
nameEx("get_predictors")
### * get_predictors

flush(stderr()); flush(stdout())

### Name: get_predictors
### Title: Get the data from model predictors
### Aliases: get_predictors

### ** Examples

m <- lm(mpg ~ wt + cyl + vs, data = mtcars)
head(get_predictors(m))



cleanEx()
nameEx("get_priors")
### * get_priors

flush(stderr()); flush(stdout())

### Name: get_priors
### Title: Get summary of priors used for a model
### Aliases: get_priors

### ** Examples

## Not run: 
##D library(rstanarm)
##D model <- stan_glm(Sepal.Width ~ Species * Petal.Length, data = iris)
##D get_priors(model)
## End(Not run)




cleanEx()
nameEx("get_random")
### * get_random

flush(stderr()); flush(stdout())

### Name: get_random
### Title: Get the data from random effects
### Aliases: get_random

### ** Examples

library(lme4)
data(sleepstudy)
# prepare some data...
sleepstudy$mygrp <- sample(1:5, size = 180, replace = TRUE)
sleepstudy$mysubgrp <- NA
for (i in 1:5) {
  filter_group <- sleepstudy$mygrp == i
  sleepstudy$mysubgrp[filter_group] <-
    sample(1:30, size = sum(filter_group), replace = TRUE)
}

m <- lmer(
  Reaction ~ Days + (1 | mygrp / mysubgrp) + (1 | Subject),
  data = sleepstudy
)

head(get_random(m))



cleanEx()
nameEx("get_response")
### * get_response

flush(stderr()); flush(stdout())

### Name: get_response
### Title: Get the values from the response variable
### Aliases: get_response

### ** Examples

library(lme4)
data(cbpp)
data(mtcars)
cbpp$trials <- cbpp$size - cbpp$incidence

m <- glm(cbind(incidence, trials) ~ period, data = cbpp, family = binomial)
head(get_response(m))
get_response(m, select = "incidence")

m <- lm(mpg ~ wt + cyl + vs, data = mtcars)
get_response(m)



cleanEx()
nameEx("get_statistic")
### * get_statistic

flush(stderr()); flush(stdout())

### Name: get_statistic
### Title: Get statistic associated with estimates
### Aliases: get_statistic get_statistic.default get_statistic.glmmTMB
###   get_statistic.clm2 get_statistic.betamfx get_statistic.logitmfx
###   get_statistic.emmGrid get_statistic.gee get_statistic.betareg
###   get_statistic.DirichletRegModel

### ** Examples

data(mtcars)
m <- lm(mpg ~ wt + cyl + vs, data = mtcars)
get_statistic(m)



cleanEx()
nameEx("get_varcov")
### * get_varcov

flush(stderr()); flush(stdout())

### Name: get_varcov
### Title: Get variance-covariance matrix from models
### Aliases: get_varcov get_varcov.betareg get_varcov.DirichletRegModel
###   get_varcov.clm2 get_varcov.truncreg get_varcov.gamlss
###   get_varcov.hurdle get_varcov.zcpglm get_varcov.MixMod
###   get_varcov.glmmTMB get_varcov.brmsfit get_varcov.betamfx
###   get_varcov.mixor

### ** Examples

data(mtcars)
m <- lm(mpg ~ wt + cyl + vs, data = mtcars)
get_varcov(m)



cleanEx()
nameEx("get_variance")
### * get_variance

flush(stderr()); flush(stdout())

### Name: get_variance
### Title: Get variance components from random effects models
### Aliases: get_variance get_variance_residual get_variance_fixed
###   get_variance_random get_variance_distribution get_variance_dispersion
###   get_variance_intercept get_variance_slope
###   get_correlation_slope_intercept

### ** Examples

## Not run: 
##D library(lme4)
##D data(sleepstudy)
##D m <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)
##D 
##D get_variance(m)
##D get_variance_fixed(m)
##D get_variance_residual(m)
## End(Not run)



cleanEx()
nameEx("get_weights")
### * get_weights

flush(stderr()); flush(stdout())

### Name: get_weights
### Title: Get the values from model weights
### Aliases: get_weights

### ** Examples

data(mtcars)
mtcars$weight <- rnorm(nrow(mtcars), 1, .3)
m <- lm(mpg ~ wt + cyl + vs, data = mtcars, weights = weight)
get_weights(m)



cleanEx()
nameEx("has_intercept")
### * has_intercept

flush(stderr()); flush(stdout())

### Name: has_intercept
### Title: Checks if model has an intercept
### Aliases: has_intercept

### ** Examples

model <- lm(mpg ~ 0 + gear, data = mtcars)
has_intercept(model)

model <- lm(mpg ~ gear, data = mtcars)
has_intercept(model)

library(lme4)
model <- lmer(Reaction ~ 0 + Days + (Days | Subject), data = sleepstudy)
has_intercept(model)

model <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
has_intercept(model)



cleanEx()
nameEx("is_model")
### * is_model

flush(stderr()); flush(stdout())

### Name: is_model
### Title: Checks if an object is a regression model object
### Aliases: is_model

### ** Examples

data(mtcars)
m <- lm(mpg ~ wt + cyl + vs, data = mtcars)

is_model(m)
is_model(mtcars)



cleanEx()
nameEx("is_model_supported")
### * is_model_supported

flush(stderr()); flush(stdout())

### Name: is_model_supported
### Title: Checks if an object is a regression model object supported in
###   'insight' package.
### Aliases: is_model_supported supported_models

### ** Examples

data(mtcars)
m <- lm(mpg ~ wt + cyl + vs, data = mtcars)

is_model_supported(m)
is_model_supported(mtcars)



cleanEx()
nameEx("is_multivariate")
### * is_multivariate

flush(stderr()); flush(stdout())

### Name: is_multivariate
### Title: Checks if an object stems from a multivariate response model
### Aliases: is_multivariate

### ** Examples

## Not run: 
##D library(rstanarm)
##D data("pbcLong")
##D model <- stan_mvmer(
##D   formula = list(
##D     logBili ~ year + (1 | id),
##D     albumin ~ sex + year + (year | id)
##D   ),
##D   data = pbcLong,
##D   chains = 1, cores = 1, seed = 12345, iter = 1000
##D )
##D 
##D f <- find_formula(model)
##D is_multivariate(model)
##D is_multivariate(f)
## End(Not run)



cleanEx()
nameEx("is_nullmodel")
### * is_nullmodel

flush(stderr()); flush(stdout())

### Name: is_nullmodel
### Title: Checks if model is a null-model (intercept-only)
### Aliases: is_nullmodel

### ** Examples

model <- lm(mpg ~ 1, data = mtcars)
is_nullmodel(model)

model <- lm(mpg ~ gear, data = mtcars)
is_nullmodel(model)

library(lme4)
model <- lmer(Reaction ~ 1 + (Days | Subject), data = sleepstudy)
is_nullmodel(model)

model <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
is_nullmodel(model)



cleanEx()
nameEx("link_function")
### * link_function

flush(stderr()); flush(stdout())

### Name: link_function
### Title: Get link-function from model object
### Aliases: link_function link_function.betamfx link_function.gamlss
###   link_function.betareg link_function.DirichletRegModel

### ** Examples

# example from ?stats::glm
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
m <- glm(counts ~ outcome + treatment, family = poisson())

link_function(m)(.3)
# same as
log(.3)



cleanEx()
nameEx("link_inverse")
### * link_inverse

flush(stderr()); flush(stdout())

### Name: link_inverse
### Title: Get link-inverse function from model object
### Aliases: link_inverse link_inverse.betareg
###   link_inverse.DirichletRegModel link_inverse.betamfx
###   link_inverse.gamlss

### ** Examples

# example from ?stats::glm
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
m <- glm(counts ~ outcome + treatment, family = poisson())

link_inverse(m)(.3)
# same as
exp(.3)



cleanEx()
nameEx("model_info")
### * model_info

flush(stderr()); flush(stdout())

### Name: model_info
### Title: Access information from model objects
### Aliases: model_info model_info.default

### ** Examples

ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive = 20 - numdead)
dat <- data.frame(ldose, sex, SF, stringsAsFactors = FALSE)
m <- glm(SF ~ sex * ldose, family = binomial)

model_info(m)
## Not run: 
##D library(glmmTMB)
##D data("Salamanders")
##D m <- glmmTMB(
##D   count ~ spp + cover + mined + (1 | site),
##D   ziformula = ~ spp + mined,
##D   dispformula = ~DOY,
##D   data = Salamanders,
##D   family = nbinom2
##D )
## End(Not run)

model_info(m)



cleanEx()
nameEx("n_obs")
### * n_obs

flush(stderr()); flush(stdout())

### Name: n_obs
### Title: Get number of observations from a model
### Aliases: n_obs n_obs.svyolr n_obs.stanmvreg

### ** Examples

data(mtcars)
m <- lm(mpg ~ wt + cyl + vs, data = mtcars)
n_obs(m)



cleanEx()
nameEx("null_model")
### * null_model

flush(stderr()); flush(stdout())

### Name: null_model
### Title: Compute intercept-only model for regression models
### Aliases: null_model

### ** Examples

if (require("lme4")) {
  data(sleepstudy)
  m <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)
  summary(m)
  summary(null_model(m))
}



cleanEx()
nameEx("print_color")
### * print_color

flush(stderr()); flush(stdout())

### Name: print_color
### Title: Coloured console output
### Aliases: print_color print_colour

### ** Examples

print_color("I'm blue dabedi dabedei", "blue")



cleanEx()
nameEx("print_parameters")
### * print_parameters

flush(stderr()); flush(stdout())

### Name: print_parameters
### Title: Prepare summary statistics of model parameters for printing
### Aliases: print_parameters

### ** Examples

## Not run: 
##D library(bayestestR)
##D model <- download_model("brms_zi_2")
##D x <- hdi(model, effects = "all", component = "all")
##D 
##D # hdi() returns a data frame; here we use only the informaton on
##D # parameter names and HDI values
##D tmp <- as.data.frame(x)[, 1:4]
##D tmp
##D 
##D # Based on the "split_by" argument, we get a list of data frames that
##D # is split into several parts that reflect the model components.
##D print_parameters(model, tmp)
##D 
##D # This is the standard print()-method for "bayestestR::hdi"-objects.
##D # For printing methods, it is easy to print complex summary statistics
##D # in a clean way to the console by splitting the information into
##D # different model components.
##D x
## End(Not run)




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
