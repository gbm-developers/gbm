pkgname <- "vip"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('vip')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("gen_friedman")
### * gen_friedman

flush(stderr()); flush(stdout())

### Name: gen_friedman
### Title: Friedman benchmark data
### Aliases: gen_friedman

### ** Examples

gen_friedman()



cleanEx()
nameEx("list_metrics")
### * list_metrics

flush(stderr()); flush(stdout())

### Name: list_metrics
### Title: List metrics
### Aliases: list_metrics

### ** Examples

(metrics <- list_metrics())
metrics[metrics$Task == "Multiclass classification", ]



cleanEx()
nameEx("metrics")
### * metrics

flush(stderr()); flush(stdout())

### Name: metric_mse
### Title: Model metrics
### Aliases: metric_mse metric_rmse metric_sse metric_mae metric_rsquared
###   metric_accuracy metric_error metric_auc metric_logLoss metric_mauc

### ** Examples

x <- rnorm(10)
y <- rnorm(10)
metric_mse(x, y)
metric_rsquared(x, y)



cleanEx()
nameEx("vi")
### * vi

flush(stderr()); flush(stdout())

### Name: vi
### Title: Variable importance
### Aliases: vi vi.default vi.model_fit vi.WrappedModel vi.Learner

### ** Examples

#
# A projection pursuit regression example
#

# Load the sample data
data(mtcars)

# Fit a projection pursuit regression model
mtcars.ppr <- ppr(mpg ~ ., data = mtcars, nterms = 1)

# Compute variable importance scores
vi(mtcars.ppr, method = "firm", ice = TRUE)
vi(mtcars.ppr, method = "firm", ice = TRUE,
   var_fun = list("con" = mad, "cat" = function(x) diff(range(x)) / 4))

# Plot variable importance scores
vip(mtcars.ppr, method = "firm", ice = TRUE)



cleanEx()
nameEx("vi_permute")
### * vi_permute

flush(stderr()); flush(stdout())

### Name: vi_permute
### Title: Permutation-based variable importance
### Aliases: vi_permute vi_permute.default

### ** Examples

## Not run: 
##D # Load required packages
##D library(ggplot2)  # for ggtitle() function
##D library(nnet)     # for fitting neural networks
##D 
##D # Simulate training data
##D trn <- gen_friedman(500, seed = 101)  # ?vip::gen_friedman
##D 
##D # Inspect data
##D tibble::as.tibble(trn)
##D 
##D # Fit PPR and NN models (hyperparameters were chosen using the caret package
##D # with 5 repeats of 5-fold cross-validation)
##D pp <- ppr(y ~ ., data = trn, nterms = 11)
##D set.seed(0803) # for reproducibility
##D nn <- nnet(y ~ ., data = trn, size = 7, decay = 0.1, linout = TRUE,
##D            maxit = 500)
##D 
##D # Plot VI scores
##D set.seed(2021)  # for reproducibility
##D p1 <- vip(pp, method = "permute", target = "y", metric = "rsquared",
##D           pred_wrapper = predict) + ggtitle("PPR")
##D p2 <- vip(nn, method = "permute", target = "y", metric = "rsquared",
##D           pred_wrapper = predict) + ggtitle("NN")
##D grid.arrange(p1, p2, ncol = 2)
##D 
##D # Mean absolute error
##D mae <- function(actual, predicted) {
##D   mean(abs(actual - predicted))
##D }
##D 
##D # Permutation-based VIP with user-defined MAE metric
##D set.seed(1101)  # for reproducibility
##D vip(pp, method = "permute", target = "y", metric = mae,
##D     smaller_is_better = TRUE,
##D     pred_wrapper = function(object, newdata) predict(object, newdata)
##D ) + ggtitle("PPR")
## End(Not run)



cleanEx()
nameEx("vint")
### * vint

flush(stderr()); flush(stdout())

### Name: vint
### Title: Interaction effects
### Aliases: vint

### ** Examples

## Not run: 
##D #
##D # The Friedman 1 benchmark problem
##D #
##D 
##D # Load required packages
##D library(gbm)
##D library(ggplot2)
##D library(mlbench)
##D 
##D # Simulate training data
##D trn <- gen_friedman(500, seed = 101)  # ?vip::gen_friedman
##D 
##D #
##D # NOTE: The only interaction that actually occurs in the model from which
##D # these data are generated is between x.1 and x.2!
##D #
##D 
##D # Fit a GBM to the training data
##D set.seed(102)  # for reproducibility
##D fit <- gbm(y ~ ., data = trn, distribution = "gaussian", n.trees = 1000,
##D            interaction.depth = 2, shrinkage = 0.01, bag.fraction = 0.8,
##D            cv.folds = 5)
##D best_iter <- gbm.perf(fit, plot.it = FALSE, method = "cv")
##D 
##D # Quantify relative interaction strength
##D all_pairs <- combn(paste0("x.", 1:10), m = 2)
##D res <- NULL
##D for (i in seq_along(all_pairs)) {
##D   interact <- vint(fit, feature_names = all_pairs[, i], n.trees = best_iter)
##D   res <- rbind(res, interact)
##D }
##D 
##D # Plot top 20 results
##D top_20 <- res[1L:20L, ]
##D ggplot(top_20, aes(x = reorder(Variables, Interaction), y = Interaction)) +
##D   geom_col() +
##D   coord_flip() +
##D   xlab("") +
##D   ylab("Interaction strength")
## End(Not run)



cleanEx()
nameEx("vip")
### * vip

flush(stderr()); flush(stdout())

### Name: vip
### Title: Variable importance plots
### Aliases: vip vip.default vip.model_fit

### ** Examples

#
# A projection pursuit regression example
#

# Load the sample data
data(mtcars)

# Fit a projection pursuit regression model
model <- ppr(mpg ~ ., data = mtcars, nterms = 1)

# Construct variable importance plot
vip(model, method = "firm")

# Better yet, store the variable importance scores and then plot
vi_scores <- vi(model, method = "firm")
vip(vi_scores, geom = "point", horiz = FALSE)
vip(vi_scores, geom = "point", horiz = FALSE, aesthetics = list(size = 3))

# The `%T>%` operator is imported for convenience; see ?magrittr::`%T>%`
# for details
vi_scores <- model %>%
  vi(method = "firm") %T>%
  {print(vip(.))}
vi_scores

# Permutation scores (barplot w/ raw values and jittering)
pfun <- function(object, newdata) predict(object, newdata = newdata)
vip(model, method = "permute", train = mtcars, target = "mpg", nsim = 10,
    metric = "rmse", pred_wrapper = pfun,
    aesthetics = list(color = "grey50", fill = "grey50"),
    all_permutations = TRUE, jitter = TRUE)

# Permutation scores (boxplot)
vip(model, method = "permute", train = mtcars, target = "mpg", nsim = 10,
    metric = "rmse", pred_wrapper = pfun, geom = "boxplot")

# Permutation scores (boxplot colored by feature)
library(ggplot2)  # for `aes_string()` function
vip(model, method = "permute", train = mtcars, target = "mpg", nsim = 10,
    metric = "rmse", pred_wrapper = pfun, geom = "boxplot",
    all_permutations = TRUE, mapping = aes_string(fill = "Variable"),
    aesthetics = list(color = "grey35", size = 0.8))



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
