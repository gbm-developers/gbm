#' Cross-validate a gbm
#' 
#' Functions for cross-validating gbm. These functions are used internally and
#' are not intended for end-user direct usage.
#' 
#' These functions are not intended for end-user direct usage, but are used
#' internally by \code{gbm}.
#' 
#' @aliases gbmCrossVal gbmCrossValModelBuild gbmDoFold gbmCrossValErr
#' gbmCrossValPredictions
#' @param cv.folds The number of cross-validation folds.
#' @param nTrain The number of training samples.
#' @param n.cores The number of cores to use.
#' @param class.stratify.cv Whether or not stratified cross-validation samples
#' are used.
#' @param data The data.
#' @param x The model matrix.
#' @param y The response variable.
#' @param offset The offset.
#' @param distribution The type of loss function. See \code{\link{gbm}}.
#' @param w Observation weights.
#' @param var.monotone See \code{\link{gbm}}.
#' @param n.trees The number of trees to fit.
#' @param interaction.depth The degree of allowed interactions. See
#' \code{\link{gbm}}.
#' @param n.minobsinnode See \code{\link{gbm}}.
#' @param shrinkage See \code{\link{gbm}}.
#' @param bag.fraction See \code{\link{gbm}}.
#' @param var.names See \code{\link{gbm}}.
#' @param response.name See \code{\link{gbm}}.
#' @param group Used when \code{distribution = "pairwise"}. See
#' \code{\link{gbm}}.
#' @param i.train Items in the training set.
#' @param cv.models A list containing the models for each fold.
#' @param cv.group A vector indicating the cross-validation fold for each
#' member of the training set.
#' @param best.iter.cv The iteration with lowest cross-validation error.
#' @param X Index (cross-validation fold) on which to subset.
#' @param s Random seed.
#' @return A list containing the cross-validation error and predictions.
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @seealso \code{\link{gbm}}
#' @references J.H. Friedman (2001). "Greedy Function Approximation: A Gradient
#' Boosting Machine," Annals of Statistics 29(5):1189-1232.
#' 
#' L. Breiman (2001).
#' \url{https://www.stat.berkeley.edu/users/breiman/randomforest2001.pdf}.
#' @keywords models

# Perform gbm cross-validation
#
# This function has far too many arguments, but there isn't the
# abstraction in gbm to lose them.
#' @rdname gbmCrossVal
#' @export
gbmCrossVal <- function(cv.folds, nTrain, n.cores,
                        class.stratify.cv, data,
                        x, y, offset, distribution, w, var.monotone,
                        n.trees, interaction.depth, n.minobsinnode,
                        shrinkage, bag.fraction,
                        var.names, response.name, group) {
  i.train <- 1:nTrain
  cv.group <- getCVgroup(distribution, class.stratify.cv, y,
                         i.train, cv.folds, group)
  ## build the models
  cv.models <- gbmCrossValModelBuild(cv.folds, cv.group, n.cores,
                                     i.train, x, y, offset,
                                     distribution, w, var.monotone,
                                     n.trees, interaction.depth,
                                     n.minobsinnode, shrinkage,
                                     bag.fraction, var.names,
                                     response.name, group)

  ## get the errors
  cv.error  <- gbmCrossValErr(cv.models, cv.folds, cv.group, nTrain, n.trees)
  best.iter.cv <- which.min(cv.error)

  ## get the predictions
  predictions <- gbmCrossValPredictions(cv.models, cv.folds, cv.group,
                                        best.iter.cv, distribution,
                                        data[i.train, ], y)
  list(error = cv.error, predictions = predictions)
}


# Get the gbm cross-validation error
#' @rdname gbmCrossVal
#' @export
gbmCrossValErr <- function(cv.models, cv.folds, cv.group, nTrain, n.trees) {
  in.group <- tabulate(cv.group, nbins = cv.folds)
  cv.error <- vapply(1:cv.folds,
                     function(index) {
                       model <- cv.models[[index]]
                       model$valid.error * in.group[[index]]
                     }, double(n.trees))
  ## this is now a (n.trees, cv.folds) matrix
  
  ## and now a n.trees vector
  rowSums(cv.error) / nTrain
}


#' @rdname gbmCrossVal
#' @export
gbmCrossValPredictions <- function(cv.models, cv.folds, cv.group,
                                   best.iter.cv, distribution, data, y) {
  
  # Get the predictions for GBM cross validation. This function is not as nice 
  # as it could be (i.e., leakage of y)
  
  # Test that cv.group and data match
  if (nrow(data) != length(cv.group)) {
    stop("Mismatch between `data` and `cv.group`.")
  }
  
  # This is a little complicated due to multinomial distribution
  num.cols <- if (distribution$name == "multinomial") {
    nlevels(factor(y))
  } else {
    1
  }
  
  # Initialize results matrix
  res <- matrix(nrow = nrow(data), ncol = num.cols)
  
  # There's no real reason to do this as other than a for loop
  data.names <- names(data)  # column names
  for (ind in 1:cv.folds) {
    
    # These are the particular elements
    flag <- cv.group == ind
    model <- cv.models[[ind]]
    
    # The %in% here is to handle coxph
    # my.data  <- data[flag, !(data.names %in% model$response.name)]
    my.data  <- data[flag, model$var.names]
    predictions <- predict(model, newdata = my.data, n.trees = best.iter.cv)  # FIXME
    predictions <- matrix(predictions, ncol = num.cols)
    res[flag, ] <- predictions
  }
  
  # Handle multinomial case
  if (distribution$name != "multinomial") {
    res <- as.numeric(res)
  }
  
  # Return the result
  res
  
}


# Perform gbm cross-validation
#
# This function has far too many arguments.
#' @rdname gbmCrossVal
#' @export
gbmCrossValModelBuild <- function(cv.folds, cv.group, n.cores, i.train, x, y, 
                                  offset, distribution, w, var.monotone, 
                                  n.trees, interaction.depth, n.minobsinnode,
                                  shrinkage, bag.fraction, var.names, 
                                  response.name, group) {
  
  # Set random seeds
  seeds <- as.integer(runif(cv.folds, -(2^31 - 1), 2^31))

  # Perform cross-validation model builds
  if (!is.null(n.cores) && n.cores == 1) {  
    lapply(1:cv.folds, FUN = gbmDoFold, i.train, x, 
           y, offset, distribution, w, var.monotone, n.trees,
           interaction.depth, n.minobsinnode, shrinkage, 
           bag.fraction, cv.group, var.names, response.name, group, 
           seeds)
  } else {
    # Set up cluster and add finalizer
    cluster <- gbmCluster(n.cores)
    on.exit(parallel::stopCluster(cluster))
    parallel::parLapply(cl = cluster, X = 1:cv.folds, fun = gbmDoFold, i.train, x, 
                        y, offset, distribution, w, var.monotone, n.trees,
                        interaction.depth, n.minobsinnode, shrinkage, 
                        bag.fraction, cv.group, var.names, response.name, group, 
                        seeds)
  }

}


#' @rdname gbmCrossVal
#' @export
gbmDoFold <- function(X, i.train, x, y, offset, distribution, w, var.monotone, 
                      n.trees, interaction.depth, n.minobsinnode, shrinkage, 
                      bag.fraction, cv.group, var.names, response.name, group, 
                      s) {
  
  # Do specified cross-validation fold---a self-contained function for passing 
  # to individual cores.
  
  # Load required packages for core
  library(gbm, quietly = TRUE)
  
  # Print CV information
  cat("CV:", X, "\n")
  
  # Setup
  set.seed(s[[X]])
  i <- order(cv.group == X)
  x <- x[i.train, , drop = TRUE][i, , drop = FALSE]
  y <- y[i.train][i]
  offset <- offset[i.train][i]
  nTrain <- length(which(cv.group != X))
  group <- group[i.train][i]
  
  # Return a fitted GBM
  gbm.fit(x = x, y = y, offset = offset, distribution = distribution,
          w = w, var.monotone = var.monotone, n.trees = n.trees,
          interaction.depth = interaction.depth,
          n.minobsinnode = n.minobsinnode,
          shrinkage = shrinkage, bag.fraction = bag.fraction,
          nTrain = nTrain, keep.data = FALSE, verbose = FALSE, 
          response.name = response.name, group = group)
  
}
