#' Arrange multiple grobs on a page
#'
#' See \code{\link[gridExtra]{grid.arrange}} for more details.
#'
#' @name grid.arrange
#' @rdname grid.arrange
#' @keywords internal
#' @export
#' @importFrom gridExtra grid.arrange
#' @usage grid.arrange(..., newpage = TRUE)
NULL


#' @keywords internal
getAvailableDistributions <- function() {
  c("adaboost", "bernoulli", "coxph", "gaussian", "huberized", "laplace", 
    "multinomial", "pairwise", "poisson", "quantile", "tdist")
}


#' @keywords internal
guess_error_method <- function(object) {
  if (has_train_test_split(object)) {
    "test"
  } else if (has_cross_validation(object)) {
    "cv"
  } else {
    "OOB"
  }
}


#' @keywords internal
has_train_test_split <- function(object) {
  object$train.fraction < 1
}


#' @keywords internal
has_cross_validation <- function(object) {
  !is.null(object$cv.error)
}


#' @keywords internal
best_iter <- function(object, method) {
  check_if_gbm_fit(object)
  if (method == "OOB") {
    best_iter_out_of_bag(object)
  } else if (method == "test") {
    best_iter_test(object)
  } else if (method == "cv") {
    best_iter_cv(object)
  } else {
    stop("method must be one of \"cv\", \"test\", or \"OOB\"")
  }
}


#' @keywords internal
best_iter_test <- function(object) {
  check_if_gbm_fit(object)
  best_iter_test <- which.min(object$valid.error)
  return(best_iter_test)
}


#' @keywords internal
best_iter_cv <- function(object) {
  check_if_gbm_fit(object)
  if(!has_cross_validation(object)) {
    stop('In order to use method="cv" gbm must be called with cv_folds>1.')
  }
  best_iter_cv <- which.min(object$cv.error)
  return(best_iter_cv)
}


#' @keywords internal
best_iter_out_of_bag <- function(object) {
  check_if_gbm_fit(object)
  if(object$bag.fraction == 1) {
    stop("Cannot compute OOB estimate or the OOB curve when bag_fraction=1.")
  }
  if(all(!is.finite(object$oobag.improve))) {
    stop("Cannot compute OOB estimate or the OOB curve. No finite OOB ", 
         "estimates of improvement.")
  }
  message("OOB generally underestimates the optimal number of iterations ", 
          "although predictive performance is reasonably competitive. Using ", 
          "cv_folds>1 when calling gbm usually results in improved predictive ",
          "performance.")
  smoother <- generate_smoother_oobag(object)
  best_iter_oob <- smoother$x[which.min(-cumsum(smoother$y))]
  attr(best_iter_oob, "smoother") <- smoother
  return(best_iter_oob)
}


#' @keywords internal
generate_smoother_oobag <- function(object) {
  check_if_gbm_fit(object)
  x <- seq_len(object$n.trees)
  smoother <- loess(object$oobag.improve ~ x,
                    enp.target = min(max(4, length(x) / 10), 50))
  smoother$y <- smoother$fitted
  smoother$x <- x
  return(smoother)
}


#' @keywords internal
check_if_gbm_fit <- function(object) {
  if (!inherits(object, "gbm")) {
    stop(deparse(substitute(object)), " is not a valid \"gbm\" object.")
  }
}


#' @keywords internal
get_ylab <- function(object) {
  check_if_gbm_fit(object)
  if (object$distribution$name != "pairwise") {
    switch(substring(object$distribution$name, 1, 2), 
           ga = "Squared error loss", 
           be = "Bernoulli deviance", 
           po = "Poisson deviance", 
           ad = "AdaBoost exponential bound", 
           co = "Cox partial deviance", 
           la = "Absolute loss", 
           qu = "Quantile loss", 
           mu = "Multinomial deviance", 
           td = "t-distribution deviance")
  } else {
    switch(object$distribution$metric, 
           conc = "Fraction of concordant pairs", 
           ndcg = "Normalized discounted cumulative gain", 
           map = "Mean average precision",
           mrr = "Mean reciprocal rank")
  }
}


#' @keywords internal
get_ylim <- function(object, method) {
  check_if_gbm_fit(object)
  if(object$train.fraction == 1) {
    if ( method=="cv" ) { 
      range(object$train.error, object$cv.error) 
    } else if ( method == "test" ) { 
      range( object$train.error, object$valid.error) 
    } else { 
      range(object$train.error) 
    }
  } else {
    range(object$train.error, object$valid.error)
  }
}
