#' gbm internal functions
#' 
#' Helper functions for preprocessing data prior to building a \code{"gbm"}
#' object.
#'
#' @param y The response variable.
#' @param d,distribution The distribution, either specified by the user or
#' implied.
#' @param class.stratify.cv Whether or not to stratify, if provided by the user.
#' @param i.train Computed internally by \code{gbm}.
#' @param group The group, if using \code{distibution = "pairwise"}.
#' @param strat Whether or not to stratify.
#' @param cv.folds The number of cross-validation folds.
#' @param x The design matrix.
#' @param id The interaction depth.
#' @param w The weights.
#' @param n The number of cores to use in the cluster.
#' @param o The offset.
#' 
#' @details 
#' These are functions used internally by \code{gbm} and not intended for direct 
#' use by the user.
#' 
#' @aliases guessDist getStratify getCVgroup checkMissing checkID checkWeights
#' checkOffset getVarNames gbmCluster
#' 
#' @rdname gbm-internals
#' @export
guessDist <- function(y){
  # If distribution is not given, try to guess it
  if (length(unique(y)) == 2){ d <- "bernoulli" }
  else if (class(y) == "Surv" ){ d <- "coxph" }
  else if (is.factor(y)){ d <- "multinomial" }
  else{ d <- "gaussian" }
  cat(paste("Distribution not specified, assuming", d, "...\n"))
  list(name=d)
}


#' @rdname gbm-internals
#' @export
getCVgroup <- function(distribution, class.stratify.cv, y, i.train, cv.folds, 
                       group) {
    # Construct cross-validation groups depending on the type of model to be fit
    if (distribution$name %in% c( "bernoulli", "multinomial" ) & class.stratify.cv ){
      nc <- table(y[i.train]) # Number in each class
      uc <- names(nc)
      if (min(nc) < cv.folds){
        stop( paste("The smallest class has only", min(nc), "objects in the training set. Can't do", cv.folds, "fold cross-validation."))
      }
      cv.group <- vector(length = length(i.train))
      for (i in 1:length(uc)){
        cv.group[y[i.train] == uc[i]] <- sample(rep(1:cv.folds , length = nc[i]))
      }
    } # Close if
    else if (distribution$name == "pairwise") {
      # Split into CV folds at group boundaries
      s <- sample(rep(1:cv.folds, length=nlevels(group)))
      cv.group <- s[as.integer(group[i.train])]
    }
    else {
      cv.group <- sample(rep(1:cv.folds, length=length(i.train)))
    }
    cv.group
  }


#' @rdname gbm-internals
#' @export
getStratify <- function(strat, d){
  if (is.null(strat)){
    if (d$name == "multinomial" ){ strat <- TRUE }
    else { strat <- FALSE }
  }
  else {
    if (!is.element(d$name, c( "bernoulli", "multinomial"))){
      warning("You can only use class.stratify.cv when distribution is bernoulli or multinomial. Ignored.")
      strat <- FALSE
    }
  } # Close else
  strat
}


#' @rdname gbm-internals
#' @export
checkMissing <- function(x, y){
  nms <- getVarNames(x)
  #### Check for NaNs in x and NAs in response
  j <- apply(x, 2, function(z) any(is.nan(z)))
  if(any(j)) {
    stop("Use NA for missing values. NaN found in predictor variables:",
         paste(nms[j],collapse=","))
  }
  if(any(is.na(y))) stop("Missing values are not allowed in the response")
  invisible(NULL)
}


#' @rdname gbm-internals
#' @export
checkWeights <- function(w, n){
  # Logical checks on weights
  if(length(w)==0) { w <- rep(1, n) }
  else if(any(w < 0)) stop("negative weights not allowed")
  w
}


#' @rdname gbm-internals
#' @export
checkID <- function(id){
  # Check for disallowed interaction.depth
  if(id < 1) {
    stop("interaction.depth must be at least 1.")
  }
  else if(id > 49) {
    stop("interaction.depth must be less than 50. You should also ask yourself why you want such large interaction terms. A value between 1 and 5 should be sufficient for most applications.")
  }
  invisible(id)
}


#' @rdname gbm-internals
#' @export
checkOffset <- function(o, y){
  # Check offset
  if(is.null(o) | all(o==0)) { o <- NA  }
  else if(length(o) != length(y))   {
    stop("The length of offset does not equal the length of y.")
  }
  o
}


#' @rdname gbm-internals
#' @export
getVarNames <- function(x){
  if(is.matrix(x)) { var.names <- colnames(x) }
  else if(is.data.frame(x)) { var.names <- names(x) }
  else { var.names <- paste("X",1:ncol(x),sep="") }
  var.names
}


#' @rdname gbm-internals
#' @export
gbmCluster <- function(n) {
  # If number of cores (n) not given, try to work it out from the number
  # that appear to be available and the number of CV folds.
  if (is.null(n)) {
    n <- parallel::detectCores()
  }
  parallel::makeCluster(n)
}
