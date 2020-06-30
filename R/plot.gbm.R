#' Marginal plots of fitted gbm objects
#' 
#' Plots the marginal effect of the selected variables by "integrating" out the
#' other variables.
#' 
#' \code{plot.gbm} produces low dimensional projections of the
#' \code{\link{gbm.object}} by integrating out the variables not included in
#' the \code{i.var} argument. The function selects a grid of points and uses
#' the weighted tree traversal method described in Friedman (2001) to do the
#' integration. Based on the variable types included in the projection,
#' \code{plot.gbm} selects an appropriate display choosing amongst line plots,
#' contour plots, and \code{\link[lattice:Lattice]{lattice}} plots. If the default
#' graphics are not sufficient the user may set \code{return.grid = TRUE}, store
#' the result of the function, and develop another graphic display more
#' appropriate to the particular example.
#' 
#' @param x A \code{\link{gbm.object}} that was fit using a call to 
#' \code{\link{gbm}}.
#' 
#' @param i.var Vector of indices or the names of the variables to plot. If
#' using indices, the variables are indexed in the same order that they appear
#' in the initial \code{gbm} formula. If \code{length(i.var)} is between 1 and
#' 3 then \code{plot.gbm} produces the plots. Otherwise, \code{plot.gbm}
#' returns only the grid of evaluation points and their average predictions
#' 
#' @param n.trees Integer specifying the number of trees to use to generate the 
#' plot. Default is to use \code{x$n.trees} (i.e., the entire ensemble).
#' 
#' @param continuous.resolution Integer specifying the number of equally space 
#' points at which to evaluate continuous predictors.
#' 
#' @param return.grid Logical indicating whether or not to produce graphics 
#' \code{FALSE} or only return the grid of evaluation points and their average
#' predictions \code{TRUE}. This is useful for customizing the graphics for 
#' special variable types, or for higher dimensional graphs.
#' 
#' @param type Character string specifying the type of prediction to plot on the 
#' vertical axis. See \code{\link{predict.gbm}} for details.
#' 
#' @param level.plot Logical indicating whether or not to use a false color 
#' level plot (\code{TRUE}) or a 3-D surface (\code{FALSE}). Default is 
#' \code{TRUE}.
#'
#' @param contour Logical indicating whether or not to add contour lines to the
#' level plot. Only used when \code{level.plot = TRUE}. Default is \code{FALSE}.
#'
#' @param number Integer specifying the number of conditional intervals to use
#' for the continuous panel variables. See \code{\link[graphics:coplot]{co.intervals}}
#' and \code{\link[lattice:shingles]{equal.count}} for further details.
#'
#' @param overlap The fraction of overlap of the conditioning variables. See
#' \code{\link[graphics:coplot]{co.intervals}} and \code{\link[lattice:shingles]{equal.count}}
#' for further details.
#'
#' @param col.regions Color vector to be used if \code{level.plot} is
#' \code{TRUE}. Defaults to the wonderful Matplotlib 'viridis' color map
#' provided by the \code{viridis} package. See \code{\link[viridis:reexports]{viridis}}
#' for details.
#' 
#' @param ... Additional optional arguments to be passed onto 
#' \code{\link[graphics:plot.default]{plot}}.
#' 
#' @return If \code{return.grid = TRUE}, a grid of evaluation points and their 
#' average predictions. Otherwise, a plot is returned.
#' 
#' @note More flexible plotting is available using the 
#' \code{\link[pdp]{partial}} and \code{\link[pdp]{plotPartial}} functions.
#' 
#' @seealso \code{\link[pdp]{partial}}, \code{\link[pdp]{plotPartial}}, 
#' \code{\link{gbm}}, and \code{\link{gbm.object}}.
#'
#' @references J. H. Friedman (2001). "Greedy Function Approximation: A Gradient
#' Boosting Machine," Annals of Statistics 29(4).
#' 
#' @references B. M. Greenwell (2017). "pdp: An R Package for Constructing 
#' Partial Dependence Plots," The R Journal 9(1), 421--436. 
#' \url{https://journal.r-project.org/archive/2017/RJ-2017-016/index.html}.
#' 
#' @export plot.gbm
#' @export
plot.gbm <- function(x, i.var = 1, n.trees = x$n.trees, 
                     continuous.resolution = 100, return.grid = FALSE, 
                     type = c("link", "response"), level.plot = TRUE, 
                     contour = FALSE, number = 4, overlap = 0.1,
                     col.regions = viridis::viridis, ...) {
  
  # Match type argument
  type <- match.arg(type)
  
  # Sanity checks
  if(all(is.character(i.var))) {
    i <- match(i.var, x$var.names)
    if(any(is.na(i))) {
      stop("Requested variables not found in ", deparse(substitute(x)), ": ", 
           i.var[is.na(i)])
    } else {
      i.var <- i
    }
  }
  if((min(i.var) < 1) || (max(i.var) > length(x$var.names))) {
    warning("i.var must be between 1 and ", length(x$var.names))
  }
  if(n.trees > x$n.trees) {
    warning(paste("n.trees exceeds the number of tree(s) in the model: ",
                  x$n.trees, ". Using ", x$n.trees, 
                  " tree(s) instead.", sep = ""))
    n.trees <- x$n.trees
  }
  
  if(length(i.var) > 3) {
    warning("plot.gbm() will only create up to (and including) 3-way ", 
            "interaction plots.\nBeyond that, plot.gbm() will only return ",
            "the plotting data structure.")
    return.grid <- TRUE
  }
  
  # Generate grid of predictor values on which to compute the partial 
  # dependence values
  grid.levels <- vector("list", length(i.var))
  for(i in 1:length(i.var)) {
    if(is.numeric(x$var.levels[[i.var[i]]])) {  # continuous
      grid.levels[[i]] <- seq(from = min(x$var.levels[[i.var[i]]]), 
                              to = max(x$var.levels[[i.var[i]]]),
                              length = continuous.resolution)
    } else {  # categorical
      grid.levels[[i]] <- 
        as.numeric(factor(x$var.levels[[i.var[i]]], 
                          levels = x$var.levels[[i.var[i]]])) - 1
    }
  }
  X <- expand.grid(grid.levels)
  names(X) <- paste("X", 1:length(i.var), sep = "")
  
  # For compatibility with gbm version 1.6
  if (is.null(x$num.classes)) {
    x$num.classes <- 1
  }
  
  # Compute partial dependence values
  y <- .Call("gbm_plot", X = as.double(data.matrix(X)), 
             cRows = as.integer(nrow(X)), cCols = as.integer(ncol(X)),
             n.class = as.integer(x$num.classes), 
             i.var = as.integer(i.var - 1), n.trees = as.integer(n.trees),
             initF = as.double(x$initF), trees = x$trees, 
             c.splits = x$c.splits, var.type = as.integer(x$var.type),
             PACKAGE = "gbm")
  
  if (x$distribution$name == "multinomial") {  # reshape into matrix
    X$y <- matrix(y, ncol = x$num.classes)
    colnames(X$y) <- x$classes
    
    # Convert to class probabilities (if requested)
    if (type == "response") {
      X$y <- exp(X$y)
      X$y <- X$y / matrix(rowSums(X$y), ncol = ncol(X$y), nrow = nrow(X$y))
    }
  } else if(is.element(x$distribution$name, c("bernoulli", "pairwise")) && 
            type == "response") {
    X$y <- 1 / (1 + exp(-y))
  } else if ((x$distribution$name == "poisson") && (type == "response")) {
    X$y <- exp(y)
  } else if (type == "response"){
    warning("`type = \"response\"` only implemented for \"bernoulli\", ",
            "\"poisson\", \"multinomial\", and \"pairwise\" distributions. ",
            "Ignoring." )
  } else { 
    X$y <- y 
  }
  
  # Transform categorical variables back to factors
  f.factor <- rep(FALSE, length(i.var))
  for(i in 1:length(i.var)) {
    if(!is.numeric(x$var.levels[[i.var[i]]])) {
      X[,i] <- factor(x$var.levels[[i.var[i]]][X[, i] + 1],
                      levels = x$var.levels[[i.var[i]]])
      f.factor[i] <- TRUE
    }
  }
  
  # Return original variable names
  names(X)[1:length(i.var)] <- x$var.names[i.var]
  
  # Return grid only (if requested)
  if(return.grid) {
    return(X)
  }
  
  # Determine number of predictors
  nx <- length(i.var)
  
  # Determine which type of plot to draw based on the number of predictors
  if (nx == 1L) {
    
    # Single predictor
    plotOnePredictorPDP(X, ...)
    
  } else if (nx == 2) {
    
    # Two predictors
    plotTwoPredictorPDP(X, level.plot = level.plot, contour = contour,
                        col.regions = col.regions, ...)
    
  } else {
    
    # Three predictors (paneled version of plotTwoPredictorPDP)
    plotThreePredictorPDP(X, nx = nx, level.plot = level.plot,
                          contour = contour, col.regions = col.regions,
                          number = number, overlap = overlap, ...)
    
  }
  
}


#' @keywords internal
plotOnePredictorPDP <- function(X, ...) {
  
  # Use the first column to determine which type of plot to construct
  if (is.numeric(X[[1L]])) {
    
    # Draw a line plot
    lattice::xyplot(stats::as.formula(paste("y ~", names(X)[1L])), 
                    data = X, type = "l", ...)
    
  } else {
    
    # Draw a Cleveland dot plot
    lattice::dotplot(stats::as.formula(paste("y ~", names(X)[1L])), 
                     data = X, xlab = names(X)[1L], ...)
    
  }
}


#' @keywords internal
plotTwoPredictorPDP <- function(X, level.plot, contour, col.regions, ...) {
  
  # Use the first two columns to determine which type of plot to construct
  if (is.factor(X[[1L]]) && is.factor(X[[2L]])) {
    
    # Draw a Cleveland dot plot
    lattice::dotplot(stats::as.formula(
      paste("y ~", paste(names(X)[1L:2L], collapse = "|"))
    ), data = X, xlab = names(X)[1L], ...)
    
  } else if (is.factor(X[[1L]]) || is.factor(X[[2L]])) {
    
    # Lattice plot formula
    form <- if (is.factor(X[[1L]])) {
      stats::as.formula(paste("y ~", paste(names(X)[2L:1L], collapse = "|")))
    } else {
      stats::as.formula(paste("y ~", paste(names(X)[1L:2L], collapse = "|")))
    }
    
    # Draw a paneled line plot
    lattice::xyplot(form, data = X, type = "l", ...)
    
  } else {
    
    # Lattice plot formula
    form <- stats::as.formula(
      paste("y ~", paste(names(X)[1L:2L], collapse = "*"))
    )
    
    # Draw a three-dimensional surface
    if (level.plot) {
      
      # Draw a false color level plot
      lattice::levelplot(form, data = X, col.regions = col.regions, 
                         contour = contour, ...)
      
    } else {
      
      # Draw a wireframe plot
      lattice::wireframe(form, data = X, ...)
      
    }
    
  }
}


#' @keywords internal
plotThreePredictorPDP <- function(X, nx, level.plot, contour, col.regions, 
                                  number, overlap, ...) {
  
  # Factor, numeric, numeric
  if (is.factor(X[[1L]]) && !is.factor(X[[2L]]) && !is.factor(X[[3L]])) {
    X[, 1L:3L] <- X[, c(2L, 3L, 1L)]
  }
  
  # Numeric, factor, numeric
  if (!is.factor(X[[1L]]) && is.factor(X[[2L]]) && !is.factor(X[[3L]])) {
    X[, 1L:3L] <- X[, c(1L, 3L, 2L)]
  }
  
  # Factor, factor, numeric
  if (is.factor(X[[1L]]) && is.factor(X[[2L]]) && !is.factor(X[[3L]])) {
    X[, 1L:3L] <- X[, c(3L, 1L, 2L)]
  }
  
  # Factor, numeric, factor
  if (is.factor(X[[1L]]) && !is.factor(X[[2L]]) && is.factor(X[[3L]])) {
    X[, 1L:3L] <- X[, c(2L, 1L, 3L)]
  }
  
  # Convert third predictor to a factor using the equal count algorithm
  if (is.numeric(X[[3L]])) {
    X[[3L]] <- equal.count(X[[3L]], number = number, overlap = overlap)
  }
  
  if (is.factor(X[[1L]]) && is.factor(X[[2L]])) {
    
    # Lattice plot formula
    form <- stats::as.formula(
      paste("y ~", names(X)[1L], "|", paste(names(X)[2L:nx], collapse = "*"))
    )
    
    # Produce a paneled dotplot
    lattice::dotplot(form, data = X, xlab = names(X)[1L], ...)
    
  } else if (is.numeric(X[[1L]]) && is.factor(X[[2L]])) {
    
    # Lattice plot formula
    form <- stats::as.formula(
      paste("y ~", names(X)[1L], "|", paste(names(X)[2L:nx], collapse = "*"))
    )
  
    # Produce a paneled lineplot
    lattice::xyplot(form, data = X, type = "l", ...)
    
  } else {
    
    # Lattice plot formula
    form <- stats::as.formula(
      paste("y ~", paste(names(X)[1L:2L], collapse = "*"), "|",
            paste(names(X)[3L:nx], collapse = "*"))
    )
    
    # Draw a three-dimensional surface
    if (level.plot) {
      
      # Draw a false color level plot
      lattice::levelplot(form, data = X, col.regions = col.regions, 
                         contour = contour, ...)
      
    } else {
      
      # Draw a wireframe plot
      lattice::wireframe(form, data = X, ...)
      
    }
    
  }
  
}