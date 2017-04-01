#' GBM performance
#' 
#' Estimates the optimal number of boosting iterations for a \code{gbm} object
#' and optionally plots various performance measures
#' 
#' @param object A \code{\link{gbm.object}} created from an initial call to
#' \code{\link{gbm}}.
#'   
#' @param plot.it An indicator of whether or not to plot the performance
#' measures. Setting \code{plot.it = TRUE} creates two plots. The first plot
#' plots \code{object$train.error} (in black) and \code{object$valid.error} 
#' (in red) versus the iteration number. The scale of the error measurement, 
#' shown on the left vertical axis, depends on the \code{distribution} 
#' argument used in the initial call to \code{\link{gbm}}.
#'   
#' @param oobag.curve Indicates whether to plot the out-of-bag performance
#' measures in a second plot.
#'   
#' @param overlay If TRUE and oobag.curve=TRUE then a right y-axis is added to
#' the training and test error plot and the estimated cumulative improvement 
#' in the loss function is plotted versus the iteration number.
#'   
#' @param method Indicate the method used to estimate the optimal number of
#' boosting iterations. \code{method = "OOB"} computes the out-of-bag estimate
#' and \code{method = "test"} uses the test (or validation) dataset to compute 
#' an out-of-sample estimate. \code{method = "cv"} extracts the optimal number 
#' of iterations using cross-validation if \code{gbm} was called with
#' \code{cv.folds} > 1.
#'   
#' @return \code{gbm.perf} Returns the estimated optimal number of iterations.
#'   The method of computation depends on the \code{method} argument.
#'   
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' 
#' @seealso \code{\link{gbm}}, \code{\link{gbm.object}}
#' 
#' @keywords nonlinear survival nonparametric tree
#' 
#' @export
gbm.perf <- function(object, plot.it = TRUE, oobag.curve = FALSE,
                     overlay = TRUE, method) {
  
  # Determine method, if missing
  if (missing(method)) {
    method <- guess_error_method(object)
  }
  
  # Determine "optimal" number of iterations
  best.iter <- best_iter(object, method = method)
  
  # Determine an appropriate y-axis label
  ylab <- get_ylab(object)
  
  # Determine an appropriate range for the y-axis
  ylim <- get_ylim(object, method = method)
  
  # Plot results
  plot(object$train.error, ylim = ylim, type = "l", xlab = "Iteration", 
       ylab = ylab)
  
  if(object$train.fraction!=1) {
    lines(object$valid.error,col="red")
  }
  if(method=="cv") {
    lines(object$cv.error,col="green")
  }
  if(!is.na(best.iter)) {
    abline(v=best.iter,col="blue",lwd=2,lty=2)
  }
  if(oobag.curve) {
    if(overlay) {
      smoother <- attr(best.iter, "smoother")
      par(new = TRUE)
      plot(smoother$x,
           cumsum(smoother$y),
           col="blue",
           type="l",
           xlab="",ylab="",
           axes=FALSE)
      axis(4,srt=0)
      at <- mean(range(smoother$y))
      mtext(paste("OOB improvement in",ylab),side=4,srt=270,line=2)
      abline(h=0,col="blue",lwd=2)
    }
    
    plot(object$oobag.improve,type="l",
         xlab="Iteration",
         ylab=paste("OOB change in",ylab))
    lines(smoother,col="red",lwd=2)
    abline(h=0,col="blue",lwd=1)
    
    abline(v=best.iter,col="blue",lwd=1)
  }
  return(best.iter)
}
