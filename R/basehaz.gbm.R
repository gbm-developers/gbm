# rd2rox <- function(path = file.choose()) {
#   info <- Rd2roxygen::parse_file(path)
#   cat(Rd2roxygen::create_roxygen(info), sep = "\n")
# }


#' Baseline hazard function
#' 
#' Computes the Breslow estimator of the baseline hazard function for a
#' proportional hazard regression model.
#' 
#' The proportional hazard model assumes h(t|x)=lambda(t)*exp(f(x)).
#' \code{\link{gbm}} can estimate the f(x) component via partial likelihood.
#' After estimating f(x), \code{basehaz.gbm} can compute the a nonparametric
#' estimate of lambda(t).
#' 
#' @param t The survival times.
#' @param delta The censoring indicator.
#' @param f.x The predicted values of the regression model on the log hazard
#'   scale.
#' @param t.eval Values at which the baseline hazard will be evaluated.
#' @param smooth If \code{TRUE} \code{basehaz.gbm} will smooth the estimated
#'   baseline hazard using Friedman's super smoother \code{\link{supsmu}}.
#' @param cumulative If \code{TRUE} the cumulative survival function will be
#'   computed.
#' @return A vector of length equal to the length of t (or of length
#'   \code{t.eval} if \code{t.eval} is not \code{NULL}) containing the baseline
#'   hazard evaluated at t (or at \code{t.eval} if \code{t.eval} is not
#'  \code{NULL}). If \code{cumulative} is set to \code{TRUE} then the returned
#'   vector evaluates the cumulative hazard function at those values.
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' @seealso \code{\link[survival]{survfit}}, \code{\link{gbm}}
#' @references 
#' N. Breslow (1972). "Discussion of `Regression Models and
#' Life-Tables' by D.R. Cox," Journal of the Royal Statistical Society, Series
#' B, 34(2):216-217.
#' 
#' N. Breslow (1974). "Covariance analysis of censored survival data,"
#' Biometrics 30:89-99.
#' @keywords methods survival
#' @export
basehaz.gbm <- function(t,delta, f.x, t.eval = NULL, smooth = FALSE,
                        cumulative = TRUE) {
  
  t.unique <- sort(unique(t[delta==1]))
  alpha <- length(t.unique)
  for(i in 1:length(t.unique)) {
    alpha[i] <- sum(t[delta==1]==t.unique[i])/
      sum(exp(f.x[t>=t.unique[i]]))
  }
  
  if(!smooth && !cumulative) {
    if(!is.null(t.eval)) {
      stop("Cannot evaluate unsmoothed baseline hazard at t.eval.")
    }
  } else {
    if(smooth && !cumulative) {
      lambda.smooth <- supsmu(t.unique,alpha)
    } else {
      if(smooth && cumulative)
      {
        lambda.smooth <- supsmu(t.unique, cumsum(alpha))
      } else {  # (!smooth && cumulative) - THE DEFAULT
        lambda.smooth <- list(x = t.unique, y = cumsum(alpha))
      }
    }
  }

  
  obj <- if(!is.null(t.eval)) {
    approx(lambda.smooth$x, lambda.smooth$y, xout = t.eval)$y
  } else {
    approx(lambda.smooth$x, lambda.smooth$y, xout = t)$y
  }
  
  return(obj)
  
}
