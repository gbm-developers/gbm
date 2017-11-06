#' Generalized Boosted Regression Models (GBMs)
#' 
#' This package implements extensions to Freund and Schapire's AdaBoost
#' algorithm and J. Friedman's gradient boosting machine. Includes regression
#' methods for least squares, absolute loss, logistic, Poisson, Cox
#' proportional hazards partial likelihood, multinomial, t-distribution,
#' AdaBoost exponential loss, Learning to Rank, and Huberized hinge loss.
#' 
#' Further information is available in vignette: 
#' \code{browseVignettes(package = "gbm")}
#' 
#' @import lattice
#' 
#' @importFrom grDevices rainbow
#' @importFrom graphics abline axis barplot lines mtext par plot polygon rug 
#' @importFrom graphics segments title
#' @importFrom stats approx binomial delete.response gaussian glm loess 
#' @importFrom stats model.extract model.frame model.offset model.response 
#' @importFrom stats model.weights na.pass poisson predict quantile rbinom 
#' @importFrom stats reformulate reorder rexp rnorm runif sd supsmu terms var
#' @importFrom stats weighted.mean
#' @importFrom survival Surv
#' 
#' @useDynLib gbm, .registration = TRUE
#' 
#' @name gbm-package
#' 
#' @docType package
#' 
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com} with contributions by
#' Daniel Edwards, Brian Kriegler, Stefan Schroedl and Harry Southworth.
#' 
#' @references 
#' Y. Freund and R.E. Schapire (1997) \dQuote{A decision-theoretic
#' generalization of on-line learning and an application to boosting,}
#' \emph{Journal of Computer and System Sciences,} 55(1):119-139.
#' 
#' G. Ridgeway (1999). \dQuote{The state of boosting,} \emph{Computing Science
#' and Statistics} 31:172-181.
#' 
#' J.H. Friedman, T. Hastie, R. Tibshirani (2000). \dQuote{Additive Logistic
#' Regression: a Statistical View of Boosting,} \emph{Annals of Statistics}
#' 28(2):337-374.
#' 
#' J.H. Friedman (2001). \dQuote{Greedy Function Approximation: A Gradient
#' Boosting Machine,} \emph{Annals of Statistics} 29(5):1189-1232.
#' 
#' J.H. Friedman (2002). \dQuote{Stochastic Gradient Boosting,}
#' \emph{Computational Statistics and Data Analysis} 38(4):367-378.
#' 
#' The \url{http://statweb.stanford.edu/~jhf/R-MART} website.
#' 
#' @keywords package
NULL