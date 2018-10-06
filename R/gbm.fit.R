#' Generalized Boosted Regression Modeling (GBM)
#' 
#' Workhorse function providing the link between R and the C++ gbm engine.
#' \code{gbm} is a front-end to \code{gbm.fit} that uses the familiar R
#' modeling formulas. However, \code{\link[stats]{model.frame}} is very slow if
#' there are many predictor variables. For power-users with many variables use
#' \code{gbm.fit}. For general practice \code{gbm} is preferable.
#'
#' @param x A data frame or matrix containing the predictor variables. The 
#' number of rows in \code{x} must be the same as the length of \code{y}. 
#' 
#' @param y A vector of outcomes. The number of rows in \code{x} must be the 
#' same as the length of \code{y}.
#' 
#' @param offset A vector of offset values.
#' 
#' @param misc An R object that is simply passed on to the gbm engine. It can be 
#' used for additional data for the specific distribution. Currently it is only 
#' used for passing the censoring indicator for the Cox proportional hazards 
#' model.
#' 
#' @param distribution Either a character string specifying the name of the
#' distribution to use or a list with a component \code{name} specifying the
#' distribution and any additional parameters needed. If not specified,
#' \code{gbm} will try to guess: if the response has only 2 unique values,
#' bernoulli is assumed; otherwise, if the response is a factor, multinomial is
#' assumed; otherwise, if the response has class \code{"Surv"}, coxph is 
#' assumed; otherwise, gaussian is assumed.
#' 
#' Currently available options are \code{"gaussian"} (squared error), 
#' \code{"laplace"} (absolute loss), \code{"tdist"} (t-distribution loss), 
#' \code{"bernoulli"} (logistic regression for 0-1 outcomes), 
#' \code{"huberized"} (huberized hinge loss for 0-1 outcomes), classes), 
#' \code{"adaboost"} (the AdaBoost exponential loss for 0-1 outcomes),
#' \code{"poisson"} (count outcomes), \code{"coxph"} (right censored 
#' observations), \code{"quantile"}, or \code{"pairwise"} (ranking measure 
#' using the LambdaMart algorithm).
#' 
#' If quantile regression is specified, \code{distribution} must be a list of
#' the form \code{list(name = "quantile", alpha = 0.25)} where \code{alpha} is 
#' the quantile to estimate. The current version's quantile regression method 
#' does not handle non-constant weights and will stop.
#' 
#' If \code{"tdist"} is specified, the default degrees of freedom is 4 and 
#' this can be controlled by specifying 
#' \code{distribution = list(name = "tdist", df = DF)} where \code{DF} is your 
#' chosen degrees of freedom.
#' 
#' If "pairwise" regression is specified, \code{distribution} must be a list of
#' the form \code{list(name="pairwise",group=...,metric=...,max.rank=...)}
#' (\code{metric} and \code{max.rank} are optional, see below). \code{group} is
#' a character vector with the column names of \code{data} that jointly
#' indicate the group an instance belongs to (typically a query in Information
#' Retrieval applications). For training, only pairs of instances from the same
#' group and with different target labels can be considered. \code{metric} is
#' the IR measure to use, one of 
#' \describe{ 
#'   \item{list("conc")}{Fraction of concordant pairs; for binary labels, this 
#'         is equivalent to the Area under the ROC Curve}
#'   \item{:}{Fraction of concordant pairs; for binary labels, this is 
#'            equivalent to the Area under the ROC Curve} 
#'   \item{list("mrr")}{Mean reciprocal rank of the highest-ranked positive 
#'         instance}
#'   \item{:}{Mean reciprocal rank of the highest-ranked positive instance}
#'   \item{list("map")}{Mean average precision, a generalization of \code{mrr} 
#'         to multiple positive instances}\item{:}{Mean average precision, a
#'         generalization of \code{mrr} to multiple positive instances}
#'   \item{list("ndcg:")}{Normalized discounted cumulative gain. The score is 
#'         the weighted sum (DCG) of the user-supplied target values, weighted 
#'         by log(rank+1), and normalized to the maximum achievable value. This 
#'         is the default if the user did not specify a metric.} 
#' }
#' 
#' \code{ndcg} and \code{conc} allow arbitrary target values, while binary
#' targets {0,1} are expected for \code{map} and \code{mrr}. For \code{ndcg}
#' and \code{mrr}, a cut-off can be chosen using a positive integer parameter
#' \code{max.rank}. If left unspecified, all ranks are taken into account.
#' 
#' Note that splitting of instances into training and validation sets follows
#' group boundaries and therefore only approximates the specified
#' \code{train.fraction} ratio (the same applies to cross-validation folds).
#' Internally, queries are randomly shuffled before training, to avoid bias.
#' 
#' Weights can be used in conjunction with pairwise metrics, however it is
#' assumed that they are constant for instances from the same group.
#' 
#' For details and background on the algorithm, see e.g. Burges (2010).
#' 
#' @param w A vector of weights of the same length as the \code{y}.
#'
#' @param var.monotone an optional vector, the same length as the number of
#' predictors, indicating which variables have a monotone increasing (+1),
#' decreasing (-1), or arbitrary (0) relationship with the outcome.
#' 
#' @param n.trees the total number of trees to fit. This is equivalent to the
#' number of iterations and the number of basis functions in the additive
#' expansion.
#' 
#' @param interaction.depth The maximum depth of variable interactions. A value
#' of 1 implies an additive model, a value of 2 implies a model with up to 2-way 
#' interactions, etc. Default is \code{1}.
#' 
#' @param n.minobsinnode Integer specifying the minimum number of observations 
#' in the trees terminal nodes. Note that this is the actual number of 
#' observations not the total weight.
#' 
#' @param shrinkage The shrinkage parameter applied to each tree in the
#' expansion. Also known as the learning rate or step-size reduction; 0.001 to 
#' 0.1 usually work, but a smaller learning rate typically requires more trees.
#' Default is \code{0.1}.
#' 
#' @param bag.fraction The fraction of the training set observations randomly
#' selected to propose the next tree in the expansion. This introduces
#' randomnesses into the model fit. If \code{bag.fraction} < 1 then running the
#' same model twice will result in similar but different fits. \code{gbm} uses
#' the R random number generator so \code{set.seed} can ensure that the model
#' can be reconstructed. Preferably, the user can save the returned
#' \code{\link{gbm.object}} using \code{\link{save}}. Default is \code{0.5}.
#' 
#' @param nTrain An integer representing the number of cases on which to train.
#' This is the preferred way of specification for \code{gbm.fit}; The option
#' \code{train.fraction} in \code{gbm.fit} is deprecated and only maintained
#' for backward compatibility. These two parameters are mutually exclusive. If
#' both are unspecified, all data is used for training.
#' 
#' @param train.fraction The first \code{train.fraction * nrows(data)}
#' observations are used to fit the \code{gbm} and the remainder are used for
#' computing out-of-sample estimates of the loss function.
#' 
#' @param keep.data Logical indicating whether or not to keep the data and an 
#' index of the data stored with the object. Keeping the data and index makes 
#' subsequent calls to \code{\link{gbm.more}} faster at the cost of storing an 
#' extra copy of the dataset.
#' 
#' @param verbose Logical indicating whether or not to print out progress and 
#' performance indicators (\code{TRUE}). If this option is left unspecified for 
#' \code{gbm.more}, then it uses \code{verbose} from \code{object}. Default is
#' \code{FALSE}.
#' 
#' @param var.names Vector of strings of length equal to the number of columns 
#' of \code{x} containing the names of the predictor variables.
#' 
#' @param response.name Character string label for the response variable.
#' 
#' @param group The \code{group} to use when \code{distribution = "pairwise"}.
#' 
#' @return A \code{\link{gbm.object}} object.
#' 
#' @details 
#' This package implements the generalized boosted modeling framework. Boosting
#' is the process of iteratively adding basis functions in a greedy fashion so
#' that each additional basis function further reduces the selected loss
#' function. This implementation closely follows Friedman's Gradient Boosting
#' Machine (Friedman, 2001).
#' 
#' In addition to many of the features documented in the Gradient Boosting
#' Machine, \code{gbm} offers additional features including the out-of-bag
#' estimator for the optimal number of iterations, the ability to store and
#' manipulate the resulting \code{gbm} object, and a variety of other loss
#' functions that had not previously had associated boosting algorithms,
#' including the Cox partial likelihood for censored data, the poisson
#' likelihood for count outcomes, and a gradient boosting implementation to
#' minimize the AdaBoost exponential loss function.
#' 
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' 
#' Quantile regression code developed by Brian Kriegler
#' \email{bk@@stat.ucla.edu}
#' 
#' t-distribution, and multinomial code developed by Harry Southworth and
#' Daniel Edwards
#' 
#' Pairwise code developed by Stefan Schroedl \email{schroedl@@a9.com}
#' 
#' @seealso \code{\link{gbm.object}}, \code{\link{gbm.perf}}, 
#' \code{\link{plot.gbm}}, \code{\link{predict.gbm}}, \code{\link{summary.gbm}}, 
#' and \code{\link{pretty.gbm.tree}}.
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
#' B. Kriegler (2007). Cost-Sensitive Stochastic Gradient Boosting Within a 
#' Quantitative Regression Framework. Ph.D. Dissertation. University of 
#' California at Los Angeles, Los Angeles, CA, USA. Advisor(s) Richard A. Berk. 
#' url{https://dl.acm.org/citation.cfm?id=1354603}.
#' 
#' C. Burges (2010). \dQuote{From RankNet to LambdaRank to LambdaMART: An
#' Overview,} Microsoft Research Technical Report MSR-TR-2010-82.
#' 
#' @export
gbm.fit <- function(x, y, offset = NULL, misc = NULL, distribution = "bernoulli",
                    w = NULL, var.monotone = NULL, n.trees = 100, 
                    interaction.depth = 1, n.minobsinnode = 10, 
                    shrinkage = 0.001, bag.fraction = 0.5, nTrain = NULL,
                    train.fraction = NULL, keep.data = TRUE, verbose = TRUE,
                    var.names = NULL, response.name = "y", group = NULL) {
  
  # Reformat distribution into a named list
  if(is.character(distribution)) { 
    distribution <- list(name = distribution) 
  }
  
  # Dimensions of predictor data
  cRows <- nrow(x)
  cCols <- ncol(x)
  
  if(nrow(x) != ifelse(class(y) == "Surv", nrow(y), length(y))) {
    stop("The number of rows in x does not equal the length of y.")
  }
  
  # The preferred way to specify the number of training instances is via the 
  # parameter `nTrain`. The parameter `train.fraction` is only maintained for 
  # back compatibility.
  if(!is.null(nTrain) && !is.null(train.fraction)) {
    stop("Parameters `nTrain` and `train.fraction` cannot both be specified.")
  } else if(!is.null(train.fraction)) {
    warning("Parameter `train.fraction` is deprecated, please specify ",
            "`nTrain` instead.")
    nTrain <- floor(train.fraction*cRows)
  } else if(is.null(nTrain)) {
    nTrain <- cRows  # both undefined, use all training data
  }
  if (is.null(train.fraction)){
    train.fraction <- nTrain / cRows
  }
  
  # Extract var.names if NULL
  if(is.null(var.names)) {
    var.names <- getVarNames(x)
  }
  
  # Check size of data
  if(nTrain * bag.fraction <= 2 * n.minobsinnode + 1) {
    stop("The data set is too small or the subsampling rate is too large: ",
         "`nTrain * bag.fraction <= n.minobsinnode`")
  }
  
  if (distribution$name != "pairwise") {
    w <- w * length(w) / sum(w)  # normalize to N
  }
  
  # Sanity checks
  ch <- checkMissing(x, y)
  interaction.depth <- checkID(interaction.depth)
  w <- checkWeights(w, length(y))
  offset <- checkOffset(offset, y)
  
  Misc <- NA
  
  # setup variable types
  var.type <- rep(0,cCols)
  var.levels <- vector("list",cCols)
  for(i in 1:length(var.type))
  {
    if(all(is.na(x[,i])))
    {
      stop("variable ",i,": ",var.names[i]," has only missing values.")
    }
    if(is.ordered(x[,i]))
    {
      var.levels[[i]] <- levels(factor(x[,i]))
      x[,i] <- as.numeric(factor(x[,i]))-1
      var.type[i] <- 0
    }
    else if(is.factor(x[,i]))
    {
      if(length(levels(x[,i]))>1024)
        stop("gbm does not currently handle categorical variables with more than 1024 levels. Variable ",i,": ",var.names[i]," has ",length(levels(x[,i]))," levels.")
      var.levels[[i]] <- levels(factor(x[,i]))
      x[,i] <- as.numeric(factor(x[,i]))-1
      var.type[i] <- max(x[,i],na.rm=TRUE)+1
    }
    else if(is.numeric(x[,i]))
    {
      var.levels[[i]] <- quantile(x[,i],prob=(0:10)/10,na.rm=TRUE)
    }
    else
    {
      stop("variable ",i,": ",var.names[i]," is not of type numeric, ordered, or factor.")
    }
    
    # check for some variation in each variable
    if(length(unique(var.levels[[i]])) == 1)
    {
      warning("variable ",i,": ",var.names[i]," has no variation.")
    }
  }
  
  nClass <- 1
  
  if(!("name" %in% names(distribution))) {
    stop("The distribution is missing a `name` component; for example, ",
         "distribution = list(name = \"gaussian\").")
  }
  supported.distributions <- getAvailableDistributions()
  distribution.call.name <- distribution$name
  
  # Check for potential problems with the distribution
  if(!is.element(distribution$name, supported.distributions)) {
    stop("Distribution ",distribution$name," is not supported")
  }
  if((distribution$name == "bernoulli") && !all(is.element(y,0:1))) {
    stop("Bernoulli requires the response to be in {0,1}")
    if (is.factor(y)) {
      y <- as.integer(y) - 1
    }
  }
  if((distribution$name == "huberized") && !all(is.element(y,0:1))) {
    stop("Huberized square hinged loss requires the response to be in {0,1}")
    if (is.factor(y)) {
      y <- as.integer(y) - 1
    }
  }
  if((distribution$name == "poisson") && any(y<0)) {
    stop("Poisson requires the response to be positive")
  }
  if((distribution$name == "poisson") && any(y != trunc(y))) {
    stop("Poisson requires the response to be a positive integer")
  }
  if((distribution$name == "adaboost") && !all(is.element(y,0:1))) {
    stop("This version of AdaBoost requires the response to be in {0,1}")
    if (is.factor(y)) {
      y <- as.integer(y) - 1
    }
  }
  if(distribution$name == "quantile") {
    if(length(unique(w)) > 1) {
      stop("This version of gbm for the quantile regression lacks a weighted quantile. For now the weights must be constant.")
    }
    if(is.null(distribution$alpha)) {
      stop("For quantile regression, the distribution parameter must be a list with a parameter 'alpha' indicating the quantile, for example list(name=\"quantile\",alpha=0.95).")
    } else {
      if((distribution$alpha < 0) || (distribution$alpha > 1)) {
        stop("alpha must be between 0 and 1.")
      }
    }
    Misc <- c(alpha=distribution$alpha)
  }
  if(distribution$name == "coxph") {
    if(class(y)!="Surv") {
      stop("Outcome must be a survival object Surv(time,failure)")
    }
    if(attr(y,"type")!="right") {
      stop("gbm() currently only handles right censored observations")
    }
    Misc <- y[,2]
    y <- y[,1]
    
    # reverse sort the failure times to compute risk sets on the fly
    i.train <- order(-y[1:nTrain])
    n.test <- cRows - nTrain
    if(n.test > 0) {
      i.test <- order(-y[(nTrain+1):cRows]) + nTrain
    }
    else {
      i.test <- NULL
    }
    i.timeorder <- c(i.train,i.test)
    
    y <- y[i.timeorder]
    Misc <- Misc[i.timeorder]
    x <- x[i.timeorder,,drop=FALSE]
    w <- w[i.timeorder]
    if(!is.na(offset)) offset <- offset[i.timeorder]
  }
  if(distribution$name == "tdist") {
    if (is.null(distribution$df) || !is.numeric(distribution$df)){
      Misc <- 4
    }
    else {
      Misc <- distribution$df[1]
    }
  }
  if (distribution$name == "multinomial") {
    ## Ensure that the training set contains all classes
    classes <- attr(factor(y), "levels")
    nClass <- length(classes)
    
    if (nClass > nTrain) {
      stop(paste("Number of classes (", nClass, ") must be less than the",
                 " size of the training set (", nTrain, ").", sep = ""))
    }
    
    new.idx <- as.vector(sapply(classes, function(a,x){ min((1:length(x))[x==a]) }, y))
    
    all.idx <- 1:length(y)
    new.idx <- c(new.idx, all.idx[!(all.idx %in% new.idx)])
    
    y <- y[new.idx]
    x <- x[new.idx, ]
    w <- w[new.idx]
    if (!is.null(offset)) {
      offset <- offset[new.idx]
    }
    
    ## Get the factors
    y <- as.numeric(as.vector(outer(y, classes, "==")))
    
    ## Fill out the weight and offset
    w <- rep(w, nClass)
    if (!is.null(offset)) {
      offset <- rep(offset, nClass)
    }
  } # close if (dist... == "multinomial"
  
  if(distribution$name == "pairwise") {
    distribution.metric <- distribution[["metric"]]
    if (!is.null(distribution.metric)) {
      distribution.metric <- tolower(distribution.metric)
      supported.metrics <- c("conc", "ndcg", "map", "mrr")
      if (!is.element(distribution.metric, supported.metrics)) {
        stop("Metric '", distribution.metric, "' is not supported, use either 'conc', 'ndcg', 'map', or 'mrr'")
      }
      metric <- distribution.metric
    } else {
      warning("No metric specified, using 'ndcg'")
      metric <- "ndcg" # default
      distribution[["metric"]] <- metric
    }
    
    if (any(y<0)) {
      stop("targets for 'pairwise' should be non-negative")
    }
    
    if (is.element(metric, c("mrr", "map")) && (!all(is.element(y, 0:1)))) {
      stop("Metrics 'map' and 'mrr' require the response to be in {0,1}")
    }
    
    # Cut-off rank for metrics
    # Default of 0 means no cutoff
    
    max.rank <- 0
    if (!is.null(distribution[["max.rank"]]) && distribution[["max.rank"]] > 0) {
      if (is.element(metric, c("ndcg", "mrr"))) {
        max.rank <- distribution[["max.rank"]]
      }
      else {
        stop("Parameter 'max.rank' cannot be specified for metric '", distribution.metric, "', only supported for 'ndcg' and 'mrr'")
      }
    }
    
    # We pass the cut-off rank to the C function as the last element in the Misc vector
    Misc <- c(group, max.rank)
    
    distribution.call.name <- sprintf("pairwise_%s", metric)
  } # close if (dist... == "pairwise"
  
  # create index upfront... subtract one for 0 based order
  x.order <- apply(x[1:nTrain,,drop=FALSE],2,order,na.last=FALSE)-1
  
  x <- as.vector(data.matrix(x))
  predF <- rep(0,length(y))
  train.error <- rep(0,n.trees)
  valid.error <- rep(0,n.trees)
  oobag.improve <- rep(0,n.trees)
  
  if(is.null(var.monotone)) {
    var.monotone <- rep(0,cCols)
  } else if(length(var.monotone)!=cCols) {
    stop("Length of var.monotone != number of predictors")
  } else if(!all(is.element(var.monotone,-1:1))) {
    stop("var.monotone must be -1, 0, or 1")
  }
  fError <- FALSE
  
  gbm.obj <- .Call("gbm_fit",
                   Y=as.double(y),
                   Offset=as.double(offset),
                   X=as.double(x),
                   X.order=as.integer(x.order),
                   weights=as.double(w),
                   Misc=as.double(Misc),
                   cRows=as.integer(cRows),
                   cCols=as.integer(cCols),
                   var.type=as.integer(var.type),
                   var.monotone=as.integer(var.monotone),
                   distribution=as.character(distribution.call.name),
                   n.trees=as.integer(n.trees),
                   interaction.depth=as.integer(interaction.depth),
                   n.minobsinnode=as.integer(n.minobsinnode),
                   n.classes = as.integer(nClass),
                   shrinkage=as.double(shrinkage),
                   bag.fraction=as.double(bag.fraction),
                   nTrain=as.integer(nTrain),
                   fit.old=as.double(NA),
                   n.cat.splits.old=as.integer(0),
                   n.trees.old=as.integer(0),
                   verbose=as.integer(verbose),
                   PACKAGE = "gbm")
  
  names(gbm.obj) <- c("initF","fit","train.error","valid.error",
                      "oobag.improve","trees","c.splits")
  
  gbm.obj$bag.fraction <- bag.fraction
  gbm.obj$distribution <- distribution
  gbm.obj$interaction.depth <- interaction.depth
  gbm.obj$n.minobsinnode <- n.minobsinnode
  gbm.obj$num.classes <- nClass
  gbm.obj$n.trees <- length(gbm.obj$trees) / nClass
  gbm.obj$nTrain <- nTrain
  gbm.obj$train.fraction <- train.fraction
  gbm.obj$response.name <- response.name
  gbm.obj$shrinkage <- shrinkage
  gbm.obj$var.levels <- var.levels
  gbm.obj$var.monotone <- var.monotone
  gbm.obj$var.names <- var.names
  gbm.obj$var.type <- var.type
  gbm.obj$verbose <- verbose
  gbm.obj$Terms <- NULL
  
  if(distribution$name == "coxph") {
    gbm.obj$fit[i.timeorder] <- gbm.obj$fit
  }
  ## If K-Classification is used then split the fit and tree components
  if (distribution$name == "multinomial") {
    gbm.obj$fit <- matrix(gbm.obj$fit, ncol = nClass)
    dimnames(gbm.obj$fit)[[2]] <- classes
    gbm.obj$classes <- classes
    
    ## Also get the class estimators
    exp.f <- exp(gbm.obj$fit)
    denom <- matrix(rep(rowSums(exp.f), nClass), ncol = nClass)
    gbm.obj$estimator <- exp.f/denom
  }
  
  if(keep.data) {
    if(distribution$name == "coxph") {
      # Put the observations back in order
      gbm.obj$data <- list(
        y = y,
        x = x,
        x.order = x.order,
        offset = offset,
        Misc = Misc,
        w = w,
        i.timeorder = i.timeorder
      )
    }
    else if ( distribution$name == "multinomial" ) {
      # Restore original order of the data
      new.idx <- order(new.idx)
      gbm.obj$data <- list(
        y = as.vector(matrix(y, ncol = length(classes), byrow = FALSE)[new.idx, ]),
        x = as.vector(matrix(x, ncol = length(var.names), byrow = FALSE)[new.idx, ]),
        x.order = x.order,
        offset = offset[new.idx],
        Misc = Misc, 
        w = w[new.idx] 
      )
    } else {
      gbm.obj$data <- list(
        y = y,
        x = x,
        x.order = x.order,
        offset = offset,
        Misc = Misc,
        w = w
      )
    }
  }
  else {
    gbm.obj$data <- NULL
  }
  
  # Reuturn object of class "gbm"
  class(gbm.obj) <- "gbm"
  gbm.obj
  
}

