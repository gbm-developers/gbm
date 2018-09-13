#' Generalized Boosted Regression Modeling (GBM)
#' 
#' Adds additional trees to a \code{\link{gbm.object}} object.
#' 
#' @param object A \code{\link{gbm.object}} object created from an initial call 
#' to \code{\link{gbm}}.
#' 
#' @param n.new.trees Integer specifying the number of additional trees to add 
#' to \code{object}. Default is 100.
#' 
#' @param data An optional data frame containing the variables in the model. By
#' default the variables are taken from \code{environment(formula)}, typically
#' the environment from which \code{gbm} is called. If \code{keep.data=TRUE} in
#' the initial call to \code{gbm} then \code{gbm} stores a copy with the
#' object. If \code{keep.data=FALSE} then subsequent calls to
#' \code{\link{gbm.more}} must resupply the same dataset. It becomes the user's
#' responsibility to resupply the same data at this point.
#' 
#' @param weights An optional vector of weights to be used in the fitting
#' process. Must be positive but do not need to be normalized. If
#' \code{keep.data=FALSE} in the initial call to \code{gbm} then it is the
#' user's responsibility to resupply the weights to \code{\link{gbm.more}}.
#' 
#' @param offset A vector of offset values.
#' 
#' @param verbose Logical indicating whether or not to print out progress and 
#' performance indicators (\code{TRUE}). If this option is left unspecified for 
#' \code{gbm.more}, then it uses \code{verbose} from \code{object}. Default is
#' \code{FALSE}.
#' 
#' @return A \code{\link{gbm.object}} object.
#'
#' @export
#' 
#' @examples
#' #
#' # A least squares regression example 
#' #
#' 
#' # Simulate data
#' set.seed(101)  # for reproducibility
#' N <- 1000
#' X1 <- runif(N)
#' X2 <- 2 * runif(N)
#' X3 <- ordered(sample(letters[1:4], N, replace = TRUE), levels = letters[4:1])
#' X4 <- factor(sample(letters[1:6], N, replace = TRUE))
#' X5 <- factor(sample(letters[1:3], N, replace = TRUE))
#' X6 <- 3 * runif(N) 
#' mu <- c(-1, 0, 1, 2)[as.numeric(X3)]
#' SNR <- 10  # signal-to-noise ratio
#' Y <- X1 ^ 1.5 + 2 * (X2 ^ 0.5) + mu
#' sigma <- sqrt(var(Y) / SNR)
#' Y <- Y + rnorm(N, 0, sigma)
#' X1[sample(1:N,size=500)] <- NA  # introduce some missing values
#' X4[sample(1:N,size=300)] <- NA  # introduce some missing values
#' data <- data.frame(Y, X1, X2, X3, X4, X5, X6)
#' 
#' # Fit a GBM
#' set.seed(102)  # for reproducibility
#' gbm1 <- gbm(Y ~ ., data = data, var.monotone = c(0, 0, 0, 0, 0, 0),
#'             distribution = "gaussian", n.trees = 100, shrinkage = 0.1,             
#'             interaction.depth = 3, bag.fraction = 0.5, train.fraction = 0.5,  
#'             n.minobsinnode = 10, cv.folds = 5, keep.data = TRUE, 
#'             verbose = FALSE, n.cores = 1)  
#' 
#' # Check performance using the out-of-bag (OOB) error; the OOB error typically
#' # underestimates the optimal number of iterations
#' best.iter <- gbm.perf(gbm1, method = "OOB")
#' print(best.iter)
#' 
#' # Check performance using the 50% heldout test set
#' best.iter <- gbm.perf(gbm1, method = "test")
#' print(best.iter)
#' 
#' # Check performance using 5-fold cross-validation
#' best.iter <- gbm.perf(gbm1, method = "cv")
#' print(best.iter)
#' 
#' # Plot relative influence of each variable
#' par(mfrow = c(1, 2))
#' summary(gbm1, n.trees = 1)          # using first tree
#' summary(gbm1, n.trees = best.iter)  # using estimated best number of trees
#' 
#' # Compactly print the first and last trees for curiosity
#' print(pretty.gbm.tree(gbm1, i.tree = 1))
#' print(pretty.gbm.tree(gbm1, i.tree = gbm1$n.trees))
#' 
#' # Simulate new data
#' set.seed(103)  # for reproducibility
#' N <- 1000
#' X1 <- runif(N)
#' X2 <- 2 * runif(N)
#' X3 <- ordered(sample(letters[1:4], N, replace = TRUE))
#' X4 <- factor(sample(letters[1:6], N, replace = TRUE))
#' X5 <- factor(sample(letters[1:3], N, replace = TRUE))
#' X6 <- 3 * runif(N) 
#' mu <- c(-1, 0, 1, 2)[as.numeric(X3)]
#' Y <- X1 ^ 1.5 + 2 * (X2 ^ 0.5) + mu + rnorm(N, 0, sigma)
#' data2 <- data.frame(Y, X1, X2, X3, X4, X5, X6)
#' 
#' # Predict on the new data using the "best" number of trees; by default,
#' # predictions will be on the link scale
#' Yhat <- predict(gbm1, newdata = data2, n.trees = best.iter, type = "link")
#' 
#' # least squares error
#' print(sum((data2$Y - Yhat)^2))
#' 
#' # Construct univariate partial dependence plots
#' p1 <- plot(gbm1, i.var = 1, n.trees = best.iter)
#' p2 <- plot(gbm1, i.var = 2, n.trees = best.iter)
#' p3 <- plot(gbm1, i.var = "X3", n.trees = best.iter)  # can use index or name
#' grid.arrange(p1, p2, p3, ncol = 3)
#' 
#' # Construct bivariate partial dependence plots
#' plot(gbm1, i.var = 1:2, n.trees = best.iter)
#' plot(gbm1, i.var = c("X2", "X3"), n.trees = best.iter)
#' plot(gbm1, i.var = 3:4, n.trees = best.iter)
#' 
#' # Construct trivariate partial dependence plots
#' plot(gbm1, i.var = c(1, 2, 6), n.trees = best.iter, 
#'      continuous.resolution = 20)
#' plot(gbm1, i.var = 1:3, n.trees = best.iter)
#' plot(gbm1, i.var = 2:4, n.trees = best.iter)
#' plot(gbm1, i.var = 3:5, n.trees = best.iter)
#' 
#' # Add more (i.e., 100) boosting iterations to the ensemble
#' gbm2 <- gbm.more(gbm1, n.new.trees = 100, verbose = FALSE)
gbm.more <- function(object,
                     n.new.trees = 100,
                     data = NULL,
                     weights = NULL,
                     offset = NULL,
                     verbose = NULL) {
   theCall <- match.call()
   nTrain  <- object$nTrain

   if (object$distribution$name != "pairwise")
   {
      distribution.call.name <- object$distribution$name
   }
   else
   {
      distribution.call.name <- sprintf("pairwise_%s", object$distribution$metric)
   }

   if(is.null(object$Terms) && is.null(object$data))
   {
      stop("The gbm model was fit using gbm.fit (rather than gbm) and keep.data was set to FALSE. gbm.more cannot locate the dataset.")
   }
   else if(is.null(object$data) && is.null(data))
   {
      stop("keep.data was set to FALSE on original gbm call and argument 'data' is NULL")
   }
   else if(is.null(object$data))
   {
      m <- eval(object$m, parent.frame())

      Terms <- attr(m, "terms")
      a <- attributes(Terms)

      y <- as.vector(model.extract(m, "response"))
      offset <- model.extract(m,offset)
      x <- model.frame(delete.response(Terms),
                       data,
                       na.action=na.pass)

      w <- weights
      if(length(w)==0) w <- rep(1, nrow(x))

      if (object$distribution$name != "pairwise")
      {
         w <- w*length(w)/sum(w) # normalize to N
      }

      if(is.null(offset) || (offset==0))
      {
         offset <- NA
      }
      Misc <- NA

      if(object$distribution$name == "coxph")
      {
         Misc <- as.numeric(y)[-(1:cRows)]
         y <- as.numeric(y)[1:cRows]

         # reverse sort the failure times to compute risk sets on the fly
         i.train <- order(-y[1:nTrain])
         i.test <- order(-y[(nTrain+1):cRows]) + nTrain
         i.timeorder <- c(i.train,i.test)

         y <- y[i.timeorder]
         Misc <- Misc[i.timeorder]
         x <- x[i.timeorder,,drop=FALSE]
         w <- w[i.timeorder]
         if(!is.na(offset)) offset <- offset[i.timeorder]
         object$fit <- object$fit[i.timeorder]
      }
      else if(object$distribution$name == "tdist" ){
         Misc <- object$distribution$df
      }
      else if (object$distribution$name == "pairwise"){

         # Check if group names are valid
         distribution.group <- object$distribution$group
         i <- match(distribution.group, colnames(data))
         if (any(is.na(i)))
         {
            stop("Group column does not occur in data: ", distribution.group[is.na(i)])
         }

         # construct group index
         group <- factor(do.call(paste, c(data[,distribution.group, drop=FALSE], sep=":")))

         # Check that weights are constant across groups
         if ((!missing(weights)) && (!is.null(weights)))
         {
            w.min <- tapply(w, INDEX=group, FUN=min)
            w.max <- tapply(w, INDEX=group, FUN=max)

            if (any(w.min != w.max))
            {
               stop("For distribution 'pairwise', all instances for the same group must have the same weight")
            }

            # Normalize across groups
            w <- w * length(w.min) / sum(w.min)
         }

         # Shuffle groups, to remove bias when splitting into train/test set and/or CV folds
         perm.levels  <- levels(group)[sample(1:nlevels(group))]
         group        <- factor(group, levels=perm.levels)

         # The C function expects instances to be sorted by group and descending by target

         ord.group    <- object$ord.group
         group        <- group[ord.group]
         y            <- y[ord.group]
         x            <- x[ord.group,,drop=FALSE]
         w            <- x[ord.group]
         object$fit   <- object$fit[ord.group] # object$fit is stored in the original order

         # Split into train and validation set, at group boundary
         num.groups.train <- max(1, round(object$train.fraction * nlevels(group)))

         # include all groups up to the num.groups.train
         nTrain           <- max(which(group==levels(group)[num.groups.train]))

         metric <- object$distribution[["metric"]]

         if (is.element(metric, c("mrr", "map")) && (!all(is.element(y, 0:1))))
         {
            stop("Metrics 'map' and 'mrr' require the response to be in {0,1}")
         }

         # Cut-off rank for metrics
         # We pass this argument as the last element in the Misc vector
         # Default of 0 means no cutoff

         max.rank <- 0
         if (!is.null(object$distribution[["max.rank"]]) && object$distribution[["max.rank"]] > 0)
         {
            if (is.element(metric, c("ndcg", "mrr")))
            {
               max.rank <- object$distribution[["max.rank"]]
            }
            else
            {
               stop("Parameter 'max.rank' cannot be specified for metric '", metric, "', only supported for 'ndcg' and 'mrr'")
            }
         }

         Misc <- c(group, max.rank)

      }

      # create index upfront... subtract one for 0 based order
      x.order <- apply(x[1:nTrain,,drop=FALSE],2,order,na.last=FALSE)-1
      x <- data.matrix(x)
      cRows <- nrow(x)
      cCols <- ncol(x)
   }
   else
   {
      y       <- object$data$y
      x       <- object$data$x
      x.order <- object$data$x.order
      offset  <- object$data$offset
      Misc    <- object$data$Misc
      w       <- object$data$w
      nTrain  <- object$nTrain
      cRows   <- length(y)
      cCols   <- length(x)/cRows
      if(object$distribution$name == "coxph")
      {
         i.timeorder <- object$data$i.timeorder
         object$fit  <- object$fit[i.timeorder]
      }
      if (object$distribution$name == "pairwise") 
      {
         object$fit   <- object$fit[object$ord.group] # object$fit is stored in the original order
      }
   }

   if(is.null(verbose))
   {
      verbose <- object$verbose
   }
   x <- as.vector(x)

   gbm.obj <- .Call("gbm_fit",
                    Y = as.double(y),
                    Offset = as.double(offset),
                    X = as.double(x),
                    X.order = as.integer(x.order),
                    weights = as.double(w),
                    Misc = as.double(Misc),
                    cRows = as.integer(cRows),
                    cCols = as.integer(cCols),
                    var.type = as.integer(object$var.type),
                    var.monotone = as.integer(object$var.monotone),
                    distribution = as.character(distribution.call.name),
                    n.trees = as.integer(n.new.trees),
                    interaction.depth = as.integer(object$interaction.depth),
                    n.minobsinnode = as.integer(object$n.minobsinnode),
                    n.classes = as.integer(object$num.classes),
                    shrinkage = as.double(object$shrinkage),
                    bag.fraction = as.double(object$bag.fraction),
                    train.fraction = as.integer(nTrain),
                    fit.old = as.double(object$fit),
                    n.cat.splits.old = as.integer(length(object$c.splits)),
                    n.trees.old = as.integer(object$n.trees),
                    verbose = as.integer(verbose),
                    PACKAGE = "gbm")
   names(gbm.obj) <- c("initF","fit","train.error","valid.error",
                       "oobag.improve","trees","c.splits")

   gbm.obj$initF         <- object$initF
   gbm.obj$train.error   <- c(object$train.error, gbm.obj$train.error)
   gbm.obj$valid.error   <- c(object$valid.error, gbm.obj$valid.error)
   gbm.obj$oobag.improve <- c(object$oobag.improve, gbm.obj$oobag.improve)
   gbm.obj$trees         <- c(object$trees, gbm.obj$trees)
   gbm.obj$c.splits      <- c(object$c.splits, gbm.obj$c.splits)

   # cv.error not updated when using gbm.more
   gbm.obj$cv.error      <- object$cv.error
   gbm.obj$cv.folds      <- object$cv.folds

   gbm.obj$n.trees        <- length(gbm.obj$trees)
   gbm.obj$distribution   <- object$distribution
   gbm.obj$train.fraction <- object$train.fraction
   gbm.obj$shrinkage      <- object$shrinkage
   gbm.obj$bag.fraction   <- object$bag.fraction
   gbm.obj$var.type       <- object$var.type
   gbm.obj$var.monotone   <- object$var.monotone
   gbm.obj$var.names      <- object$var.names
   gbm.obj$interaction.depth <- object$interaction.depth
   gbm.obj$n.minobsinnode    <- object$n.minobsinnode
   gbm.obj$num.classes       <- object$num.classes
   gbm.obj$nTrain            <- object$nTrain
   gbm.obj$response.name     <- object$response.name
   gbm.obj$Terms             <- object$Terms
   gbm.obj$var.levels        <- object$var.levels
   gbm.obj$verbose           <- verbose

   if(object$distribution$name == "coxph")
   {
      gbm.obj$fit[i.timeorder] <- gbm.obj$fit
   }

   if (object$distribution$name == "pairwise")
   {
      # Data has been reordered according to queries.
      # We need to permute the fitted values to correspond
      # to the original order.
      gbm.obj$fit <- gbm.obj$fit[order(object$ord.group)]
      object$fit  <- object$fit[order(object$ord.group)]
      gbm.obj$ord.group <- object$ord.group
   }  

   if(!is.null(object$data))
   {
      gbm.obj$data <- object$data
   }
   else
   {
      gbm.obj$data <- NULL
   }
   gbm.obj$m <- object$m
   gbm.obj$call <- theCall

   class(gbm.obj) <- "gbm"
   return(gbm.obj)
}
