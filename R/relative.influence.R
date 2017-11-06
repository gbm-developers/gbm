#' Methods for estimating relative influence
#' 
#' Helper functions for computing the relative influence of each variable in
#' the gbm object.
#' 
#' @details
#' This is not intended for end-user use. These functions offer the different
#' methods for computing the relative influence in \code{\link{summary.gbm}}.
#' \code{gbm.loss} is a helper function for \code{permutation.test.gbm}.
#' 
#' @aliases relative.influence permutation.test.gbm gbm.loss
#' 
#' @param object a \code{gbm} object created from an initial call to
#' \code{\link{gbm}}.
#' 
#' @param n.trees the number of trees to use for computations. If not provided,
#' the the function will guess: if a test set was used in fitting, the number
#' of trees resulting in lowest test set error will be used; otherwise, if
#' cross-validation was performed, the number of trees resulting in lowest
#' cross-validation error will be used; otherwise, all trees will be used.
#' 
#' @param scale.  whether or not the result should be scaled. Defaults to
#' \code{FALSE}.
#' 
#' @param sort.  whether or not the results should be (reverse) sorted.
#' Defaults to \code{FALSE}.
#' 
#' @param y,f,w,offset,dist,baseline For \code{gbm.loss}: These components are
#' the outcome, predicted value, observation weight, offset, distribution, and
#' comparison loss function, respectively.
#' 
#' @param group,max.rank Used internally when \code{distribution =
#' \'pairwise\'}.
#' 
#' @return By default, returns an unprocessed vector of estimated relative
#' influences. If the \code{scale.} and \code{sort.} arguments are used,
#' returns a processed version of the same.
#' 
#' @author Greg Ridgeway \email{gregridgeway@@gmail.com}
#' 
#' @seealso \code{\link{summary.gbm}}
#' 
#' @references J.H. Friedman (2001). "Greedy Function Approximation: A Gradient
#' Boosting Machine," Annals of Statistics 29(5):1189-1232.
#' 
#' L. Breiman (2001).
#' \url{https://www.stat.berkeley.edu/users/breiman/randomforest2001.pdf}.
#' 
#' @keywords hplot
#' 
#' @rdname relative.influence
#' 
#' @export
relative.influence <- function(object,
                               n.trees,
                               scale. = FALSE,
                               sort. = FALSE )
{

   if( missing( n.trees ) ){
      if ( object$train.fraction < 1 ){
         n.trees <- gbm.perf( object, method="test", plot.it=FALSE )
      }
      else if ( !is.null( object$cv.error ) ){
         n.trees <- gbm.perf( object, method="cv", plot.it = FALSE )
      }
      else{
         # If dist=multinomial, object$n.trees = n.trees * num.classes
         # so use the following instead.
         n.trees <- length( object$train.error )
      }
      cat( paste( "n.trees not given. Using", n.trees, "trees.\n" ) )
      if (object$distribution == "multinomial"){
          n.trees <- n.trees * object$num.classes
      }
   }
   get.rel.inf <- function(obj)
   {
      lapply(split(obj[[6]],obj[[1]]),sum) # 6 - Improvement, 1 - var name
   }

   temp <- unlist(lapply(object$trees[1:n.trees],get.rel.inf))
   rel.inf.compact <- unlist(lapply(split(temp,names(temp)),sum))
   rel.inf.compact <- rel.inf.compact[names(rel.inf.compact)!="-1"]

   # rel.inf.compact excludes those variable that never entered the model
   # insert 0's for the excluded variables
   rel.inf <- rep(0,length(object$var.names))
   i <- as.numeric(names(rel.inf.compact))+1
   rel.inf[i] <- rel.inf.compact

   names(rel.inf) <- object$var.names

   if (scale.){
      rel.inf <- rel.inf / max(rel.inf)
   }
   if (sort.){
      rel.inf <- rev(sort(rel.inf))
   }

   return(rel.inf=rel.inf)
}


#' @rdname relative.influence
#' @export
permutation.test.gbm <- function(object,
                                 n.trees)
{
  # get variables used in the model
  i.vars <- sort(unique(unlist(lapply(object$trees[1:n.trees],
                                      function(x){unique(x[[1]])}))))
  i.vars <- i.vars[i.vars!=-1] + 1
  rel.inf <- rep(0,length(object$var.names))
  
  if(!is.null(object$data))
  {
    y            <- object$data$y
    os           <- object$data$offset
    Misc         <- object$data$Misc
    w            <- object$data$w
    x            <- matrix(object$data$x, ncol=length(object$var.names))
    object$Terms <- NULL # this makes predict.gbm take x as it is
    
    if (object$distribution$name == "pairwise")
    {
      # group and cutoff are only relevant for distribution "pairwise"
      # in this case, the last element specifies the max rank
      # max rank = 0 means no cut off
      group     <- Misc[1:length(y)]
      max.rank  <- Misc[length(y)+1]
    }
  }
  else
  {
    stop("Model was fit with keep.data=FALSE. permutation.test.gbm has not been implemented for that case.")
  }
  
  # the index shuffler
  j <- sample(1:nrow(x))
  for(i in 1:length(i.vars))
  {
    x[ ,i.vars[i]]  <- x[j,i.vars[i]]
    
    new.pred <- predict.gbm(object,newdata=x,n.trees=n.trees)
    rel.inf[i.vars[i]] <- gbm.loss(y,new.pred,w,os,
                                   object$distribution,
                                   object$train.error[n.trees],
                                   group,
                                   max.rank)
    
    x[j,i.vars[i]] <- x[ ,i.vars[i]]
  }
  
  return(rel.inf=rel.inf)
}


#' @rdname relative.influence
#' @export
gbm.loss <- function(y, f, w, offset, dist, baseline, group=NULL, max.rank=NULL)
{
  if (!is.na(offset))
  {
    f <- offset+f
  }
  
  if (dist$name != "pairwise")
  {
    switch(dist$name,
           gaussian = weighted.mean((y - f)^2,w) - baseline,
           bernoulli = -2*weighted.mean(y*f - log(1+exp(f)),w) - baseline,
           laplace = weighted.mean(abs(y-f),w) - baseline,
           adaboost = weighted.mean(exp(-(2*y-1)*f),w) - baseline,
           poisson = -2*weighted.mean(y*f-exp(f),w) - baseline,
           stop(paste("Distribution",dist$name,"is not yet supported for method=permutation.test.gbm")))
  }
  else # dist$name == "pairwise"
  {
    if (is.null(dist$metric))
    {
      stop("No metric specified for distribution 'pairwise'")
    }
    if (!is.element(dist$metric, c("conc", "ndcg", "map", "mrr")))
    {
      stop("Invalid metric '", dist$metric, "' specified for distribution 'pairwise'")
    }
    if (is.null(group))
    {
      stop("For distribution 'pairwise', parameter 'group' has to be supplied")
    }
    # Loss = 1 - utility
    (1 - perf.pairwise(y, f, group, dist$metric, w, max.rank)) - baseline
  }
}
