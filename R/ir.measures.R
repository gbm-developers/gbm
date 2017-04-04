# Functions to compute IR measures for pairwise loss for
# a single group
# Notes:
# * Inputs are passed as a 2-elemen (y,f) list, to
#   facilitate the 'by' iteration
# * Return the respective metric, or a negative value if
#   it is undefined for the given group
# * For simplicity, we have no special handling for ties;
#   instead, we break ties randomly. This is slightly
#   inaccurate for individual groups, but should have
#   a small effect on the overall measure.

#' Compute Information Retrieval measures.
#' 
#' Functions to compute Information Retrieval measures for pairwise loss for a
#' single group. The function returns the respective metric, or a negative
#' value if it is undefined for the given group.
#' 
#' @param obs Observed value.
#' @param pred Predicted value.
#' @param metric What type of performance measure to compute.
#' @param y,y.f,f,w,group,max.rank Used internally.
#' @param x ?.
#' @return The requested performance measure.
#'
#' @details
#' For simplicity, we have no special handling for ties; instead, we break ties
#' randomly. This is slightly inaccurate for individual groups, but should have
#' only a small effect on the overall measure.
#' 
#' \code{gbm.conc} computes the concordance index: Fraction of all pairs (i,j)
#' with i<j, x[i] != x[j], such that x[j] < x[i]
#' 
#' If \code{obs} is binary, then \code{gbm.roc.area(obs, pred) =
#' gbm.conc(obs[order(-pred)])}.
#' 
#' \code{gbm.conc} is more general as it allows non-binary targets, but is
#' significantly slower.
#' 
#' @aliases gbm.roc.area gbm.conc ir.measure.conc ir.measure.auc ir.measure.mrr
#' ir.measure.map ir.measure.ndcg perf.pairwise
#' 
#' @rdname gbm.roc.area
#' 
#' @author Stefan Schroedl
#' @seealso \code{\link{gbm}}
#' @references C. Burges (2010). "From RankNet to LambdaRank to LambdaMART: An
#' Overview", Microsoft Research Technical Report MSR-TR-2010-82.
#' @keywords models
#' 
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.


# Area under ROC curve = ratio of correctly ranking pairs
#' @rdname gbm.roc.area
#' @export
gbm.roc.area <- function(obs, pred) {
   n1 <- sum(obs)
   n <- length(obs)
   if (n==n1) { return(1) }
   # Fraction of concordant pairs
   # = sum_{pos}(rank-1) / #pairs with different labels
   # #pairs = n1 * (n-n1)
   return ((mean(rank(pred)[obs > 0]) - (n1 + 1)/2)/(n - n1))
}


# Concordance Index:
# Fraction of all pairs (i,j) with i<j, x[i] != x[j], such that x[j] < x[i]
# Invariant: if obs is binary, then
#      gbm.roc.area(obs, pred) = gbm.conc(obs[order(-pred)])
# gbm.conc is more general as it allows non-binary targets,
# but is significantly slower
#' @rdname gbm.roc.area
#' @export
gbm.conc <- function(x)
{
   lx <- length(x)
   return (sum(mapply(function(r) { sum(x[(r+1):lx]<x[r]) }, 1:(lx-1))))
}


#' @rdname gbm.roc.area
#' @export
ir.measure.conc <- function(y.f, max.rank=0)
{
   # Note: max.rank is meaningless for CONC

   y           <- y.f[[1]]
   f           <- y.f[[2]]

   tab         <- table(y)
   csum        <- cumsum(tab)
   total.pairs <- sum(tab * (csum - tab))

   if (total.pairs == 0)
   {
      return (-1.0)
   }
   else
   {
      return (gbm.conc(y[order(-f)]) / total.pairs)
   }
}


#' @rdname gbm.roc.area
#' @export
ir.measure.auc <- function(y.f, max.rank=0)
{
   # Note: max.rank is meaningless for AUC
   y       <- y.f[[1]]
   f       <- y.f[[2]]
   num.pos <- sum(y>0)

   if (length(f) <= 1 || num.pos == 0 || num.pos == length(f))
   {
      return (-1.0)
   }
   else
   {
      return (gbm.roc.area(obs=y, pred=f))
   }
}


#' @rdname gbm.roc.area
#' @export
ir.measure.mrr <- function(y.f, max.rank)
{
   y       <- y.f[[1]]
   f       <- y.f[[2]]
   num.pos <- sum(y>0)

   if (length(f) <= 1 || num.pos == 0 || num.pos == length(f))
   {
      return (-1.0)
   }

   ord         <- order(f, decreasing=TRUE)
   min.idx.pos <- min(which(y[ord]>0))

   if (min.idx.pos <= max.rank)
   {
      return (1.0 / min.idx.pos)
   }
   else
   {
      return (0.0)
   }
}


#' @rdname gbm.roc.area
#' @export
ir.measure.map <- function(y.f, max.rank=0)
{
   # Note: max.rank is meaningless for MAP

   y         <- y.f[[1]]
   f         <- y.f[[2]]
   ord       <- order(f, decreasing=TRUE)
   idx.pos   <- which(y[ord]>0)
   num.pos   <- length(idx.pos)

   if (length(f) <= 1 || num.pos == 0 || num.pos == length(f))
   {
      return (-1.0)
   }

   # Above and including the rank of the i-th positive result,
   # there are exactly i positives and rank(i) total results
   return (sum((1:length(idx.pos))/idx.pos) / num.pos)
}


#' @rdname gbm.roc.area
#' @export
ir.measure.ndcg <- function(y.f, max.rank)
{
   y         <- y.f[[1]]
   f         <- y.f[[2]]

   if (length(f) <= 1 || all(diff(y)==0))
   {
      return (-1.0)
   }

   num.items <- min(length(f), max.rank)
   ord       <- order(f, decreasing=TRUE)

   dcg       <- sum(y[ord][1:num.items] / log2(2:(num.items+1)))

   # The best possible DCG: order by target
   ord.max   <- order(y, decreasing=TRUE)
   dcg.max   <- sum(y[ord.max][1:num.items] / log2(2:(num.items+1)))

   # Normalize
   return (dcg / dcg.max)
}


#' @rdname gbm.roc.area
#' @export
perf.pairwise <- function(y, f, group, metric="ndcg", w=NULL, max.rank=0)
{
  func.name <- switch(metric,
                      conc = "ir.measure.conc",
                      mrr  = "ir.measure.mrr",
                      map  = "ir.measure.map",
                      ndcg = "ir.measure.ndcg",
                      stop(paste("Metric",metric,"is not supported"))
  )
  
  # Optimization: for binary targets,
  # AUC is equivalent but faster than CONC
  if (metric == "conc" && all(is.element(y, 0:1)))
  {
    func.name <- "ir.measure.auc"
  }
  
  # Max rank = 0 means no cut off
  if (max.rank <= 0)
  {
    max.rank <- length(y)+1
  }
  
  # Random tie breaking in case of duplicate scores.
  # (Without tie breaking, we would overestimate if instances are
  # sorted descending on target)
  f <- f + 1E-10 * runif(length(f), min=-0.5, max=0.5)
  
  measure.by.group <- as.matrix(by(list(y, f), INDICES=group, FUN=get(func.name), max.rank=max.rank))
  
  # Exclude groups with single result or only negative or positive instances
  idx <- which((!is.null(measure.by.group)) & measure.by.group >= 0)
  
  if (is.null(w))
  {
    return (mean(measure.by.group[idx]))
  }
  else
  {
    # Assumption: weights are constant per group
    w.by.group <- tapply(w, group, mean)
    return (weighted.mean(measure.by.group[idx], w=w.by.group[idx]))
  }
}
