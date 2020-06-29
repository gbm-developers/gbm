pkgname <- "inTrees"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('inTrees')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("GBM2List")
### * GBM2List

flush(stderr()); flush(stdout())

### Name: GBM2List
### Title: Transform gbm object to a list of trees
### Aliases: GBM2List
### Keywords: gbm

### ** Examples

    library(gbm)
    data(iris)
    X <- iris[,1:(ncol(iris)-1)]
    target <- iris[,"Species"] 
    gbmFit <- gbm(Species~ ., data=iris, n.tree = 400,
                    interaction.depth = 10,distribution="multinomial")
    treeList <- GBM2List(gbmFit,X)
    ruleExec = extractRules(treeList,X)
    ruleExec <- unique(ruleExec)
    #ruleExec <- ruleExec[1:min(2000,length(ruleExec)),,drop=FALSE]
    ruleMetric <- getRuleMetric(ruleExec,X,target)
    ruleMetric <- pruneRule(ruleMetric,X,target)
    ruleMetric <- unique(ruleMetric)
    learner <- buildLearner(ruleMetric,X,target)
    pred <- applyLearner(learner,X)
    readableLearner <- presentRules(learner,colnames(X)) # more readable format
    err <- 1-sum(pred==target)/length(pred);



cleanEx()
nameEx("Num2Level")
### * Num2Level

flush(stderr()); flush(stdout())

### Name: Num2Level
### Title: internal function
### Aliases: Num2Level
### Keywords: internal

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rfList, splitV) 
{
    for (i in 1:rfList$ntree) {
        rfList$list[[i]] <- data.frame(rfList$list[[i]])
        rfList$list[[i]][, "prediction"] <- data.frame(dicretizeVector(rfList$list[[i]][, 
            "prediction"], splitV))
        colnames(rfList$list[[i]]) <- c("left daughter", "right daughter", 
            "split var", "split point", "status", "prediction")
    }
    return(rfList)
  }



cleanEx()
nameEx("RF2List")
### * RF2List

flush(stderr()); flush(stdout())

### Name: RF2List
### Title: Transform a random forest object to a list of trees
### Aliases: RF2List
### Keywords: randomforest

### ** Examples

library(RRF)
data(iris)
X <- iris[,1:(ncol(iris)-1)]
target <- iris[,"Species"] 
rf <- RRF(X,as.factor(target),ntree=100) # build an ordinary RF 
treeList <- RF2List(rf)
ruleExec <- extractRules(treeList,X) # transform to R-executable rules



cleanEx()
nameEx("XGB2List")
### * XGB2List

flush(stderr()); flush(stdout())

### Name: XGB2List
### Title: Transform an xgboost object to a list of trees
### Aliases: XGB2List
### Keywords: xgboost

### ** Examples

	library(data.table)
	library(xgboost)
	# test data set 1: iris
	X <- within(iris,rm("Species")); Y <- iris[,"Species"]
	X <- within(iris,rm("Species")); Y <- iris[,"Species"]
	model_mat <- model.matrix(~. -1, data=X)
	xgb <- xgboost(model_mat, label = as.numeric(Y) - 1, nrounds = 20, 
		objective = "multi:softprob", num_class = 3 )
	tree_list <- XGB2List(xgb,model_mat)



cleanEx()
nameEx("applyLearner")
### * applyLearner

flush(stderr()); flush(stdout())

### Name: applyLearner
### Title: apply a simplified tree ensemble learner (STEL) to data
### Aliases: applyLearner
### Keywords: apply predict

### ** Examples

# see function "buildLearner" for examples
# pred <- applyLearner(learner,X)



cleanEx()
nameEx("buildLearner")
### * buildLearner

flush(stderr()); flush(stdout())

### Name: buildLearner
### Title: build a simplified tree ensemble learner (STEL)
### Aliases: buildLearner
### Keywords: STEL learner

### ** Examples

data(iris)
library(RRF)
X <- iris[,1:(ncol(iris)-1)]
target <- iris[,"Species"] 
rf <- RRF(X,as.factor(target),ntree=100) # build an ordinary RF 
treeList <- RF2List(rf)
ruleExec <- extractRules(treeList,X)
ruleExec <- unique(ruleExec)
ruleMetric <- getRuleMetric(ruleExec,X,target) # measure rules
ruleMetric <- pruneRule(ruleMetric,X,target) # prune each rule
#ruleMetric <- selectRuleRRF(ruleMetric,X,target) # rule selection
learner <- buildLearner(ruleMetric,X,target)
pred <- applyLearner(learner,X)
read <- presentRules(learner,colnames(X)) # more readable format

# format the rule and metrics as a table in latex code
library(xtable)
print(xtable(read), include.rownames=FALSE)
print(xtable(ruleMetric[1:2,]), include.rownames=FALSE)




cleanEx()
nameEx("computeRuleInfor")
### * computeRuleInfor

flush(stderr()); flush(stdout())

### Name: computeRuleInfor
### Title: compute rule information
### Aliases: computeRuleInfor
### Keywords: internal

### ** Examples

	# this is an internal function.



cleanEx()
nameEx("dataSimulate")
### * dataSimulate

flush(stderr()); flush(stdout())

### Name: dataSimulate
### Title: Simulate data
### Aliases: dataSimulate
### Keywords: simulate

### ** Examples

res <- dataSimulate(flag=1)
X <- res$X; 
target <- res$target



cleanEx()
nameEx("dicretizeVector")
### * dicretizeVector

flush(stderr()); flush(stdout())

### Name: dicretizeVector
### Title: discretize a variable
### Aliases: dicretizeVector
### Keywords: discretize

### ** Examples

 data(iris)
 dicretizeVector(iris[,1],3)



cleanEx()
nameEx("extractRules")
### * extractRules

flush(stderr()); flush(stdout())

### Name: extractRules
### Title: Extract rules from a list of trees
### Aliases: extractRules
### Keywords: extract

### ** Examples

    library(RRF)
    data(iris)
    X <- iris[,1:(ncol(iris)-1)]
    target <- iris[,"Species"] 
    rf <- RRF(X,as.factor(target),ntree=100) # build an ordinary RF 
    treeList <- RF2List(rf)
    ruleExec <- extractRules(treeList,X,digits=4) # transform to R-executable rules
    ruleExec <- unique(ruleExec)



cleanEx()
nameEx("formatGBM")
### * formatGBM

flush(stderr()); flush(stdout())

### Name: formatGBM
### Title: internal
### Aliases: formatGBM
### Keywords: internal

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (gbmList, splitBin,X) 
{
    for (j in 1:length(gbmList$list)) {
        a <- gbmList$list[[j]]
        rownames(a) <- 1:nrow(a)
        a$status <- a$SplitVar
        a <- a[, c("LeftNode", "RightNode", "MissingNode", "SplitVar", 
            "SplitCodePred", "status")]
        a[which(a[, "SplitVar"] >= 0), c("SplitVar", "LeftNode", 
            "RightNode", "MissingNode")] <- a[which(a[, "SplitVar"] >= 
            0), c("SplitVar", "LeftNode", "RightNode", "MissingNode")] + 
            1
        ix <- a$MissingNode[which(a$MissingNode > 0)]
        if (length(ix) > 0) 
            a$status[ix] <- 10
        a <- a[, c("LeftNode", "RightNode", "SplitVar", "SplitCodePred", 
            "status")]
        cat <- which(sapply(X, is.factor) & !sapply(X, is.ordered))
        ix <- which(a[, "SplitVar"] %in% cat)
        for (i in ix) a[i, "SplitCodePred"] <- splitBin[a[i, 
            "SplitCodePred"] + 1]
        colnames(a) <- c("left daughter", "right daughter", "split var", 
            "split point", "status")
        gbmList$list[[j]] <- a
    }
    return(gbmList)
  }



cleanEx()
nameEx("getFreqPattern")
### * getFreqPattern

flush(stderr()); flush(stdout())

### Name: getFreqPattern
### Title: calculate frequent variable interactions
### Aliases: getFreqPattern
### Keywords: variable interaction

### ** Examples

library(RRF)
library(arules)
data(iris)
X <- iris[,1:(ncol(iris)-1)]
target <- iris[,"Species"] 
rf <- RRF(X,as.factor(target),ntree=100) # build an ordinary RF 
treeList <- RF2List(rf)
ruleExec <- extractRules(treeList,X) # transform to R-executable rules
ruleMetric <- getRuleMetric(ruleExec,X,target) 
freqPattern <- getFreqPattern(ruleMetric)
freqPatternMetric <- getRuleMetric(freqPattern,X,target)



cleanEx()
nameEx("getRuleMetric")
### * getRuleMetric

flush(stderr()); flush(stdout())

### Name: getRuleMetric
### Title: Assign outcomes to a conditions, and measure the rules
### Aliases: getRuleMetric
### Keywords: measure rank

### ** Examples

library(RRF)
data(iris)
X <- iris[,1:(ncol(iris)-1)]
target <- iris[,"Species"] 
rf <- RRF(X,as.factor(target),ntree=100) # build an ordinary RF 
treeList <- RF2List(rf)
ruleExec <- extractRules(treeList,X) # transform to R-executable rules
ruleExec <- unique(ruleExec)
ruleMetric <- getRuleMetric(ruleExec,X,target) # measure rules



cleanEx()
nameEx("lookupRule")
### * lookupRule

flush(stderr()); flush(stdout())

### Name: lookupRule
### Title: internal
### Aliases: lookupRule
### Keywords: internal

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rules, strList) 
{
    ix <- grep(strList[1], rules[, "condition"])
    if (length(strList) >= 2) {
        for (i in 2:length(strList)) {
            ix2 <- grep(strList[i], rules[, "condition"])
            ix <- intersect(ix, ix2)
        }
    }
    if (length(ix) >= 1) 
        return(rules[ix, , drop = FALSE])
    if (length(ix) == 0) 
        return(NULL)
  }



cleanEx()
nameEx("measureRule")
### * measureRule

flush(stderr()); flush(stdout())

### Name: measureRule
### Title: internal
### Aliases: measureRule
### Keywords: internal

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ruleExec, X, target, pred = NULL) 
{
    len <- length(unlist(strsplit(ruleExec, split = " & ")))
    origRule <- ruleExec
    ruleExec <- paste("which(", ruleExec, ")")
    ixMatch <- eval(parse(text = ruleExec))
    if (length(ixMatch) == 0) {
        v <- c("-1", "-1", "-1", "", "")
        names(v) <- c("len", "freq", "err", "condition", "pred")
        return(v)
    }
    ys <- target[ixMatch]
    freq <- round(length(ys)/nrow(X), digits = 3)
    if (is.numeric(target)) {
        ysMost <- mean(ys)
        err <- sum((ysMost - ys)^2)/length(ys)
    }
    else {
        if (length(pred) > 0) {
            ysMost = pred
        }
        else {
            ysMost <- names(which.max(table(ys)))
        }
        conf <- round(table(ys)[ysMost]/sum(table(ys)), digits = 3)
        err <- 1 - conf
    }
    rule <- origRule
    v <- c(len, freq, err, rule, ysMost)
    names(v) <- c("len", "freq", "err", "condition", "pred")
    return(v)
  }



cleanEx()
nameEx("presentRules")
### * presentRules

flush(stderr()); flush(stdout())

### Name: presentRules
### Title: Present a learner using column names instead of X[i,]
### Aliases: presentRules
### Keywords: present

### ** Examples

 # See function "buildLearner"



cleanEx()
nameEx("pruneRule")
### * pruneRule

flush(stderr()); flush(stdout())

### Name: pruneRule
### Title: Prune irrevant variable-value pair from a rule condition
### Aliases: pruneRule
### Keywords: prune

### ** Examples

# see function "buildLearner"



cleanEx()
nameEx("pruneSingleRule")
### * pruneSingleRule

flush(stderr()); flush(stdout())

### Name: pruneSingleRule
### Title: internal
### Aliases: pruneSingleRule
### Keywords: internal

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rule, X, target, maxDecay, typeDecay) 
{
    newRuleMetric <- measureRule(rule["condition"], X, target)
    errOrig <- as.numeric(newRuleMetric["err"])
    ruleV <- unlist(strsplit(rule["condition"], split = " & "))
    pred <- rule["pred"]
    if (length(ruleV) == 1) 
        return(newRuleMetric)
    for (i in length(ruleV):1) {
        restRule <- ruleV[-i]
        restRule <- paste(restRule, collapse = " & ")
        metricTmp <- measureRule(restRule, X, target, pred)
        errNew <- as.numeric(metricTmp["err"])
        if (typeDecay == 1) {
            decay <- (errNew - errOrig)/max(errOrig, 1e-06)
        }
        else {
            decay <- (errNew - errOrig)
        }
        if (decay <= maxDecay) {
            ruleV <- ruleV[-i]
            newRuleMetric <- metricTmp
            if (length(ruleV) <= 1) 
                break
        }
    }
    return(newRuleMetric)
  }



cleanEx()
nameEx("rule2Table")
### * rule2Table

flush(stderr()); flush(stdout())

### Name: rule2Table
### Title: internal function
### Aliases: rule2Table
### Keywords: internal

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ruleExec, X, target) 
{
    I <- rep(0, nrow(X))
    ruleExec <- paste("which(", ruleExec, ")")
    ixMatch <- eval(parse(text = ruleExec))
    if (length(ixMatch) > 0) 
        I[ixMatch] <- 1
    names(I) = NULL
    return(I)
  }



cleanEx()
nameEx("ruleList2Exec")
### * ruleList2Exec

flush(stderr()); flush(stdout())

### Name: ruleList2Exec
### Title: internal
### Aliases: ruleList2Exec
### Keywords: internal

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, allRulesList) 
{
    typeX = getTypeX(X)
    ruleExec <- unique(t(sapply(allRulesList, singleRuleList2Exec, 
        typeX = typeX)))
    ruleExec <- t(ruleExec)
    colnames(ruleExec) <- "condition"
    return(ruleExec)
  }



cleanEx()
nameEx("selectRuleRRF")
### * selectRuleRRF

flush(stderr()); flush(stdout())

### Name: selectRuleRRF
### Title: select a set of relevant and non-redundant rules
### Aliases: selectRuleRRF
### Keywords: select

### ** Examples

 # See function "buildLearner:



cleanEx()
nameEx("singleRuleList2Exec")
### * singleRuleList2Exec

flush(stderr()); flush(stdout())

### Name: singleRuleList2Exec
### Title: internal
### Aliases: singleRuleList2Exec
### Keywords: internal

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ruleList, typeX) 
{
    ruleExec <- ""
    vars <- ls(ruleList)
    vars <- vars[order(as.numeric(vars))]
    for (i in 1:length(vars)) {
        if (typeX[as.numeric(vars[i])] == 2) {
            values <- paste("c(", paste(paste("'", ruleList[[vars[i]]], 
                "'", sep = ""), collapse = ","), ")", sep = "")
            tmp = paste("X[,", vars[i], "] %in% ", values, sep = "")
        }
        else {
            tmp = ruleList[[vars[i]]]
        }
        if (i == 1) 
            ruleExec <- paste(ruleExec, tmp, sep = "")
        if (i > 1) 
            ruleExec <- paste(ruleExec, " & ", tmp, sep = "")
    }
    return(c(ruleExec))
  }



cleanEx()
nameEx("sortRule")
### * sortRule

flush(stderr()); flush(stdout())

### Name: sortRule
### Title: internal
### Aliases: sortRule
### Keywords: internal

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (M, decreasing = TRUE) 
{
    qIx = order((1 - as.numeric(ruleMetric[, "err"])), as.numeric(ruleMetric[, 
        "freq"]), -as.numeric(ruleMetric[, "len"]), decreasing = decreasing)
    return(M[qIx, ])
  }



cleanEx()
nameEx("voteAllRules")
### * voteAllRules

flush(stderr()); flush(stdout())

### Name: voteAllRules
### Title: internal
### Aliases: voteAllRules
### Keywords: internal

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ruleMetric, X, type = "r", method = "median") 
{
    xVoteList = vector("list", nrow(X))
    predY <- rep("", nrow(X))
    for (i in 1:nrow(ruleMetric)) {
        ixMatch <- eval(parse(text = paste("which(", ruleMetric[i, 
            "condition"], ")")))
        if (length(ixMatch) == 0) 
            next
        for (ii in ixMatch) {
            xVoteList[[ii]] = c(xVoteList[[ii]], ruleMetric[i, 
                "pred"])
        }
    }
    for (i in 1:length(xVoteList)) {
        thisV <- xVoteList[[i]]
        if (length(thisV) == 0) 
            next
        if (type == "c") 
            predY[i] <- names(table(thisV)[which.max(table(thisV))])
        if (type == "r") {
            thisV = as.numeric(thisV)
            if (method == "median") {
                predY[i] <- median(thisV)
            }
            else {
                predY[i] <- mean(thisV)
            }
        }
    }
    if (type == "r") 
        predY <- as.numeric(predY)
    return(predY)
  }



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
