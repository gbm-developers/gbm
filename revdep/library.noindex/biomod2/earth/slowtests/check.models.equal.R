# check.models.equal.R

fuzzy.equal <- function(x, y, fuzz=1e-8, msg="")
{
    if (any(y == 0))
        diff <- x - y
    else
        diff <- (x - y) / y
    eq <- all(abs(diff) < fuzz)
    if(!eq) {
        cat("First argument for", msg, "\n")
        print(head(x))
        cat("Relative difference\n")
        print(head(abs((x - y) / y)))
    }
    eq
}

check.fuzzy.equal <- function(x, y, fuzz=1e-8, msg="", verbose=FALSE)
{
    # fuzzy wuzzy was a bear, fuzzy wuzzy had no hair,
    # so fuzzy wuzzy wasn't fuzzy, wuzzy?
    # ah, but wuzzy bare?

    if (any(y == 0))
        diff <- x - y
    else
        diff <- (x - y) / y
    if (any(abs(diff) > fuzz)) {
        cat(msg, "\n1st matrix:\n", sep="")
        print(x)
        cat("2nd matrix:\n")
        print(y)
        cat("diff:\n")
        print(diff)
        stop("check.fuzzy.equal failed for ", msg, call.=FALSE)
    }
    if (verbose)
        cat(msg, "OK\n")
}

check.same <- function(m1, m2, msg="", allow.different.names=FALSE, fuzz=0)
{
    if (!identical(m1, m2)) {
        stop.with.msg = TRUE
        if (!is.null(dim(m1)) && !is.null(dim(m2))) {
            # check if it is the column names
            m1a <- m1
            m2a <- m2
            colnames(m1a) <- NULL
            colnames(m2a) <- NULL
            if (identical(m1a, m2a)) {
                cat("\nm1", msg, "\n"); print(colnames(m1)); cat("m2", msg, "\n"); print(colnames(m2)); cat("\n")
                if (allow.different.names) {
                    warning(msg, " has different column names but is otherwise identical, see above messages", call.=FALSE)
                    stop.with.msg = FALSE
                } else
                    stop(msg, " has different column names but is otherwise identical, see above messages", call.=FALSE)
            }
            # check if it is the row names
            m1a <- m1
            m2a <- m2
            rownames(m1a) <- NULL
            rownames(m2a) <- NULL
            if (identical(m1a, m2a)) {
                cat("\nm1", msg, "\n"); print(head(m1)); cat("\nm2", msg, "\n"); print(head(m2)); cat("\n")
                if (allow.different.names) {
                    warning(msg, " has different row names but is otherwise identical, see above messages", call.=FALSE)
                    stop.with.msg = FALSE
                } else
                    stop(msg, " has different row names but is otherwise identical, see above messages", call.=FALSE)
            }
        }
        if(fuzz != 0) {
            same <- fuzzy.equal(m1, m2, fuzz=fuzz, msg=msg)
            stop.with.msg = FALSE
        }
        if (stop.with.msg) {
            cat("\nm1", msg, "\n"); print(m1);
            # cat("\nm2", msg, "\n"); print(m2);
            cat("\ndifference m1-m2", msg, "\n"); print(m1-m2);
            cat("\n")
            stop(msg, " don't match, see above messages (fuzz=", fuzz, ")", call.=FALSE)
        }
    }
}

check.models.equal <- function(m1, m2, msg="", check.subsets=TRUE, allow.different.names=FALSE)
{
    m1$call <- NULL
    m2$call <- NULL
    m1$trace <- NULL
    m2$trace <- NULL
    if (identical(m1, m2))
        cat("check.models.equal identical: ", msg, "\n", sep="")
    else {
        cat("check.models.equal not identical: ", msg, " ", sep="");
        # cat("m1\n"); print(summary(m1)); cat("m2\n"); print(summary(m2)); cat("\n")
        # TODO why do we need a fuzz here and below?
        check.same(m1$bx, m2$bx, "bx", allow.different.names=allow.different.names, fuzz=1e-14)
        check.same(m1$coefficients, m2$coefficients, "coefficients", allow.different.names=allow.different.names, fuzz=1e-14)
        check.same(m1$dirs, m2$dirs, "dirs", allow.different.names=allow.different.names)
        check.same(m1$cuts, m2$cuts, "cuts", allow.different.names=allow.different.names)
        check.same(m1$residuals, m2$residuals, "residuals", fuzz=1e-14)
        check.same(m1$selected.terms, m2$selected.terms, "selected.terms")
        if (check.subsets) {
            # leaps and xtx pruning can give different prune.terms, so skip test
            check.same(m1$prune.terms, m2$prune.terms, "prune.terms")
            check.same(m1$rss.per.response, m2$rss.per.response, "rss.per.response", fuzz=1e-14)
            check.same(m1$rsq.per.response, m2$rsq.per.response, "rsq.per.response")
            check.same(m1$gcv.per.response, m2$gcv.per.response, "gcv.per.response")
            check.same(m1$grsq.per.response, m2$grsq.per.response, "grsq.per.response")
            check.same(m1$rss.per.subset, m2$rss.per.subset, "rss.per.subset", fuzz=1e-14)
            check.same(m1$gcv.per.subset, m2$gcv.per.subset, "gcv.per.subset", fuzz=1e-14)
        }
        if (!fuzzy.equal(m1$rss, m2$rss, msg=msg))
            stop("different rss")
        if (!fuzzy.equal(m1$rsq, m2$rsq, msg=msg))
            stop("different rsq")
        if (!fuzzy.equal(m1$gcv, m2$gcv, msg=msg))
            stop("different gcv")
        if (!fuzzy.equal(m1$grsq, m2$grsq, msg=msg))
            stop("different grsq")
        if (m1$rsq != m2$rsq)
            cat("(m1$rsq ", m1$rsq, " != m1$rsq ", m1$rsq, ") ", sep="")
        cat("[but within numerical tolerances]")
    }
    cat("\n")
}
