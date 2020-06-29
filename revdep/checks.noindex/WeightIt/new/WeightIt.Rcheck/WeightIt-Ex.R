pkgname <- "WeightIt"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('WeightIt')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ESS")
### * ESS

flush(stderr()); flush(stdout())

### Name: ESS
### Title: Compute effective sample size of weighted sample
### Aliases: ESS

### ** Examples

library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ps", estimand = "ATE"))
summary(W1)
ESS(W1$weights[W1$treat == 0])
ESS(W1$weights[W1$treat == 1])



cleanEx()
nameEx("as.weightit")
### * as.weightit

flush(stderr()); flush(stdout())

### Name: as.weightit
### Title: Create a 'weightit' object manually
### Aliases: as.weightit as.weightitMSM as.weightit.default
###   as.weightitMSM.default

### ** Examples

treat <- rbinom(500, 1, .3)
weights <- rchisq(500, df = 2)
W <- as.weightit(weights= weights, treat = treat,
                 estimand = "ATE")
summary(W)



cleanEx()
nameEx("get_w_from_ps")
### * get_w_from_ps

flush(stderr()); flush(stdout())

### Name: get_w_from_ps
### Title: Get weights from propensity scores corresponding to different
###   estimands
### Aliases: get_w_from_ps

### ** Examples

library("cobalt")
data("lalonde", package = "cobalt")

ps.fit <- glm(treat ~ age + educ + race + married +
                nodegree + re74 + re75, data = lalonde,
              family = binomial)
ps <- ps.fit$fitted.values

w1 <- get_w_from_ps(ps, treat = lalonde$treat,
                    estimand = "ATT")

treatAB <- factor(ifelse(lalonde$treat == 1, "A", "B"))
w2 <- get_w_from_ps(ps, treat = treatAB,
                    estimand = "ATT", focal = "A")
all.equal(w1, w2)
w3 <- get_w_from_ps(ps, treat = treatAB,
                    estimand = "ATT", treated = "A")
all.equal(w1, w3)

#Using MMWS
w4 <- get_w_from_ps(ps, treat = lalonde$treat,
                    estimand = "ATE", subclass = 20,
                    stabilize = TRUE)

#A multi-category example using GBM predicted probabilities
library(gbm)
T3 <- factor(sample(c("A", "B", "C"), nrow(lalonde), replace = TRUE))

gbm.fit <- gbm(T3 ~ age + educ + race + married +
                 nodegree + re74 + re75, data = lalonde,
               distribution = "multinomial", n.trees = 200,
               interaction.depth = 3)
ps.multi <- drop(predict(gbm.fit, type = "response",
                         n.trees = 200))
w <- get_w_from_ps(ps.multi, T3, estimand = "ATE")




cleanEx()
nameEx("make_full_rank")
### * make_full_rank

flush(stderr()); flush(stdout())

### Name: make_full_rank
### Title: Make a matrix full rank
### Aliases: make_full_rank

### ** Examples

set.seed(1000)
c1 <- rbinom(10, 1, .4)
c2 <- 1-c1
c3 <- rnorm(10)
c4 <- 10*c3
mat <- data.frame(c1, c2, c3, c4)

make_full_rank(mat) #leaves c2 and c4

make_full_rank(mat, with.intercept = FALSE) #leaves c1, c2, and c4



cleanEx()
nameEx("method-cbps")
### * method-cbps

flush(stderr()); flush(stdout())

### Name: method_cbps
### Title: Covariate Balancing Propensity Score Weighting
### Aliases: method_cbps

### ** Examples

library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "cbps", estimand = "ATT"))
summary(W1)
bal.tab(W1)

## Not run: 
##D #Balancing covariates with respect to race (multinomial)
##D (W2 <- weightit(race ~ age + educ + married +
##D                   nodegree + re74, data = lalonde,
##D                 method = "cbps", estimand = "ATE"))
##D summary(W2)
##D bal.tab(W2)
## End(Not run)

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "cbps", over = FALSE))
summary(W3)
bal.tab(W3)



cleanEx()
nameEx("method-ebal")
### * method-ebal

flush(stderr()); flush(stdout())

### Name: method_ebal
### Title: Entropy Balancing
### Aliases: method_ebal

### ** Examples

library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal", estimand = "ATT"))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal", estimand = "ATE",
                standardize = TRUE))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal"))
summary(W3)
bal.tab(W3)




cleanEx()
nameEx("method-ebcw")
### * method-ebcw

flush(stderr()); flush(stdout())

### Name: method_ebcw
### Title: Empirical Balancing Calibration Weighting
### Aliases: method_ebcw

### ** Examples

library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebcw", estimand = "ATT"))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebcw", estimand = "ATE"))
summary(W2)
bal.tab(W2)




cleanEx()
nameEx("method-gbm")
### * method-gbm

flush(stderr()); flush(stdout())

### Name: method_gbm
### Title: Propensity Score Weighting Using Generalized Boosted Models
### Aliases: method_gbm

### ** Examples


library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "gbm", estimand = "ATE",
                stop.method = "es.max"))
summary(W1)
bal.tab(W1)

## Not run: 
##D #Balancing covariates with respect to race (multinomial)
##D (W2 <- weightit(race ~ age + educ + married +
##D                   nodegree + re74, data = lalonde,
##D                 method = "gbm", estimand = "ATT",
##D                 focal = "hispan", stop.method = "ks.mean"))
##D summary(W2)
##D bal.tab(W2)
##D 
##D #Balancing covariates with respect to re75 (continuous)
##D (W3 <- weightit(re75 ~ age + educ + married +
##D                   nodegree + re74, data = lalonde,
##D                 method = "gbm", use.kernel = TRUE,
##D                 stop.method = "p.rms", trim.at = .97))
##D summary(W3)
##D bal.tab(W3)
##D 
##D #Using a t(3) density and illustrating the search for
##D #more trees.
##D W4a <- weightit(re75 ~ age + educ + married +
##D                   nodegree + re74, data = lalonde,
##D                 method = "gbm", density = "dt_3",
##D                 stop.method = "p.max",
##D                 n.trees = 10000)
##D 
##D W4a$info$best.tree #10000; optimum hasn't been found
##D plot(W4a$info$tree.val) #decreasing at right edge
##D 
##D W4b <- weightit(re75 ~ age + educ + married +
##D                   nodegree + re74, data = lalonde,
##D                 method = "gbm", density = "dt_3",
##D                 stop.method = "p.max",
##D                 start.tree = 10000,
##D                 n.trees = 20000)
##D 
##D W4b$info$best.tree #13417; optimum has been found
##D plot(W4b$info$tree.val) #increasing at right edge
##D 
##D bal.tab(W4b)
##D 
## End(Not run)



cleanEx()
nameEx("method-npcbps")
### * method-npcbps

flush(stderr()); flush(stdout())

### Name: method_npcbps
### Title: Nonparametric Covariate Balancing Propensity Score Weighting
### Aliases: method_npcbps

### ** Examples

# Examples take a long time to run
## Not run: 
##D library("cobalt")
##D data("lalonde", package = "cobalt")
##D 
##D #Balancing covariates between treatment groups (binary)
##D (W1 <- weightit(treat ~ age + educ + married +
##D                   nodegree + re74, data = lalonde,
##D                 method = "npcbps", estimand = "ATE"))
##D summary(W1)
##D bal.tab(W1)
##D 
##D #Balancing covariates with respect to race (multinomial)
##D (W2 <- weightit(race ~ age + educ + married +
##D                   nodegree + re74, data = lalonde,
##D                 method = "npcbps", estimand = "ATE"))
##D summary(W2)
##D bal.tab(W2)
## End(Not run)



cleanEx()
nameEx("method-optweight")
### * method-optweight

flush(stderr()); flush(stdout())

### Name: method_optweight
### Title: Optimization-Based Weighting
### Aliases: method_optweight

### ** Examples

library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "optweight", estimand = "ATT",
                tols = 0))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "optweight", estimand = "ATE",
                tols = .01))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "optweight", tols = .05))
summary(W3)
bal.tab(W3)



cleanEx()
nameEx("method-ps")
### * method-ps

flush(stderr()); flush(stdout())

### Name: method_ps
### Title: Propensity Score Weighting Using Generalized Linear Models
### Aliases: method_ps

### ** Examples

library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ps", estimand = "ATT",
                link = "probit"))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ps", estimand = "ATE"))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ps", use.kernel = TRUE))
summary(W3)
bal.tab(W3)



cleanEx()
nameEx("method-super")
### * method-super

flush(stderr()); flush(stdout())

### Name: method_super
### Title: Propensity Score Weighting Using SuperLearner
### Aliases: method_super

### ** Examples

library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "super", estimand = "ATT",
                SL.library = c("SL.glm", "SL.stepAIC",
                               "SL.glm.interaction")))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "super", estimand = "ATE",
                SL.library = c("SL.glm", "SL.stepAIC",
                               "SL.glm.interaction")))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re75 (continuous)
#assuming t(8) conditional density for treatment
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "super", density = "dt_8",
                SL.library = c("SL.glm", "SL.ridge",
                               "SL.glm.interaction")))
summary(W3)
bal.tab(W3)

#Balancing covariates between treatment groups (binary)
# using balance SuperLearner to minimize the average
# KS statistic
(W4 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "super", estimand = "ATT",
                SL.library = c("SL.glm", "SL.stepAIC",
                               "SL.lda"),
                SL.method = "method.balance",
                stop.method = "ks.mean"))
summary(W4)
bal.tab(W4)



cleanEx()
nameEx("method-twang")
### * method-twang

flush(stderr()); flush(stdout())

### Name: method_twang
### Title: Propensity Score Weighting Using Generalized Boosted Models with
###   TWANG
### Aliases: method_twang

### ** Examples

# Examples take a while to run
## Not run: 
##D library("cobalt")
##D data("lalonde", package = "cobalt")
##D 
##D #Balancing covariates between treatment groups (binary)
##D (W1 <- weightit(treat ~ age + educ + married +
##D                   nodegree + re74, data = lalonde,
##D                 method = "twang", estimand = "ATT",
##D                 stop.method = "es.max"))
##D summary(W1)
##D bal.tab(W1)
##D 
##D #Balancing covariates with respect to race (multinomial)
##D (W2 <- weightit(race ~ age + educ + married +
##D                   nodegree + re74, data = lalonde,
##D                 method = "twang", estimand = "ATE",
##D                 stop.method = "ks.mean"))
##D summary(W2)
##D bal.tab(W2)
##D 
##D #Balancing covariates with respect to re75 (continuous)
##D (W3 <- weightit(re75 ~ age + educ + married +
##D                   nodegree + re74, data = lalonde,
##D                 method = "twang", use.kernel = TRUE,
##D                 stop.method = "p.max"))
##D summary(W3)
##D bal.tab(W3)
## End(Not run)



cleanEx()
nameEx("method-user")
### * method-user

flush(stderr()); flush(stdout())

### Name: method_user
### Title: User-Defined Functions for Estimating Weights
### Aliases: method_user

### ** Examples


library("cobalt")
data("lalonde", package = "cobalt")

#A user-defined version of method = "ps"
my.ps <- function(treat, covs, estimand, focal = NULL) {
  covs <- make_full_rank(covs)
  d <- data.frame(treat, covs)
  f <- formula(d)
  ps <- glm(f, data = d, family = "binomial")$fitted
  w <- get_w_from_ps(ps, treat = treat, estimand = estimand,
                     focal = focal)

  return(list(w = w, ps = ps))
}

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = my.ps, estimand = "ATT"))
summary(W1)
bal.tab(W1)

#Balancing covariates for longitudinal treatments
# my.ps is used at each time point.
library("twang")
data("iptwExWide", package = "twang")
(W2 <- weightitMSM(list(tx1 ~ age + gender + use0,
                        tx2 ~ tx1 + use1 + age + gender + use0,
                        tx3 ~ tx2 + use2 + tx1 + use1 + age + gender + use0),
                   data = iptwExWide,
                   method = my.ps))
summary(W2)
bal.tab(W2)

# Kernel balancing using the KBAL package, available
# using devtools::install_github("chadhazlett/KBAL").
# Only the ATT and ATC are available. Use 'kbal.method'
# instead of 'method' in weightit() to choose between
# "ebal" and "el".

## Not run: 
##D kbal.fun <- function(treat, covs, estimand, focal, ...) {
##D     args <- list(...)
##D     if (is_not_null(focal))
##D         treat <- as.numeric(treat == focal)
##D     else if (estimand != "ATT")
##D         stop("estimand must be 'ATT' or 'ATC'.", call. = FALSE)
##D     if ("kbal.method" %in% names(args)) {
##D         names(args)[names(args) == "kbal.method"] <- "method"
##D     }
##D     args[names(args) %nin% setdiff(names(formals(KBAL::kbal)),
##D         c("X", "D"))] <- NULL
##D     k.out <- do.call(KBAL::kbal, c(list(X = covs, D = treat),
##D         args))
##D     w <- k.out$w
##D     return(list(w = w))
##D }
##D 
##D (Wk <- weightit(treat ~ age + educ + married +
##D                 nodegree + re74, data = lalonde,
##D                 method = kbal.fun, estimand = "ATT",
##D                 kbal.method = "ebal"))
##D summary(Wk)
##D bal.tab(Wk, disp.ks = TRUE)
## End(Not run)



cleanEx()
nameEx("ps.cont")
### * ps.cont

flush(stderr()); flush(stdout())

### Name: ps.cont
### Title: Generalized Propensity Score Estimation using GBM
### Aliases: ps.cont gbm.cont summary.ps.cont plot.ps.cont boxplot.ps.cont

### ** Examples

# Examples take a long time
## Not run: 
##D library("cobalt")
##D data("lalonde", package = "cobalt")
##D 
##D #Balancing covariates with respect to re75
##D psc.out <- ps.cont(re75 ~ age + educ + married +
##D                 nodegree + race + re74, data = lalonde,
##D                 stop.method = c("p.mean", "p.max"),
##D                 use.optimize = 2)
##D summary(psc.out)
##D twang::bal.table(psc.out) #twang's bal.table
## End(Not run)



cleanEx()
nameEx("sbps")
### * sbps

flush(stderr()); flush(stdout())

### Name: sbps
### Title: Subgroup Balancing Propensity Score
### Aliases: sbps SBPS summary.weightit.sbps print.weightit.sbps
###   print.summary.weightit.sbps

### ** Examples

library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups within races
(W1 <- weightit(treat ~ age + educ + married +
                nodegree + race + re74, data = lalonde,
                method = "ps", estimand = "ATT"))

(W2 <- weightit(treat ~ age + educ + married +
                nodegree + race + re74, data = lalonde,
                method = "ps", estimand = "ATT",
                by = "race"))
S <- sbps(W1, W2)
print(S)
summary(S)
bal.tab(S, cluster = "race")

#Could also have run
#  sbps(W1, moderator = "race")

S_ <- sbps(W1, W2, smooth = TRUE)
print(S_)
summary(S_)
bal.tab(S_, cluster = "race")




cleanEx()
nameEx("summary.weightit")
### * summary.weightit

flush(stderr()); flush(stdout())

### Name: summary.weightit
### Title: Print and Summarize Output
### Aliases: print.summary.weightit plot.summary.weightit summary.weightit
###   print.summary.weightitMSM summary.weightitMSM

### ** Examples

# See example at ?weightit or ?weightitMSM



cleanEx()
nameEx("trim")
### * trim

flush(stderr()); flush(stdout())

### Name: trim
### Title: Trim Large Weights
### Aliases: trim trim.weightit trim.numeric

### ** Examples

library("cobalt")
data("lalonde", package = "cobalt")

(W <- weightit(treat ~ age + educ + married +
                 nodegree + re74, data = lalonde,
               method = "ps", estimand = "ATT"))
summary(W)

#Trimming the top and bottom 5 weights
trim(W, at = 5, lower = TRUE)

#Trimming at 90th percentile
(W.trim <- trim(W, at = .9))

summary(W.trim)
#Note that only the control weights were trimmed

#Trimming a numeric vector of weights
weights <- cobalt::get.w(W)

all.equal(trim(weights, at = .9, treat = lalonde$treat),
          W.trim$weights)

#Using made up data and as.weightit()
treat <- rbinom(500, 1, .3)
weights <- rchisq(500, df = 2)
W <- as.weightit(weights= weights, treat = treat,
                 estimand = "ATE")
summary(W)
summary(trim(W, at = .95))



cleanEx()
nameEx("weightit")
### * weightit

flush(stderr()); flush(stdout())

### Name: weightit
### Title: Generate Balancing Weights
### Aliases: weightit WeightIt print.weightit

### ** Examples

library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ps", estimand = "ATT"))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multi-category)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal", estimand = "ATE"))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "cbps", over = FALSE))
summary(W3)
bal.tab(W3)



cleanEx()
nameEx("weightitMSM")
### * weightitMSM

flush(stderr()); flush(stdout())

### Name: weightitMSM
### Title: Generate Balancing Weights
### Aliases: weightitMSM print.weightitMSM

### ** Examples

library("twang")
# library("cobalt")

data("iptwExWide", package = "twang")
(W <- weightitMSM(list(tx1 ~ age + gender + use0,
                      tx2 ~ tx1 + use1 + age + gender + use0,
                      tx3 ~ tx2 + use2 + tx1 + use1 + age + gender + use0),
                 data = iptwExWide,
                 method = "ps"))
summary(W)
# bal.tab(W)

##Going from long to wide data
data("iptwExLong", package = "twang")
wide_data <- reshape(iptwExLong$covariates,    #long data
                     timevar = "time",         #time variable
                     v.names = c("use", "tx"), #time-varying
                     idvar = "ID",             #time-stable
                     direction = "wide",
                     sep = "")

(W2 <- weightitMSM(list(tx1 ~ age + gender + use1,
                      tx2 ~ tx1 + use2 + age + gender + use1,
                      tx3 ~ tx2 + use3 + tx1 + use2 + age +
                          gender + use1),
                 data = wide_data,
                 method = "ps"))
summary(W2)

all.equal(W$weights, W2$weights)



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
