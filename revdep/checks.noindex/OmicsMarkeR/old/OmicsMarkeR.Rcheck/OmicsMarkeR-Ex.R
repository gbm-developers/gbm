pkgname <- "OmicsMarkeR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('OmicsMarkeR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("RPT")
### * RPT

flush(stderr()); flush(stdout())

### Name: RPT
### Title: Robustness-Performance Trade-Off
### Aliases: RPT

### ** Examples

# RPT demo
RPT(stability=0.85, performance=0.90, beta=1)



cleanEx()
nameEx("aggregation")
### * aggregation

flush(stderr()); flush(stdout())

### Name: aggregation
### Title: Feature Aggregation
### Aliases: aggregation

### ** Examples

# test data
ranks <- replicate(5, sample(seq(50), 50))
row.names(ranks) <- paste0("V", seq(50))

aggregation(ranks, "CLA")



cleanEx()
nameEx("canberra")
### * canberra

flush(stderr()); flush(stdout())

### Name: canberra
### Title: Canberra Distance
### Aliases: canberra

### ** Examples

# Canberra demo
v1 <- seq(10)
v2 <- sample(v1, 10)
canberra(v1, v2)

canberra_stability(v1, v2)



cleanEx()
nameEx("canberra_stability")
### * canberra_stability

flush(stderr()); flush(stdout())

### Name: canberra_stability
### Title: Canberra Stability
### Aliases: canberra_stability

### ** Examples

# Canberra demo
v1 <- seq(10)
v2 <- sample(v1, 10)
canberra(v1, v2)

canberra_stability(v1, v2)



cleanEx()
nameEx("create.corr.matrix")
### * create.corr.matrix

flush(stderr()); flush(stdout())

### Name: create.corr.matrix
### Title: Correlated Multivariate Data Generator
### Aliases: create.corr.matrix

### ** Examples

# Create Multivariate Matrices

# Random Multivariate Matrix

# 50 variables, 100 samples, 1 standard devation, 0.2 noise factor

rand.mat <- create.random.matrix(nvar = 50, 
                                 nsamp = 100, 
                                 st.dev = 1, 
                                 perturb = 0.2)


# Induce correlations in a numeric matrix

# Default settings
# minimum and maximum block sizes (min.block.size = 2, max.block.size = 5)
# default correlation purturbation (k=4)
# see ?create.corr.matrix for citation for methods

corr.mat <- create.corr.matrix(rand.mat)


# Induce Discriminatory Variables

# 10 discriminatory variables (D = 10)
# default discrimination level (l = 1.5)
# default number of groups (num.groups=2)
# default correlation purturbation (k = 4)

dat.discr <- create.discr.matrix(corr.mat, D=10)



cleanEx()
nameEx("create.discr.matrix")
### * create.discr.matrix

flush(stderr()); flush(stdout())

### Name: create.discr.matrix
### Title: Discriminatory Multivariate Data Generator
### Aliases: create.discr.matrix

### ** Examples

# Create Multivariate Matrices

# Random Multivariate Matrix

# 50 variables, 100 samples, 1 standard devation, 0.2 noise factor

rand.mat <- create.random.matrix(nvar = 50, 
                                 nsamp = 100, 
                                 st.dev = 1, 
                                 perturb = 0.2)


# Induce correlations in a numeric matrix

# Default settings
# minimum and maximum block sizes (min.block.size = 2, max.block.size = 5)
# default correlation purturbation (k=4)
# see ?create.corr.matrix for citation for methods

corr.mat <- create.corr.matrix(rand.mat)


# Induce Discriminatory Variables

# 10 discriminatory variables (D = 10)
# default discrimination level (l = 1.5)
# default number of groups (num.groups=2)
# default correlation purturbation (k = 4)

dat.discr <- create.discr.matrix(corr.mat, D=10)



cleanEx()
nameEx("create.random.matrix")
### * create.random.matrix

flush(stderr()); flush(stdout())

### Name: create.random.matrix
### Title: Random Multivariate Data Generator
### Aliases: create.random.matrix

### ** Examples

# Create Multivariate Matrices

# Random Multivariate Matrix

# 50 variables, 100 samples, 1 standard devation, 0.2 noise factor

rand.mat <- create.random.matrix(nvar = 50, 
                                 nsamp = 100, 
                                 st.dev = 1, 
                                 perturb = 0.2)


# Induce correlations in a numeric matrix

# Default settings
# minimum and maximum block sizes (min.block.size = 2, max.block.size = 5)
# default correlation purturbation (k=4)
# see ?create.corr.matrix for citation for methods

corr.mat <- create.corr.matrix(rand.mat)


# Induce Discriminatory Variables

# 10 discriminatory variables (D = 10)
# default discrimination level (l = 1.5)
# default number of groups (num.groups=2)
# default correlation purturbation (k = 4)

dat.discr <- create.discr.matrix(corr.mat, D=10)



cleanEx()
nameEx("denovo.grid")
### * denovo.grid

flush(stderr()); flush(stdout())

### Name: denovo.grid
### Title: Denovo Grid Generation
### Aliases: denovo.grid

### ** Examples


# random test data
set.seed(123)
dat.discr <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50, 
                             nsamp = 100, 
                             st.dev = 1, 
                             perturb = 0.2)),
    D = 10
)

df <- data.frame(dat.discr$discr.mat, .classes = dat.discr$classes)

# create tuning grid
denovo.grid(df, "gbm", 3)



cleanEx()
nameEx("feature.table")
### * feature.table

flush(stderr()); flush(stdout())

### Name: feature.table
### Title: Feature Consistency Table
### Aliases: feature.table

### ** Examples

dat.discr <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50, 
                             nsamp = 100, 
                             st.dev = 1, 
                             perturb = 0.2)),
    D = 10
)

vars <- dat.discr$discr.mat
groups <- dat.discr$classes

fits <- fs.stability(vars, 
                     groups, 
                     method = c("plsda", "rf"), 
                     f = 10, 
                     k = 3, 
                     k.folds = 10, 
                     verbose = 'none')

feature.table(fits, "plsda")



cleanEx()
nameEx("fit.only.model")
### * fit.only.model

flush(stderr()); flush(stdout())

### Name: fit.only.model
### Title: Fit Models without Feature Selection
### Aliases: fit.only.model

### ** Examples

dat.discr <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50, 
                             nsamp = 100, 
                             st.dev = 1, 
                             perturb = 0.2)),
    D = 10
)

vars <- dat.discr$discr.mat
groups <- dat.discr$classes

fit <- fit.only.model(X=vars, 
                      Y=groups, 
                      method="plsda", 
                      p = 0.9)



cleanEx()
nameEx("fs.ensembl.stability")
### * fs.ensembl.stability

flush(stderr()); flush(stdout())

### Name: fs.ensembl.stability
### Title: Ensemble Classification & Feature Selection
### Aliases: fs.ensembl.stability

### ** Examples

## Not run: 
##D fits <- fs.ensembl.stability(vars, 
##D groups, 
##D method = c("plsda", "rf"), 
##D f = 10,
##D k = 3, 
##D k.folds = 10, 
##D verbose = 'none')
## End(Not run)



cleanEx()
nameEx("fs.stability")
### * fs.stability

flush(stderr()); flush(stdout())

### Name: fs.stability
### Title: Classification & Feature Selection
### Aliases: fs.stability

### ** Examples

dat.discr <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50, 
                             nsamp = 100, 
                             st.dev = 1, 
                             perturb = 0.2)),
    D = 10
)

vars <- dat.discr$discr.mat
groups <- dat.discr$classes

fits <- fs.stability(vars, 
                     groups, 
                     method = c("plsda", "rf"), 
                     f = 10, 
                     k = 3, 
                     k.folds = 10, 
                     verbose = 'none')



cleanEx()
nameEx("jaccard")
### * jaccard

flush(stderr()); flush(stdout())

### Name: jaccard
### Title: Jaccard Index
### Aliases: jaccard

### ** Examples

# Jaccard demo
v1 <- paste("Metabolite", seq(10), sep="_")
v2 <- sample(v1, 10)
jaccard(v1, v2)



cleanEx()
nameEx("kuncheva")
### * kuncheva

flush(stderr()); flush(stdout())

### Name: kuncheva
### Title: Kuncheva's Index
### Aliases: kuncheva

### ** Examples

# Kuncheva demo
# Assuming 50 metabolites were measured
# But only 10 were found significant

# For demonstration purposes only!!!
some.numbers <- seq(20)

# Metabolites identified from one run
v1 <- paste("Metabolite", sample(some.numbers, 10), sep="_")
# Metabolites identifed from second run
v2 <- paste("Metabolite", sample(some.numbers, 10), sep="_")
kuncheva(v1, v2, 50)



cleanEx()
nameEx("modelList")
### * modelList

flush(stderr()); flush(stdout())

### Name: modelList
### Title: Model List
### Aliases: modelList

### ** Examples

modelList()



cleanEx()
nameEx("ochiai")
### * ochiai

flush(stderr()); flush(stdout())

### Name: ochiai
### Title: Ochiai's Index
### Aliases: ochiai

### ** Examples

# Ochiai demo
v1 <- paste("Metabolite", seq(10), sep="_")
v2 <- sample(v1, 10)
ochiai(v1, v2)



cleanEx()
nameEx("pairwise.model.stability")
### * pairwise.model.stability

flush(stderr()); flush(stdout())

### Name: pairwise.model.stability
### Title: Pairwise Model Stability Metrics
### Aliases: pairwise.model.stability

### ** Examples

# pairwise.model.stability demo
# For demonstration purposes only!!!
some.numbers <- seq(20)

# A list containing the metabolite matrices for each algorithm
# As an example, let's say we have the output from two different models
# such as plsda and random forest.
# matrix of Metabolites identified (e.g. 5 trials)
plsda <- 
    replicate(5, paste("Metabolite", sample(some.numbers, 10), sep="_"))
rf <-
    replicate(5, paste("Metabolite", sample(some.numbers, 10), sep="_"))

features <- list(plsda=plsda, rf=rf)

# nc may be omitted unless using kuncheva
pairwise.model.stability(features, "kuncheva", nc=20)



cleanEx()
nameEx("pairwise.stability")
### * pairwise.stability

flush(stderr()); flush(stdout())

### Name: pairwise.stability
### Title: Pairwise Stability Metrics
### Aliases: pairwise.stability

### ** Examples

# pairwise.stability demo

# For demonstration purposes only!!!
some.numbers <- seq(20)

# matrix of Metabolites identified (e.g. 5 trials)
features <- 
    replicate(5, paste("Metabolite", sample(some.numbers, 10), sep="_"))

# nc may be omitted unless using kuncheva
pairwise.stability(features, "jaccard")



cleanEx()
nameEx("params")
### * params

flush(stderr()); flush(stdout())

### Name: params
### Title: Model Parameters and Properties
### Aliases: params

### ** Examples

params("plsda")



cleanEx()
nameEx("performance.metrics")
### * performance.metrics

flush(stderr()); flush(stdout())

### Name: performance.metrics
### Title: Performance Metrics of fs.stability or fs.ensembl.stability
###   object
### Aliases: performance.metrics

### ** Examples

dat.discr <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50, 
                             nsamp = 100, 
                             st.dev = 1, 
                             perturb = 0.2)),
    D = 10
)

vars <- dat.discr$discr.mat
groups <- dat.discr$classes

fits <- fs.stability(vars, 
                     groups, 
                     method = c("plsda", "rf"), 
                     f = 10, 
                     k = 3, 
                     k.folds = 10, 
                     verbose = 'none')

performance.metrics(fits)



cleanEx()
nameEx("perm.class")
### * perm.class

flush(stderr()); flush(stdout())

### Name: perm.class
### Title: Monte Carlo Permutation of Model Performance
### Aliases: perm.class

### ** Examples

dat.discr <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50, 
                             nsamp = 100, 
                             st.dev = 1, 
                             perturb = 0.2)),
    D = 10
)

vars <- dat.discr$discr.mat
groups <- dat.discr$classes

fits <- fs.stability(vars, 
                     groups, 
                     method = c("plsda", "rf"), 
                     f = 10, 
                     k = 3, 
                     k.folds = 10, 
                     verbose = 'none')

perm.class(fits, vars, groups, "rf", k.folds=5,
           metric="Accuracy", nperm=10)



cleanEx()
nameEx("perm.features")
### * perm.features

flush(stderr()); flush(stdout())

### Name: perm.features
### Title: Feature Selection via Monte Carlo Permutation
### Aliases: perm.features

### ** Examples

dat.discr <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50, 
                             nsamp = 100, 
                             st.dev = 1, 
                             perturb = 0.2)),
    D = 10
)

vars <- dat.discr$discr.mat
groups <- dat.discr$classes

fits <- fs.stability(vars, 
                     groups, 
                     method = c("plsda", "rf"), 
                     f = 10, 
                     k = 3, 
                     k.folds = 10, 
                     verbose = 'none')

# permute variables/features
perm.features(fits, vars, groups, "rf",
              sig.level = .05, nperm = 10)



cleanEx()
nameEx("pof")
### * pof

flush(stderr()); flush(stdout())

### Name: pof
### Title: Percentage of Overlapping Features
### Aliases: pof

### ** Examples

# Percent-Overlapping Features demo
v1 <- paste("Metabolite", seq(10), sep="_")
v2 <- sample(v1, 10)
pof(v1, v2)



cleanEx()
nameEx("predictNewClasses")
### * predictNewClasses

flush(stderr()); flush(stdout())

### Name: predictNewClasses
### Title: Class Prediction
### Aliases: predictNewClasses

### ** Examples

dat.discr <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50, 
                             nsamp = 100, 
                             st.dev = 1, 
                             perturb = 0.2)),
    D = 10
)

vars <- dat.discr$discr.mat
groups <- dat.discr$classes

fits <- fs.stability(vars, 
                     groups, 
                     method = c("plsda", "rf"), 
                     f = 10, 
                     k = 3, 
                     k.folds = 10, 
                     verbose = 'none')

newdata <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50, 
                             nsamp = 100, 
                             st.dev = 1, 
                             perturb = 0.2)),
    D = 10
)$discr.mat

orig.df <- data.frame(vars, groups)

# see what the PLSDA predicts for the new data
# NOTE, newdata does not require a .classes column
predictNewClasses(fits, "plsda", orig.df, newdata)



cleanEx()
nameEx("sorensen")
### * sorensen

flush(stderr()); flush(stdout())

### Name: sorensen
### Title: Dice-Sorensen's Index
### Aliases: sorensen

### ** Examples

# Dice-Sorensen demo
v1 <- paste("Metabolite", seq(10), sep="_")
v2 <- sample(v1, 10)
sorensen(v1, v2)



cleanEx()
nameEx("spearman")
### * spearman

flush(stderr()); flush(stdout())

### Name: spearman
### Title: Spearman Rank Correlation Coefficient
### Aliases: spearman

### ** Examples

# Spearman demo
v1 <- seq(10)
v2 <- sample(v1, 10)
spearman(v1, v2)



cleanEx()
nameEx("svmrfeFeatureRanking")
### * svmrfeFeatureRanking

flush(stderr()); flush(stdout())

### Name: svmrfeFeatureRanking
### Title: SVM Recursive Feature Extraction (Binary)
### Aliases: svmrfeFeatureRanking

### ** Examples

dat.discr <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50, 
                             nsamp = 100, 
                             st.dev = 1, 
                             perturb = 0.2)),
    D = 10
)

vars <- dat.discr$discr.mat
groups <- dat.discr$classes

# binary class feature ranking
svmrfeFeatureRanking(x = vars,
                     y = groups, 
                     c = 0.1,
                     perc.rem = 10)



cleanEx()
nameEx("svmrfeFeatureRankingForMulticlass")
### * svmrfeFeatureRankingForMulticlass

flush(stderr()); flush(stdout())

### Name: svmrfeFeatureRankingForMulticlass
### Title: SVM Recursive Feature Extraction (Multiclass)
### Aliases: svmrfeFeatureRankingForMulticlass

### ** Examples

dat.discr <- create.discr.matrix(
    create.corr.matrix(
        create.random.matrix(nvar = 50, 
                             nsamp = 100, 
                             st.dev = 1, 
                             perturb = 0.2)),
    D = 10,
    num.groups=4
)

vars <- dat.discr$discr.mat
groups <- dat.discr$classes

# multiclass
svmrfeFeatureRankingForMulticlass(x = vars,
                                  y = groups, 
                                  c = 0.1,
                                  perc.rem = 10)



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
