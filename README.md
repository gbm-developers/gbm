gbm
===

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/gbm)](https://cran.r-project.org/package=gbm)
[![Build
Status](https://travis-ci.org/gbm-developers/gbm.svg?branch=master)](https://travis-ci.org/gbm-developers/gbm)
[![Codecov test
coverage](https://codecov.io/gh/gbm-developers/gbm/branch/master/graph/badge.svg)](https://codecov.io/gh/gbm-developers/gbm?branch=master)
[![Downloads](http://cranlogs.r-pkg.org/badges/gbm)](http://cranlogs.r-pkg.org/badges/gbm)
[![Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/gbm)](http://cranlogs.r-pkg.org/badges/grand-total/gbm)

Overview
--------

The [gbm](https://cran.r-project.org/package=gbm) package, which stands
for **g**eneralized **b**oosted **m**odels, provides extensions to
Freund and Schapire’s AdaBoost algorithm and [Friedman’s gradient
boosting machine](http://projecteuclid.org/euclid.aos/1013203451). It
includes regression methods for least squares, absolute loss,
*t*-distribution loss, quantile regression, logistic, multinomial
logistic, Poisson, Cox proportional hazards partial likelihood, AdaBoost
exponential loss, Huberized hinge loss, and Learning to Rank measures
(i.e.,
[LambdaMart](https://www.microsoft.com/en-us/research/publication/from-ranknet-to-lambdarank-to-lambdamart-an-overview/)).

Installation
------------

``` r
# The easiest way to get gbm is to it install from CRAN:
install.packages("gbm")

# Alternatively, you can install the development version from GitHub:
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("gbm-developers/gbm")
```

Lifecycle
---------

[![lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)](https://www.tidyverse.org/lifecycle/#retired)

The [gbm](https://cran.r-project.org/package=gbm) package is retired and
no longer under active development. We will only make the necessary
changes to ensure that [gbm](https://cran.r-project.org/package=gbm)
remains on CRAN. For the most part, no new features will be added, and
only the most critical of bugs will be fixed.

This is a maintained version of
[gbm](https://cran.r-project.org/package=gbm) back compatible to CRAN
versions of [gbm](https://cran.r-project.org/package=gbm) 2.1.x. It
exists mainly for the purpose of reproducible research and data analyses
performed with the 2.1.x versions of
[gbm](https://cran.r-project.org/package=gbm). ~~For newer development,
and a more consistent API, try out the
[gbm3](https://github.com/gbm-developers/gbm3) package!~~
