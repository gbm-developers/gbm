# Setup

## Platform

|setting  |value                        |
|:--------|:----------------------------|
|version  |R version 3.4.4 (2018-03-15) |
|system   |x86_64, darwin15.6.0         |
|ui       |RStudio (1.2.907)            |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|tz       |America/New_York             |
|date     |2018-09-14                   |

## Packages

|package |*  |version |date       |source                        |
|:-------|:--|:-------|:----------|:-----------------------------|
|gbm     |*  |2.1.4   |2018-09-14 |local (gbm-developers/gbm@NA) |
|RUnit   |   |0.4.32  |2018-05-18 |cran (@0.4.32)                |

# Check results

5 packages with problems

|package      |version | errors| warnings| notes|
|:------------|:-------|------:|--------:|-----:|
|EnsembleBase |1.0.2   |      1|        0|     0|
|lilikoi      |0.1.0   |      1|        0|     0|
|mboost       |2.9-1   |      1|        0|     0|
|SuperLearner |2.0-24  |      1|        0|     1|
|twang        |1.5     |      1|        0|     0|

## EnsembleBase (1.0.2)
Maintainer: Alireza S. Mahani <alireza.s.mahani@gmail.com>

1 error  | 0 warnings | 0 notes

```
checking whether package ‘EnsembleBase’ can be installed ... ERROR
Installation failed.
See ‘/Users/bgreenwell/Dropbox/devel/gbm/revdep/checks/EnsembleBase.Rcheck/00install.out’ for details.
```

## lilikoi (0.1.0)
Maintainer: Fadhl Alakwaa <falakwaa@hawaii.edu>  
Bug reports: https://github.com/lanagarmire/lilikoi/issues

1 error  | 0 warnings | 0 notes

```
checking whether package ‘lilikoi’ can be installed ... ERROR
Installation failed.
See ‘/Users/bgreenwell/Dropbox/devel/gbm/revdep/checks/lilikoi.Rcheck/00install.out’ for details.
```

## mboost (2.9-1)
Maintainer: Benjamin Hofner <benjamin.hofner@pei.de>  
Bug reports: https://github.com/boost-R/mboost/issues

1 error  | 0 warnings | 0 notes

```
checking whether package ‘mboost’ can be installed ... ERROR
Installation failed.
See ‘/Users/bgreenwell/Dropbox/devel/gbm/revdep/checks/mboost.Rcheck/00install.out’ for details.
```

## SuperLearner (2.0-24)
Maintainer: Eric Polley <polley.eric@mayo.edu>

1 error  | 0 warnings | 1 note 

```
checking tests ... ERROR
  Running ‘testthat.R’ [51s/51s]
Running the tests in ‘tests/testthat.R’ failed.
Last 13 lines of output:
  lstat       -0.019759 0.003639 -5.429 8.912e-08 ***
  
  ------------------------------------------------------------------- 
  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
  ---
  Residual standard error: 0.340537 on 492 degrees of freedom;
  observations: 506;  R^2: 0.519;  adjusted R^2: 0.506;
  F-statistic: 40.86 on 13 and 492 df;  p-value: 0.
  ══ testthat results  ══════════════════════════════════════════
  OK: 98 SKIPPED: 0 FAILED: 2
  1. Error: (unknown) (@test-bartMachine.R#2) 
  2. Error: (unknown) (@test-extraTrees.R#2) 
  
  Error: testthat unit tests failed
  Execution halted

checking package dependencies ... NOTE
Packages suggested but not available for checking: ‘genefilter’ ‘sva’
```

## twang (1.5)
Maintainer: Lane Burgette <burgette@rand.org>

1 error  | 0 warnings | 0 notes

```
checking whether package ‘twang’ can be installed ... ERROR
Installation failed.
See ‘/Users/bgreenwell/Dropbox/devel/gbm/revdep/checks/twang.Rcheck/00install.out’ for details.
```

