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

42 packages

|package             |version | errors| warnings| notes|
|:-------------------|:-------|------:|--------:|-----:|
|aurelius            |0.8.4   |      0|        0|     0|
|AzureML             |0.2.14  |      0|        0|     0|
|BigTSP              |1.0     |      0|        0|     3|
|BiodiversityR       |2.10-1  |      0|        0|     0|
|biomod2             |3.3-7   |      0|        0|     2|
|bst                 |0.3-15  |      0|        0|     0|
|bujar               |0.2-3   |      0|        0|     0|
|caretEnsemble       |2.0.0   |      0|        0|     0|
|crimelinkage        |0.0.4   |      0|        0|     1|
|DALEX               |0.2.4   |      0|        0|     0|
|dismo               |1.1-4   |      0|        0|     1|
|ecospat             |3.0     |      0|        0|     0|
|EnsembleBase        |1.0.2   |      1|        0|     0|
|featurefinder       |1.0     |      0|        0|     0|
|fscaret             |0.9.4.4 |      0|        0|     0|
|gbm2sas             |2.1     |      0|        0|     0|
|gbts                |1.2.0   |      0|        0|     0|
|horserule           |1.0.0   |      0|        0|     0|
|inTrees             |1.2     |      0|        0|     0|
|IPMRF               |1.2     |      0|        0|     0|
|lilikoi             |0.1.0   |      1|        0|     0|
|mboost              |2.9-1   |      1|        0|     0|
|mlr                 |2.13    |      0|        0|     2|
|mma                 |6.1-0   |      0|        0|     0|
|mvtboost            |0.5.0   |      0|        0|     0|
|opera               |1.0     |      0|        0|     0|
|pdp                 |0.7.0   |      0|        0|     0|
|personalized        |0.1.5   |      0|        0|     0|
|Plasmode            |0.1.0   |      0|        0|     0|
|plotmo              |3.5.0   |      0|        0|     0|
|pmml                |1.5.5   |      0|        0|     0|
|preprosim           |0.2.0   |      0|        0|     0|
|scorecardModelUtils |0.0.0.9 |      0|        0|     0|
|SDMPlay             |1.2     |      0|        0|     0|
|spm                 |1.1.1   |      0|        0|     0|
|SSDM                |0.2.4   |      0|        0|     0|
|subsemble           |0.0.9   |      0|        0|     3|
|SuperLearner        |2.0-24  |      1|        0|     1|
|tsensembler         |0.0.4   |      0|        0|     0|
|twang               |1.5     |      1|        0|     0|
|vip                 |0.1.0   |      0|        0|     0|
|WeightIt            |0.4.0   |      0|        0|     0|

## aurelius (0.8.4)
Maintainer: Steven Mortimer <reportmort@gmail.com>  
Bug reports: https://github.com/opendatagroup/hadrian/issues

0 errors | 0 warnings | 0 notes

## AzureML (0.2.14)
Maintainer: Rich Calaway <richcala@microsoft.com>  
Bug reports: https://github.com/RevolutionAnalytics/AzureML/issues

0 errors | 0 warnings | 0 notes

## BigTSP (1.0)
Maintainer: Xiaolin Yang <xyang@stat.cmu.edu>

0 errors | 0 warnings | 3 notes

```
checking dependencies in R code ... NOTE
Packages in Depends field not imported from:
  ‘gbm’ ‘glmnet’ ‘randomForest’ ‘tree’
  These packages need to be imported from (in the NAMESPACE file)
  for when this namespace is loaded but not attached.

checking R code for possible problems ... NOTE
LDCAEst: no visible global function definition for ‘glmnet’
cv.LDCA: no visible global function definition for ‘cv.glmnet’
predict.LDCA: no visible global function definition for ‘predict’
predict.cv.LDCA: no visible global function definition for ‘predict’
predict.tsp.gbm: no visible global function definition for ‘predict’
predict.tsp.randomForest: no visible global function definition for
  ‘predict’
predict.tsp.tree: no visible global function definition for ‘predict’
tsp.gbm: no visible global function definition for ‘gbm.fit’
tsp.randomForest: no visible global function definition for
  ‘randomForest’
tsp.tree: no visible global function definition for ‘tree.control’
tsp.tree: no visible global function definition for ‘tree’
Undefined global functions or variables:
  cv.glmnet gbm.fit glmnet predict randomForest tree tree.control
Consider adding
  importFrom("stats", "predict")
to your NAMESPACE file.

checking Rd line widths ... NOTE
Rd file 'predict.LDCA.Rd':
  \usage lines wider than 90 characters:
     predict(object, newx, s = NULL, type = c("link", "response", "coefficients", "nonzero", "class"), exact = FALSE, offset, ...)

Rd file 'predict.tsp.randomForest.Rd':
  \usage lines wider than 90 characters:
     predict(object, newdata, type = "response", norm.votes = TRUE, predict.all = FALSE, proximity = FALSE, nodes = FALSE, cutoff, ...)

Rd file 'predict.tsp.tree.Rd':
... 6 lines ...
     tsp.gbm(x, y, offset = NULL, misc = NULL, distribution = "bernoulli", w = NULL, var.monotone = NULL, n.trees = 100, interaction.depth = ... [TRUNCATED]

Rd file 'tsp.randomForest.Rd':
  \usage lines wider than 90 characters:
     tsp.randomForest(x, y = NULL, xtest = NULL, ytest = NULL, ntree = 500, type = "classification", mtry = if (!is.null(y) && !is.factor(y) ... [TRUNCATED]

Rd file 'tsp.tree.Rd':
  \usage lines wider than 90 characters:
     tsp.tree(X, response, control = tree.control(dim(X)[1], ...), method = "recursive.partition", split = c("deviance", "gini"), x = FALSE, ... [TRUNCATED]

These lines will be truncated in the PDF manual.
```

## BiodiversityR (2.10-1)
Maintainer: Roeland Kindt <R.KINDT@CGIAR.ORG>

0 errors | 0 warnings | 0 notes

## biomod2 (3.3-7)
Maintainer: Damien Georges <damien.georges2@gmail.com>  
Bug reports: <https://r-forge.r-project.org/R/?group_id=302>

0 errors | 0 warnings | 2 notes

```
checking DESCRIPTION meta-information ... NOTE
BugReports field is not a suitable URL but contains an email address
  which will be used as from R 3.4.0

checking dependencies in R code ... NOTE
Missing or unexported object: ‘gam::step.gam’
```

## bst (0.3-15)
Maintainer: Zhu Wang <zwang@connecticutchildrens.org>

0 errors | 0 warnings | 0 notes

## bujar (0.2-3)
Maintainer: Zhu Wang <zwang@connecticutchildrens.org>

0 errors | 0 warnings | 0 notes

## caretEnsemble (2.0.0)
Maintainer: Zachary A. Deane-Mayer <zach.mayer@gmail.com>  
Bug reports: https://github.com/zachmayer/caretEnsemble/issues

0 errors | 0 warnings | 0 notes

## crimelinkage (0.0.4)
Maintainer: Michael Porter <mporter@cba.ua.edu>  
Bug reports: <mporter@cba.ua.edu>

0 errors | 0 warnings | 1 note 

```
checking DESCRIPTION meta-information ... NOTE
BugReports field is not a suitable URL but contains an email address
  which will be used as from R 3.4.0
```

## DALEX (0.2.4)
Maintainer: Przemyslaw Biecek <przemyslaw.biecek@gmail.com>  
Bug reports: https://github.com/pbiecek/DALEX/issues

0 errors | 0 warnings | 0 notes

## dismo (1.1-4)
Maintainer: Robert J. Hijmans <r.hijmans@gmail.com>

0 errors | 0 warnings | 1 note 

```
checking compiled code ... NOTE
File ‘dismo/libs/dismo.so’:
  Found no calls to: ‘R_registerRoutines’, ‘R_useDynamicSymbols’

It is good practice to register native routines and to disable symbol
search.

See ‘Writing portable packages’ in the ‘Writing R Extensions’ manual.
```

## ecospat (3.0)
Maintainer: Olivier Broennimann <olivier.broennimann@unil.ch>  
Bug reports: https://github.com/ecospat/ecospat

0 errors | 0 warnings | 0 notes

## EnsembleBase (1.0.2)
Maintainer: Alireza S. Mahani <alireza.s.mahani@gmail.com>

1 error  | 0 warnings | 0 notes

```
checking whether package ‘EnsembleBase’ can be installed ... ERROR
Installation failed.
See ‘/Users/bgreenwell/Dropbox/devel/gbm/revdep/checks/EnsembleBase.Rcheck/00install.out’ for details.
```

## featurefinder (1.0)
Maintainer: Richard Davis <davisconsulting@gmail.com>

0 errors | 0 warnings | 0 notes

## fscaret (0.9.4.4)
Maintainer: Jakub Szlek <j.szlek@uj.edu.pl>

0 errors | 0 warnings | 0 notes

## gbm2sas (2.1)
Maintainer: John R. Dixon <gbm2sas@gmail.com>

0 errors | 0 warnings | 0 notes

## gbts (1.2.0)
Maintainer: Waley W. J. Liang <wliang10@gmail.com>

0 errors | 0 warnings | 0 notes

## horserule (1.0.0)
Maintainer: Malte Nalenz <malte.nlz@googlemail.com>

0 errors | 0 warnings | 0 notes

## inTrees (1.2)
Maintainer: Houtao Deng <softwaredeng@gmail.com>  
Bug reports: https://github.com/softwaredeng/inTrees/issues

0 errors | 0 warnings | 0 notes

## IPMRF (1.2)
Maintainer: Irene Epifanio <epifanio@uji.es>

0 errors | 0 warnings | 0 notes

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

## mlr (2.13)
Maintainer: Bernd Bischl <bernd_bischl@gmx.net>  
Bug reports: https://github.com/mlr-org/mlr/issues

0 errors | 0 warnings | 2 notes

```
checking installed package size ... NOTE
  installed size is  5.5Mb
  sub-directories of 1Mb or more:
    R      2.1Mb
    data   2.3Mb

checking dependencies in R code ... NOTE
Error in dyn.load(file, DLLpath = DLLpath, ...) : 
  unable to load shared object '/Library/Frameworks/R.framework/Versions/3.4/Resources/library/rgl/libs/rgl.so':
  `maximal number of DLLs reached...
```

## mma (6.1-0)
Maintainer: Qingzhao Yu <qyu@lsuhsc.edu>

0 errors | 0 warnings | 0 notes

## mvtboost (0.5.0)
Maintainer: Patrick Miller <patrick.mil10@gmail.com>  
Bug reports: https://github.com/patr1ckm/mvtboost/issues

0 errors | 0 warnings | 0 notes

## opera (1.0)
Maintainer: Pierre Gaillard <pierre@gaillard.me>  
Bug reports: https://github.com/dralliag/opera/issues

0 errors | 0 warnings | 0 notes

## pdp (0.7.0)
Maintainer: Brandon Greenwell <greenwell.brandon@gmail.com>  
Bug reports: https://github.com/bgreenwell/pdp/issues

0 errors | 0 warnings | 0 notes

## personalized (0.1.5)
Maintainer: Jared Huling <jaredhuling@gmail.com>  
Bug reports: https://github.com/jaredhuling/personalized/issues

0 errors | 0 warnings | 0 notes

## Plasmode (0.1.0)
Maintainer: Younathan Abdia <yabdia@bwh.harvard.edu>

0 errors | 0 warnings | 0 notes

## plotmo (3.5.0)
Maintainer: Stephen Milborrow <milbo@sonic.net>

0 errors | 0 warnings | 0 notes

## pmml (1.5.5)
Maintainer: Dmitriy Bolotov <dmitriy.bolotov@softwareag.com>

0 errors | 0 warnings | 0 notes

## preprosim (0.2.0)
Maintainer: Markus Vattulainen <markus.vattulainen@gmail.com>  
Bug reports: https://github.com/mvattulainen/preprosim/issues

0 errors | 0 warnings | 0 notes

## scorecardModelUtils (0.0.0.9)
Maintainer: Arya Poddar <aryapoddar290990@gmail.com>

0 errors | 0 warnings | 0 notes

## SDMPlay (1.2)
Maintainer: Guillaumot Charlene <charleneguillaumot21@gmail.com>

0 errors | 0 warnings | 0 notes

## spm (1.1.1)
Maintainer: Jin Li <jin.li@ga.gov.au>

0 errors | 0 warnings | 0 notes

## SSDM (0.2.4)
Maintainer: Sylvain Schmitt <sylvain.schmitt@agroparistech.fr>  
Bug reports: https://github.com/sylvainschmitt/SSDM/issues

0 errors | 0 warnings | 0 notes

## subsemble (0.0.9)
Maintainer: Erin LeDell <ledell@berkeley.edu>

0 errors | 0 warnings | 3 notes

```
checking dependencies in R code ... NOTE
'library' or 'require' call to ‘parallel’ in package code.
  Please use :: or requireNamespace() instead.
  See section 'Suggested packages' in the 'Writing R Extensions' manual.

checking S3 generic/method consistency ... NOTE
Found the following apparent S3 methods exported but not registered:
  predict.subsemble
See section ‘Registering S3 methods’ in the ‘Writing R Extensions’
manual.

checking R code for possible problems ... NOTE
predict.subsemble: no visible global function definition for ‘predict’
predict.subsemble: no visible binding for global variable ‘family’
subsemble: no visible global function definition for ‘gaussian’
subsemble: no visible global function definition for ‘detectCores’
subsemble : .make_Z_l: no visible global function definition for
  ‘parSapply’
subsemble : .make_Z_l: no visible global function definition for
  ‘makeCluster’
subsemble : .make_Z_l: no visible global function definition for
  ‘stopCluster’
subsemble: no visible global function definition for ‘parSapply’
subsemble: no visible global function definition for ‘makeCluster’
subsemble: no visible global function definition for ‘stopCluster’
subsemble: no visible global function definition for ‘predict’
Undefined global functions or variables:
  detectCores family gaussian makeCluster parSapply predict stopCluster
Consider adding
  importFrom("stats", "family", "gaussian", "predict")
to your NAMESPACE file.
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

## tsensembler (0.0.4)
Maintainer: Vitor Cerqueira <cerqueira.vitormanuel@gmail.com>

0 errors | 0 warnings | 0 notes

## twang (1.5)
Maintainer: Lane Burgette <burgette@rand.org>

1 error  | 0 warnings | 0 notes

```
checking whether package ‘twang’ can be installed ... ERROR
Installation failed.
See ‘/Users/bgreenwell/Dropbox/devel/gbm/revdep/checks/twang.Rcheck/00install.out’ for details.
```

## vip (0.1.0)
Maintainer: Brandon Greenwell <greenwell.brandon@gmail.com>  
Bug reports: https://github.com/koalaverse/vip/issues

0 errors | 0 warnings | 0 notes

## WeightIt (0.4.0)
Maintainer: Noah Greifer <noah.greifer@gmail.com>  
Bug reports: https://github.com/ngreifer/WeightIt/issues

0 errors | 0 warnings | 0 notes

