# gbm 2.3.0

#### Bug fixes

* Fixed predictions from models with categorical (factor) splits on platforms
  where plain `char` is unsigned (notably Linux on ARM64/aarch64). The
  categorical split-direction codes were stored in `char` and the value `-1`
  became `255`, silently routing all left-branch observations to the missing
  branch in `predict()` and `plot()`. Model training was unaffected; only
  predictions from the stored model object were wrong, and only on affected
  platforms.

* Fixed the multinomial deviance calculation, which was one iteration behind
  because it used cached class probabilities; also fixed the corresponding
  out-of-bag improvement calculation, which was double-counting updates.

* Fixed `gbm.more()` for `distribution = "multinomial"`: continuing training
  now correctly reconstructs the class ordering and data for both
  `keep.data = TRUE` and `keep.data = FALSE`, passes along the prior fit
  safely, restores the fit as an n x K matrix, and treats `n.trees` as
  boosting iterations rather than a raw class-tree count.

* Fixed two bugs in `gbm.more()` for `distribution = "coxph"`: observation
  weights could be silently replaced with the predictor matrix in certain
  cases, and `cRows`/`cCols` were not available yet when first needed.

* Fixed the out-of-bag improvement estimate for `distribution = "coxph"`,
  which ignored the current model fit; this corrupted
  `gbm.perf(method = "OOB")` for Cox models.

* Fixed observation weights not being reordered along with the data in
  `gbmCrossVal()`'s per-fold fitting, which corrupted weighted
  cross-validation.

* `InitF()` is now always called when fitting, including when extending an
  existing model via `gbm.more()`, so distributions can allocate any
  iteration-specific internal buffers they need.

* Fixed offset handling in several places: offsets are now correctly reordered
  alongside the data for `distribution = "pairwise"`; `distribution =
  "huberized"` now includes the offset in its terminal-node estimates;
  `distribution = "poisson"` now applies its prediction clamp when an offset
  is supplied; and supplying a real offset vector no longer errors in
  `gbm.fit()` and `gbm.more()` for Cox models.

* `permutation.test.gbm()` now works for all distributions, not just
  `"pairwise"`.

* Fixed undefined behavior in `distribution = "pairwise"` model fitting where
  internal buffers were reserved but not actually resized before being
  written to.

* Fixed a potential stack overflow in `plot.gbm()` for very deep trees.

* Fixed `distribution = "quantile"` to handle observation weights.

* Fixed Poisson terminal-node predictions when a node's denominator is zero.

* `gbm()` now raises an informative error when a supplied `weights` vector's
  length doesn't match the data, instead of failing obscurely later on.

* `plot.gbm()` no longer errors if the `viridis` package isn't installed;
  it now falls back to a built-in color palette.

* Added the missing y-axis label for `distribution = "huberized"` plots.

#### Behavior changes

* `gbm(distribution = "multinomial")` no longer emits the "ill-advised...
  currently broken" warning added in 2.1.6. Multinomial support has been
  fixed and is now tested. Note this removes a warning that some downstream
  packages' tests may assert on (see `pmml`, below).

#### Other improvements

* `distribution = "adaboost"` now uses the negative gradient as its working
  response, consistent with the other distributions; among other things this
  makes `var.monotone` constraints act in the intended direction for AdaBoost
  models.

* Switched the test suite from tinytest to testthat.

* Vignette corrections: Cox model formulas, the normalized discounted
  cumulative gain formula, and formulas added for the t-distribution,
  huberized hinge loss, and multinomial deviance; converted the vignette
  build from Sweave to R Markdown and added a pkgdown site.

# gbm 2.1.9

* Maintenance update to address new R standards

# gbm 2.1.8

* Removed experimental functions `shrink.gbm()` and `shrink.gbm.pred()`; the latter seemed broken anyway. Happy to accept a PR if anyone wants to fix them.


# gbm 2.1.7

* Fix `Non-file package-anchored link(s) in documentation...` warning.


# gbm 2.1.6

* Corrected the number of arguments for `gbm_shrink_gradient()` in `gbm-init.c` [(#50)](https://github.com/gbm-developers/gbm/issues/50). (Thanks to CRAN for highlighting the issue.)

* Removed unnecessary dependency on [gridExtra](https://cran.r-project.org/package=gridExtra).

* Switched to using `lapply()` instead of `parallel::parLapply()` whenever `n.cores = 1`.

* Calling `gbm()` with `distribution = "bernoulli"` will now throw an error whenever the response is non-numeric (e.g., 0/1 factors will throw an error instead of possibly crashing the session.) [(#6)](https://github.com/gbm-developers/gbm/issues/6). (Thanks to @mzoll.)

* Calling `gbm()` with `distribution = "multinomial"` now comes with a warning message; multinomial support has always been problematic and since this package is only being maintained for backwards compatibility, it likely will not be fixed unless someone makes a PR.

* Switched from [RUnit](https://cran.r-project.org/package=RUnit) to [tinytest](https://cran.r-project.org/package=tinytest) framework. The `test.gbm()`, `test.relative.influence()`, and `validate.gbm()` functions will remain for backwards compatability. This is just the start, as more tests will be added in the future [(#51)](https://github.com/gbm-developers/gbm/issues/51).


#### Bug fixes

* Fixed a long standing bug that could occur when using k-fold cross-validation with a response that's been transformed in the model formula [(#30)](https://github.com/gbm-developers/gbm/issues/30).

* Fixed a but that would crash the session when giving "bad" input for `n.trees` in the call to `predict.gbm()` [(#45)](https://github.com/gbm-developers/gbm/issues/45). (Thanks to @ngreifer.)

* Fixed a bug where calling `predict()` could throw an error in some cases when `n.trees` was not specified.


# gbm 2.1.5

* Fixed bug that occurred whenever `distribution` was a list (e.g., "pairwise" regression) [(#27)](https://github.com/gbm-developers/gbm/issues/27).

* Fixed a bug that occurred when making predictions on new data with different factor levels [(#28)](https://github.com/gbm-developers/gbm/issues/28).

* Fixed a bug that caused `relative.influence()` to give different values whenever `n.trees` was/wasn't given for multinomial distributions [(#31)](https://github.com/gbm-developers/gbm/issues/31).

* The `plot.it` argument of `gbm.perf()` is no longer ignored [(#34)](https://github.com/gbm-developers/gbm/issues/34). 

* Fixed an error that occurred in `gbm.perf()` whenever `oobag.curve = FALSE` and `overlay = FALSE`.


# gbm 2.1.4

* Switched from `CHANGES` to `NEWS` file.

* Updated links and maintainer field in `DESCRIPTION` file.

* Fixed bug caused by factors with unused levels
[(#5)](https://github.com/gbm-developers/gbm/issues/5).

* Fixed bug with axis labels in the `plot()` method for `"gbm"` objects [(#17)](https://github.com/gbm-developers/gbm/issues/17).

* The `plot()` method for `"gbm"` objects is now more consistent and always returns a `"trellis"` object [(#19)](https://github.com/gbm-developers/gbm/issues/19). Consequently, setting graphical parameters via `par` will no longer have an effect on the output from `plot.gbm()`.

* The `plot()` method for `"gbm"` objects gained five new arguments: `level.plot`, `contour`, `number`, `overlap`, and `col.regions`; see `?plot.gbm` for details.

* The default color palette for false color level plots in `plot.gbm()` has changed to the Matplotlib 'viridis' color map.

* Fixed a number of references and URLs.
