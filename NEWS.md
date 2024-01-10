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
