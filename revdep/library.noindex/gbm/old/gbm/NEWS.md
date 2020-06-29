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

* The `plot()` method for `"gbm"` objects is now more consistent and always returns a `"trellis"` object [(#19)](https://github.com/gbm-developers/gbm/issues/19). Consequently, setting graphical parameters via `par` will no longer have an effect on the output from `plot.gbm`.

* The `plot()` method for `"gbm"` objects gained five new arguments: `level.plot`, `contour`, `number`, `overlap`, and `col.regions`; see `?plot.gbm` for details.

* The default color palette for false color level plots in `plot.gbm()` has changed to the Matplotlib 'viridis' color map.

* Fixed a number of references and URLs.
