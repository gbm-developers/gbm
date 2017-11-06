# NEWS for gbm package

### Changes for version 2.1.4
* Switched from `CHANGES` to `NEWS` file.
* Updated links and maintainer field in `DESCRIPTION` file.
* Fixed bug caused by factors with unused levels [(#5)](https://github.com/gbm-developers/gbm/issues/5).
* Fixed bug with axis labels in the `plot` method for `"gbm"` objects. [(#17)](https://github.com/gbm-developers/gbm/issues/17).
* The `plot` method for `"gbm"` objects is now more consistent and always returns a `"trellis"` object [(#19)](https://github.com/gbm-developers/gbm/issues/19). Consequently, setting graphical parameters via `par` will no longer have an effect on the output from `plot.gbm`.
* The `plot` method for `"gbm"` objects gained five new arguments: `level.plot`, `contour`, `number`, `overlap`, and `col.regions`; see `?plot.gbm` for details.
* The default color palette for false color level plots in `plot.gbm` has changed to the Matplotlib 'viridis' color map.
* Fixed a number of references and URLs.
