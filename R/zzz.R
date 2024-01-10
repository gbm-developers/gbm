#' @keywords internal
.onAttach <- function(lib, pkg) {
  vers <- utils::packageVersion("gbm")
  packageStartupMessage("Loaded gbm ", vers)
  packageStartupMessage("This version of gbm is no longer under development. Consider transitioning to gbm3, https://github.com/gbm-developers/gbm3")
}
