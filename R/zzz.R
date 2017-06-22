#' @keywords internal
.onAttach <- function(lib, pkg) {
  vers <- utils::packageVersion("gbm")
  packageStartupMessage(paste("Loaded gbm", vers))
}
