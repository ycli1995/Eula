
#' @keywords internal
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  require("showtext", quietly = TRUE)
  showtext_auto()
}
