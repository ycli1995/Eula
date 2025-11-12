
#' @keywords internal
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  suppressMessages(require("extrafont", quietly = TRUE))
  suppressMessages(require("extrafontdb", quietly = TRUE))
  suppressMessages(require("Rttf2pt1", quietly = TRUE))
  options(device = cairo_pdf)
}
