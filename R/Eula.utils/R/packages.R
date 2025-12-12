#' @include verbose.R
#' @include check.R
NULL

#' @export
loadPackages <- function(packages, lib.loc = .libPaths(), verbose = TRUE) {
  verboseMsg("Checking R packages")
  checkPackages(packages = packages, required = TRUE, lib.loc = lib.loc)
  for (i in packages) {
    verboseMsg("Loading '", i, "'")
    suppressMessages(library(i, lib.loc = lib.loc, character.only = TRUE))
  }
  invisible(NULL)
}
