#' @include verbose.R
#' @include check.R
NULL

#' @export
loadPackages <- function(packages, lib.loc = .libPaths()) {
  pipeMsg("Checking R packages")
  checkPackages(packages = packages, required = TRUE, lib.loc = lib.loc)
  for (i in packages) {
    pipeMsg("Loading '", i, "'")
    suppressMessages(library(i, lib.loc = lib.loc, character.only = TRUE))
  }
  invisible(NULL)
}
