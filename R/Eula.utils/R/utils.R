
#' Get the max character length
#'
#' Get the max character length for a string vector
#'
#' @param str A string vector
#'
#' @export
maxNChar <- function(str) {
  max(nchar(as.character(str)))
}
