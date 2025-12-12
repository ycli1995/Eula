#' @include verbose.R
NULL

#' @export
checkColumns <- function(x, cols) {
  if (length(cols) == 0) {
    stop("No column to check.")
  }
  if (!all(cols %in% colnames(x))) {
    cols <- paste(cols, collapse = ", ")
    stop("Missing columns for the input data:\n  ", cols)
  }
  invisible(NULL)
}

#' @export
checkKeys <- function(x, keys) {
  if (!is.list(x)) {
    stop("'x' must be a list.")
  }
  x <- x[keys]
  keys <- keys[lengths(x) == 0]
  if (length(keys) == 0) {
    return(invisible(NULL))
  }
  keys <- paste(keys, collapse = ", ")
  stop("Missing keys for the input list:\n  ", keys)
}

#' Check if specific packages are installed
#'
#' @param packages A character vector including package names to check.
#' @param required Logical, whether or not the `packages` are necessary.
#' @param ... `r .dot_param` See [requireNamespace()].
#'
#' @returns
#' When any of `packages` is not installed, a warning or an error will be raised
#' according to `required`.
#'
#' @export
checkPackages <- function(packages = character(), required = TRUE, ...) {
  chk <- vapply(
    X = packages,
    FUN = requireNamespace,
    FUN.VALUE = logical(1L),
    quietly = TRUE,
    ...
  )
  not.found <- packages[!chk]
  if (length(not.found) == 0) {
    return(invisible(NULL))
  }
  not.found <- paste(not.found, collapse = ", ")
  e <- paste0("Required packages not found:\n ", not.found)
  if (required) {
    stop(e)
  }
  fastWarning(e)
}

#' @export
validCharacters <- function(x) {
  nzchar(x, keepNA = FALSE) & !is.na(x) & is.atomic(x)
}
