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

#' @export
checkPackages <- function(packages, required = TRUE, ...) {
  for (i in packages) {
    if (!requireNamespace(package = i, ..., quietly = TRUE)) {
      e <- sprintf("Please install package '%s' manually.", i)
      if (required) {
        stop(e)
      }
      fastWarning(e)
    }
  }
  invisible(NULL)
}

#' @export
validCharacters <- function(x) {
  nzchar(x, keepNA = FALSE) & !is.na(x) & is.atomic(x)
}

