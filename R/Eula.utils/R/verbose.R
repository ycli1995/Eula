
#' Simple verbose message wrapper
#'
#' @param ... Pass to \code{\link{message}}
#' @param verbose Whether or not to show the message. If is \code{NULL}, will
#' search \code{verbose} variable in \code{\link{parent.frame}}.
#'
#' @return Print the progress to console when \code{verbose} is \code{TRUE}.
#'
#' @export
verboseMsg <- function(..., verbose = NULL) {
  verbose <- verbose %||% parent.frame()$verbose %||% TRUE
  if (isTRUE(x = verbose)) {
    message(...)
  }
  return(invisible())
}

#' @export
fastWarning <- function(..., call. = FALSE, immediate. = TRUE) {
  warning(..., call. = call., immediate. = immediate.)
}

