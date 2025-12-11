
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
  verbose <- verbose %||% (parent.frame()$verbose %||% TRUE)
  if (isTRUE(verbose)) {
    message(...)
  }
  invisible(NULL)
}

#' @export
fastWarning <- function(..., call. = FALSE, immediate. = TRUE) {
  warning(..., call. = call., immediate. = immediate.)
}

#' @export
pipeMsg <- function(
    ...,
    type = c("Running", "Warning", "Stopped"),
    pipe.name = NULL
) {
  type <- match.arg(type)
  pipe.name <- pipe.name %||% parent.frame()$pipe.name
  if (length(pipe.name) > 0) {
    pipe.name <- paste(pipe.name, collapse = "::")
    if (length(pipe.name) > 0) {
      pipe.name <- paste0("[", pipe.name, "]")
    }
  }
  time <- as.character(Sys.time())
  message("\n[", time, "][", type, "]", pipe.name, ": ", ...)
  if (type == "Stopped") {
    q(status = 1)
  }
  invisible(NULL)
}

#' @export
captureMsg <- function(...) {
  capture.output(..., file = stderr())
}

