#' @include verbose.R
NULL

#' Filter data according to set or ranges
#'
#' @param x An object
#' @param ... `r .dot_param`
#'
#' @export filterData
#' @name filterData
filterData <- function(x, ...) {
  UseMethod("filterData", x)
}

#' @returns
#' When `x` is a logical vector, return itself.
#'
#' @rdname filterData
#' @export
#' @method filterData logical
filterData.logical <- function(x, ...) {
  x
}

#' @param include One of the following:
#' \itemize{
#' \item A character vector indicating the element set to include.
#' \item A 2-element numeric vector indicating the range to include.
#' \item A list of several 2-element numeric vectors indicating multiple range
#' to match. Any vector that doesn't have 2 elements will be ignored.
#' }
#' @param invert Logical. If `TRUE` return the logical indices for elements that
#' are not within `include`.
#'
#' @returns
#' When `x` is a character or numeric vector, return a logical indices of the
#' elements in `x` that match `include`.
#'
#' @rdname filterData
#' @export
#' @method filterData character
filterData.character <- function(x, include = NULL, invert = FALSE, ...) {
  include <- as.character(include)
  if (length(include) > 0) {
    include <- x %in% include
  } else {
    include <- rep.int(TRUE, length(x))
  }
  if (invert) {
    return(!include)
  }
  include
}

#' @rdname filterData
#' @export
#' @method filterData factor
filterData.factor <- function(x, include = NULL, invert = FALSE, ...) {
  filterData.character(x = x, include = include, invert = invert)
}

#' @rdname filterData
#' @export
#' @method filterData numeric
filterData.numeric <- function(x, include = NULL, invert = FALSE, ...) {
  include <- .get_num_ranges(include)
  if (length(include) > 0) {
    include <- Reduce("|", lapply(include, function(r) x >= r[1] & x <= r[2]))
  } else {
    include <- rep.int(TRUE, length(x))
  }
  if (invert) {
    return(!include)
  }
  include
}

#' @param column A string indicating the column names used for filtering.
#' @param verbose `r .vb_param`
#'
#' @returns
#' When `x` is a `data.frame`, return a logical indices of the rows at which the
#' `column` matches `include`.
#'
#' @rdname filterData
#' @export
#' @method filterData data.frame
filterData.data.frame <- function(
    x,
    column,
    include = NULL,
    invert = FALSE,
    verbose = TRUE,
    ...
) {
  if (!column %in% colnames(x)) {
    fastWarning("'", column, "' not found in data. No filtering.")
    return(rep.int(TRUE, nrow(x)))
  }
  verboseMsg("Filtering data by column '", column, "'.")
  filterData(x = x[[column]], include = include, invert = invert, ...)
}

.get_num_ranges <- function(ranges) {
  if (length(ranges) == 0) {
    return(NULL)
  }
  if (is.numeric(ranges)) {
    ranges <- list(ranges)
  }
  ranges <- as.list(ranges)
  ranges <- lapply(ranges, as.numeric)
  ranges <- ranges[lengths(ranges) == 2]
  if (length(ranges) == 0) {
    stop("'ranges' must be a list of 2-element numeric vectors.")
  }
  ranges
}
