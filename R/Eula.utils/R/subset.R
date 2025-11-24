#' @include verbose.R
NULL

#' @export
filterData <- function(x, ...) {
  UseMethod("filterData", x)
}

#' @export
#' @method filterData character
filterData.character <- function(x, include = NULL, exclude = NULL, ...) {
  include <- as.character(include)
  exclude <- as.character(exclude)
  if (length(include) > 0) {
    include <- x %in% include
  } else {
    include <- rep.int(TRUE, length(x))
  }
  if (length(exclude) > 0) {
    exclude <- !x %in% exclude
  } else {
    exclude <- rep.int(FALSE, length(x))
  }
  (include & !exclude)
}

#' @export
#' @method filterData factor
filterData.factor <- function(x, include = NULL, exclude = NULL, ...) {
  filterData.character(x = x, include = include, exclude = exclude)
}

#' @export
#' @method filterData logical
filterData.logical <- function(x, ...) {
  x
}

#' @export
#' @method filterData numeric
filterData.numeric <- function(x, include = NULL, exclude = NULL, ...) {
  include <- as.list(include)
  exclude <- as.list(exclude)

  include <- include[lengths(include) == 2]
  include <- lapply(include, as.numeric)
  exclude <- exclude[lengths(exclude) == 2]
  exclude <- lapply(exclude, as.numeric)

  if (length(include) > 0) {
    include <- Reduce("|", lapply(include, function(r) x >= r[1] & x <= r[2]))
  } else {
    include <- rep.int(TRUE, length(x))
  }
  if (length(exclude) > 0) {
    exclude <- Reduce("|", lapply(exclude, function(r) x >= r[1] & x <= r[2]))
  } else {
    exclude <- rep.int(FALSE, length(x))
  }
  (include & !exclude)
}

#' @export
#' @method filterData data.frame
filterData.data.frame <- function(
    x,
    column,
    include = NULL,
    exclude = NULL,
    verbose = TRUE,
    ...
) {
  if (!column %in% colnames(x)) {
    fastWarning("'", column, "' not found in data. No filtering.")
    return(rep.int(TRUE, nrow(x)))
  }
  verboseMsg("Filtering data by column '", column, "'.")
  filterData(x = x[[column]], include = include, exclude = exclude)
}
