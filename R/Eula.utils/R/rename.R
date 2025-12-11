#' @include verbose.R
#'
#' @importFrom magrittr %>%
NULL

#' Change factor levels
#'
#' @param x An object to modify
#' @param ... `r .dot_param`
#'
#' @name recodeFactor
#' @export recodeFactor
recodeFactor <- function(x, ...) {
  UseMethod("recodeFactor", x)
}

#' @param map A list mapping new levels to old levels. Must formatted as
#' following:
#' ```
#' list(
#'   new.level1 = c(old.level1, old.level3, ...),
#'   new.level2 = old.level2,
#'   new.level3 = old.level4,
#'   ...
#' )
#' ```
#' @param keep.order Logical, whether or not to keep the orignial order of old
#' levels. Default is `FALSE`, which means the new order should depend on names
#' of `map`.
#'
#' @returns
#' A factor or character vector that is recoded according to `map`.
#'
#' @importFrom forcats fct_recode fct_relevel
#'
#' @rdname recodeFactor
#' @export
#' @method recodeFactor default
recodeFactor.default <- function(x, map, keep.order = FALSE, ...) {
  map <- unlistMap(map)
  old.class <- class(x)
  x <- as.factor(x)
  map <- map[map %in% levels(x)]
  x <- fct_recode(x, !!!map)
  if (keep.order) {
    return(as(x, old.class))
  }
  x <- fct_relevel(x, unique(names(map)))
  as(x, old.class)
}

#' @param column Character or integer scalar, indicating which column to be
#' recoded.
#' @param verbose `r .vb_param`
#'
#' @rdname recodeFactor
#' @export
#' @method recodeFactor data.frame
recodeFactor.data.frame <- function(
    x,
    map,
    column,
    keep.order = FALSE,
    verbose = TRUE,
    ...
) {
  column <- column[1]
  if (!column %in% colnames(x)) {
    return(x)
  }
  verboseMsg("Replace entries in column '", column, "':")
  column.old <- paste0(column, "_old")
  x[, column.old] <- x[, column]
  x[, column] <- recodeFactor(x[, column], map = map, keep.order = keep.order)
  if (verbose) {
    captureMsg(table(x[, column], x[, column.old]))
  }
  x
}

#' @export
unlistMap <- function(map, reverse = FALSE) {
  map <- stack(map)
  if (reverse) {
    map <- setNames(as.character(map$ind), as.character(map$values))
    return(map[unique(names(map))])
  }
  map <- setNames(as.character(map$values), as.character(map$ind))
  map[!duplicated(map)]
}

#' @export
pasteFactors <- function(x, y, collapse = "_", rev = FALSE) {
  if (!is.factor(x)) {
    x <- as.factor(x)
  }
  if (!is.factor(y)) {
    y <- as.factor(y)
  }
  lv <- sapply(levels(x), paste0, collapse, levels(y), simplify = TRUE)
  if (rev) {
    lv <- t(lv)
  }
  lv <- as.vector(lv)
  str <- factor(paste0(x, collapse, y), lv)
  str
}
