#' @include verbose.R
#'
#' @importFrom magrittr %>%
NULL

#' @export recodeFactor
recodeFactor <- function(x, ...) {
  UseMethod("recodeFactor", x)
}

#' @importFrom dplyr recode recode_factor
#'
#' @export
#' @method recodeFactor default
recodeFactor.default <- function(x, map, keep.orders = FALSE, ...) {
  if (keep.orders) {
    recode.func <- dplyr::recode
  } else {
    recode.func <- dplyr::recode_factor
  }
  map <- unlistMap(map = map)
  x <- x %>%
    as.factor() %>%
    recode.func(!!!map) %>%
    as(Class = class(x))
  x
}

#' @export
#' @method recodeFactor data.frame
recodeFactor.data.frame <- function(
    x,
    map,
    column,
    keep.orders = FALSE,
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
  x[, column] <- recodeFactor(x[, column], map = map, keep.orders = keep.orders)
  if (verbose) {
    captureMsg(table(x[, column], x[, column.old]))
  }
  x
}

#' @export
unlistMap <- function(map) {
  map <- stack(map)
  map <- setNames(as.character(map$ind), nm = map$values)
  map[unique(names(map))]
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
