#' Calculate the row means of logged values
#'
#' Calculate the row means of logged values in non-log space.
#'
#' @param object An object
#' @param ... `r .dot_param`
#'
#' @returns
#' A numeric vector containing the row means of `object`. If `object` is a 1-D
#' array, return a numeric scalar.
#'
#' @name rowExpMean
#' @export rowExpMean
rowExpMean <- function(object, ...) {
  UseMethod("rowExpMean", object)
}

#' @param log1p Logical, whether to transform the result into log-space. Default
#' is `FALSE`.
#'
#' @rdname rowExpMean
#' @export
#' @method rowExpMean numeric
rowExpMean.numeric <- function(object, log1p = FALSE, ...) {
  out <- mean(expm1(object), ...)
  if (isTRUE(log1p)) {
    return(log1p(out))
  }
  out
}

#' @importFrom MatrixGenerics rowMeans
#' @importClassesFrom Matrix CsparseMatrix
#'
#' @rdname rowExpMean
#' @export
#' @method rowExpMean CsparseMatrix
rowExpMean.CsparseMatrix <- function(object, log1p = FALSE, ...) {
  out <- rowMeans(expm1(object), ...)
  if (isTRUE(log1p)) {
    return(log1p(out))
  }
  out
}

#' @rdname rowExpMean
#' @export
#' @method rowExpMean matrix
rowExpMean.matrix <- function(object, log1p = FALSE, ...) {
  rowExpMean.CsparseMatrix(object = object, log1p = log1p, ...)
}
