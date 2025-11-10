
#' @export rowExpMean
rowExpMean <- function(object, ...) {
  UseMethod("rowExpMean", object)
}

#' @export
#' @method rowExpMean CsparseMatrix
rowExpMean.CsparseMatrix <- function(object, log = FALSE, ...) {
  if (log) {
    return(log1p(Matrix::rowMeans(expm1(object), ...)))
  }
  Matrix::rowMeans(expm1(object), ...)
}

#' @export
#' @method rowExpMean matrix
rowExpMean.matrix <- function(object, log = FALSE, ...) {
  rowExpMean.CsparseMatrix(object = object, log = log, ...)
}
