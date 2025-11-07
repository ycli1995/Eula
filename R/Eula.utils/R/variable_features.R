
#' @export rowExpMean
rowExpMean <- function(object, ...) {
  UseMethod("rowExpMean", object)
}

#' @export
#' @method rowExpMean CsparseMatrix
rowExpMean.CsparseMatrix <- function(object, ...) {
  log1p(Matrix::rowMeans(expm1(object), ...))
}

#' @export
#' @method rowExpMean matrix
rowExpMean.matrix <- function(object, ...) {
  rowExpMean.CsparseMatrix(object = object, ...)
}
