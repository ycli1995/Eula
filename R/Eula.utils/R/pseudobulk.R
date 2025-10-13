
#' @export pseudobulkMatrix
pseudobulkMatrix <- function(object, ...) {
  UseMethod("pseudobulkMatrix", object)
}

#' @export
#' @method pseudobulkMatrix matrix
pseudobulkMatrix.matrix <- function(
    object,
    cell.groups,
    method = c("aggregate", "average"),
    reverse = FALSE,
    bulk = FALSE,
    ...
) {
  pseudobulkMatrix.CsparseMatrix(
    object = object,
    cell.groups = cell.groups,
    method = method,
    reverse = reverse,
    bulk = bulk,
    sparse = FALSE,
    ...
  )
}

#' @importClassesFrom Matrix CsparseMatrix
#' @export
#' @method pseudobulkMatrix CsparseMatrix
pseudobulkMatrix.CsparseMatrix <- function(
    object,
    cell.groups,
    method = c("aggregate", "average"),
    reverse = FALSE,
    bulk = FALSE,
    sparse = FALSE,
    ...
) {
  method <- match.arg(method)
  if (length(cell.groups) != ncol(object)) {
    stop("'cell.groups' must have the same length as 'ncol(object)'")
  }
  cell.groups <- as.factor(cell.groups)
  cat.mat <- categoryMatrix(cell.groups, method = method, reverse = reverse)
  mat <- object %*% cat.mat
  rownames(mat) <- rownames(object)
  colnames(mat) <- levels(cell.groups)
  if (bulk) {
    func <- switch(
      EXPR = method,
      "aggregate" = Matrix::rowSums,
      "average" = Matrix::rowMeans
    )
    bulk.mat <- func(object)
    mat <- cbind("bulk" = bulk.mat, mat)
  }
  if (sparse) {
    if (inherits(mat, "CsparseMatrix")) {
      return(mat)
    }
    return(as(mat, "CsparseMatrix"))
  }
  as.matrix(mat)
}

#' @importFrom Matrix sparse.model.matrix
#' @export
categoryMatrix <- function(
    labels,
    method = c("aggregate", "average"),
    reverse = FALSE
) {
  method <- match.arg(method)
  if (length(unique(labels)) == 1) {
    message("All grouping variables have 1 value only.")
    if (reverse) {
      stop("Cannot use 'reverse = TRUE'.")
    }
    cat.mat <- matrix(1, length(labels), 1)
    if (method == "sum") {
      return(cat.mat)
    }
    return(cat.mat / sum(cat.mat))
  }
  labels <- as.factor(labels)
  cat.mat <- sparse.model.matrix(~0 + labels)
  colnames(cat.mat) <- levels(labels)
  if (reverse) {
    cat.mat <- as(1 - cat.mat, "CsparseMatrix")
  }
  colsums <- Matrix::colSums(cat.mat)
  cat.mat <- cat.mat[, colsums > 0, drop = FALSE]
  colsums <- colsums[colsums > 0]
  if (method == "aggregate") {
    return(cat.mat)
  }
  cat.mat@x <- cat.mat@x / rep.int(colsums, diff(cat.mat@p))
  cat.mat
}
