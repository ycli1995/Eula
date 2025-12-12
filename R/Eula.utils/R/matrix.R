#' @include verbose.R
NULL

#' Convert to dense matrix
#'
#' Convert a matrix-like object into dense matrix in row chunks.
#'
#' @param object An object
#' @param ... `r .dot_param`
#'
#' @returns
#' A dense `matrix` converted from `object`.
#'
#' @name asMatrix
#' @export asMatrix
asMatrix <- function(object, ...) {
  UseMethod("asMatrix", object)
}

#' @param chunk.size An integer indicating number of rows to convert in each
#' chunk.
#' @param verbose `r .vb_param`
#'
#' @rdname asMatrix
#' @export
#' @method asMatrix default
asMatrix.default <- function(object, chunk.size = 1000, verbose = TRUE) {
  breaks <- 0:ceiling(nrow(object) / chunk.size) * chunk.size
  k <- cut(seq_len(nrow(object)), breaks = breaks)
  mat <- matrix(0, nrow(object), ncol(object), dimnames = dimnames(object))
  for (i in levels(k)) {
    verboseMsg("as.matrix : ", paste0(range(which(k == i)), collapse = "-"))
    mat[which(k == i), ] <- as.matrix(object[which(k == i), ])
  }
  mat
}

#' @rdname asMatrix
#' @export
#' @method asMatrix matrix
asMatrix.matrix <- function(object, ...) {
  object
}

#' @rdname asMatrix
#' @export
#' @method asMatrix data.frame
asMatrix.data.frame <- function(object, ...) {
  as.matrix(object)
}

#' @export
normMinMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#' @export
cutMinMax <- function(x, limits = NULL) {
  if (length(limits) != 2) {
    fastWarning("'limits' should be a 2-element numeric vector.")
    return(x)
  }
  x[x > limits[2]] <- limits[2]
  x[x < limits[1]] <- limits[1]
  x
}

#' @export rowExpMean
rowExpMean <- function(object, ...) {
  UseMethod("rowExpMean", object)
}

#' @importFrom Matrix rowMeans
#' @importClassesFrom Matrix CsparseMatrix
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

#' @export rowMeanPct
rowMeanPct <- function(object, ...) {
  UseMethod("rowMeanPct", object)
}

#' @importFrom Matrix rowSums
#' @importClassesFrom Matrix CsparseMatrix
#' @export
#' @method rowMeanPct CsparseMatrix
rowMeanPct.CsparseMatrix <- function(
    object,
    cell.groups,
    mean.fxn = rowExpMean,
    min.exp = 0,
    ...
) {
  if (!is.list(cell.groups)) {
    stop("'cell.groups' must be a list containing cells for each group.")
  }
  if (!is.function(mean.fxn)) {
    stop("'mean.fxn' must be a function.")
  }
  all.groups <- names(cell.groups)
  all.features <- rownames(object)
  n.cells <- avg.exp <- avg.pct <- matrix(
    data = 0,
    nrow = nrow(object),
    ncol = length(cell.groups),
    dimnames = list(all.features, all.groups)
  )
  for (i in seq_along(cell.groups)) {
    cells <- cell.groups[[i]]
    n.cells[, i] <- Matrix::rowSums(object[, cells, drop = FALSE] > min.exp)
    avg.pct[, i] <- round(n.cells[, i] / length(cells), digits = 3)
    avg.exp[, i] <- mean.fxn(object[, cells, drop = FALSE])
  }
  list(n.cells = n.cells, avg.pct = avg.pct, avg.exp = avg.exp)
}

#' @export
#' @method rowMeanPct matrix
rowMeanPct.matrix <- function(
    object,
    group.by,
    mean.fxn = rowExpMean,
    min.exp = 0,
    ...
) {
  rowMeanPct.CsparseMatrix(
    object = object,
    group.by = group.by,
    mean.fxn = mean.fxn,
    min.exp = min.exp,
    ...
  )
}
