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

#' @param log Logical, whether to transform the result into log-space. Default
#' is `FALSE`.
#'
#' @rdname rowExpMean
#' @export
#' @method rowExpMean numeric
rowExpMean.numeric <- function(object, log = FALSE, ...) {
  out <- mean(expm1(object))
  if (log) {
    return(log1p(out))
  }
  out
}

#' @importFrom Matrix rowMeans
#' @importClassesFrom Matrix CsparseMatrix
#'
#' @rdname rowExpMean
#' @export
#' @method rowExpMean CsparseMatrix
rowExpMean.CsparseMatrix <- function(object, log = FALSE, ...) {
  out <- Matrix::rowMeans(expm1(object), ...)
  if (log) {
    return(log1p(out))
  }
  out
}

#' @rdname rowExpMean
#' @export
#' @method rowExpMean matrix
rowExpMean.matrix <- function(object, log = FALSE, ...) {
  rowExpMean.CsparseMatrix(object = object, log = log, ...)
}

#' Calculate the means and percents of values
#'
#' Calculate the means (in log-spcae) and positive (typically non-zero) percents
#' for each row in a matrix, where columns are grouping by a factor vector.
#'
#' @param object A matrix-like object
#' @param ... `r .dot_param`
#'
#' @name rowMeanPct
#' @export rowMeanPct
rowMeanPct <- function(object, ...) {
  UseMethod("rowMeanPct", object)
}

#' @param cell.groups One of the following:
#' \itemize{
#' \item A factor vector containing group labels for columns of `object`. The
#' length must be equal to `ncol(object)`.
#' \item A list containing the column indices (integer or character if `object`
#' has column names) for each group.
#' }
#'
#' @param mean.fxn A function used to calculate the row means of `object`.
#' Default is \code{\link{rowExpMean}} in non-logged space.
#' @param min.exp A numeric scalar indicating the minimum value to calculate the
#' positive percents. Default is 0.
#'
#' @returns
#' A list containing three matrices where each row stands for a feature and each
#' column stands for a group:
#' \itemize{
#' \item `n.cells`: The numbers of positive observations for each feature within
#' corresponding groups.
#' \item `avg.exp`: The mean values calculated from `mean.fxn`.
#' \item `avg.pct`: The percents of positive observations within corresponding
#' groups.
#' }
#'
#' @importFrom Matrix rowSums
#' @importClassesFrom Matrix CsparseMatrix
#'
#' @rdname rowMeanPct
#' @export
#' @method rowMeanPct CsparseMatrix
rowMeanPct.CsparseMatrix <- function(
    object,
    cell.groups,
    mean.fxn = rowExpMean,
    min.exp = 0,
    ...
) {
  if (is.character(cell.groups)) {
    cell.groups <- as.factor(cell.groups)
  }
  if (is.factor(cell.groups)) {
    if (!length(cell.groups) == ncol(object)) {
      stop(
        "`cell.groups` is a factor whose length (", length(cell.groups), ") is",
        " not equal to column number (", ncol(object), ") of `object`."
      )
    }
    cell.groups <- split(seq_len(ncol(object)), f = cell.groups)
  }
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
  list(n.cells = n.cells, avg.exp = avg.exp, avg.pct = avg.pct)
}

#' @rdname rowMeanPct
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
