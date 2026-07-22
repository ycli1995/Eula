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
#' Default is [rowExpMean()] in non-logged space.
#' @param min.exp A numeric scalar indicating the minimum value to calculate the
#' positive percents. Default is 0.
#'
#' @details `...` are typically passed to `mean.fxn`.
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
#' @importFrom MatrixGenerics rowSums
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
  cell.groups <- as.factor(cell.groups)
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
    if (length(cell.groups[[i]]) > 0) {
      idx <- cell.groups[[i]]
      n.cells[, i] <- rowSums(object[, idx, drop = FALSE] > min.exp)
      avg.pct[, i] <- n.cells[, i] / length(idx)
      avg.exp[, i] <- mean.fxn(object[, idx, drop = FALSE], ...)
    }
  }
  list(n.cells = n.cells, avg.exp = avg.exp, avg.pct = avg.pct)
}

#' @rdname rowMeanPct
#' @export
#' @method rowMeanPct matrix
rowMeanPct.matrix <- function(
  object,
  cell.groups,
  mean.fxn = rowExpMean,
  min.exp = 0,
  ...
) {
  rowMeanPct.CsparseMatrix(
    object = object,
    cell.groups = cell.groups,
    mean.fxn = mean.fxn,
    min.exp = min.exp,
    ...
  )
}
