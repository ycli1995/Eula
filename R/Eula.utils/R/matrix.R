
#' @export
asMatrix <- function(object, max.row = 1000) {
  breaks <- 0:ceiling(nrow(object) / max.row) * max.row
  k <- cut(1:nrow(object), breaks = breaks)
  mat <- matrix(0, nrow(object), ncol(object), dimnames = dimnames(object))
  for (i in levels(k)) {
    message("as.matrix : ", paste0(range(which(k == i)), collapse = "-") )
    mat[which(k == i)] <- as.matrix(object[which(k == i), ])
  }
  mat
}

#' @export
min_max <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#' @export
min_max_cut <- function(x, limits = NULL) {
  if (length(limits) != 2) {
    fastWarning("'limits' should be a 2-element numeric vector.")
    return(x)
  }
  x[x > limits[2]] <- limits[2]
  x[x < limits[1]] <- limits[1]
  x
}

#' @export rowMeanPct
rowMeanPct <- function(object, ...) {
  UseMethod("rowMeanPct", object)
}

#' @export
#' @method rowMeanPct CsparseMatrix
rowMeanPct.CsparseMatrix <- function(
    object,
    group.by,
    mean.fxn = rowExpMean,
    min.exp = 0,
    ...
) {
  if (length(group.by) != ncol(object)) {
    stop(
      "Length of 'group.by' (", length(group.by), ") is not equal to ",
      "'ncol(object)' (", ncol(object), ")"
    )
  }
  if (!is.function(mean.fxn)) {
    stop("'mean.fxn' must be a function.")
  }
  group.by <- as.factor(group.by)
  group.by <- droplevels(group.by)
  all.groups <- levels(group.by)
  all.features <- rownames(object)

  n.cells <- avg.exp <- avg.pct <- matrix(
    data = 0,
    nrow = nrow(object),
    ncol = length(all.groups),
    dimnames = list(all.features, all.groups)
  )
  for (i in all.groups) {
    cells <- which(group.by == i)
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
