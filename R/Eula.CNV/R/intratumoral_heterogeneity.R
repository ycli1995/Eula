#' @export 
ITHScore <- function(cell.embeddings, ...) {
  keep <- .filter_cells_for_ith(cell.embeddings)
  cell.embeddings <- cell.embeddings[keep, , drop = FALSE]

  mu <- Matrix::rowMeans(cell.embeddings)
  sigma <- MatrixGenerics::rowSds(cell.embeddings)

  score.per.cell <- Matrix::rowSums((cell.embeddings - mu)^2)
  score.per.cell <- sqrt(score.per.cell)
  mean(score.per.cell)
}

.filter_cells_for_ith <- function(cell.embeddings, ...) {
  filtered <- list()
  for (i in seq_len(3)) {
    filtered[[i]] <- .filter_outliers(cell.embeddings[, i], ...)
  }
  filtered <- Reduce("&", filtered)
  filtered
}

.filter_outliers <- function(x, sd = 3, na.rm = TRUE, ...) {
  mu <- mean(x, na.rm = na.rm)
  sigma <- sd(x, na.rm = na.rm)

  keep <- x >= (mu - sd * sigma) &
    x <= (mu + 3 * sigma)

  keep[is.na(x)] <- FALSE
  keep
}
