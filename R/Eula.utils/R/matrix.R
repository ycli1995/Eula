
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
