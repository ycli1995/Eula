
random_matrix <- function(nrow, ncol, density = 0.5, max.val = 10) {
  m <- matrix(
    rbinom(nrow * ncol, 1, density) * 
      sample.int(max.val, nrow * ncol, replace = TRUE),
    nrow = nrow
  )
  m <- as(Matrix::Matrix(m, sparse = TRUE), "CsparseMatrix")
  colnames(m) <- paste0("col", seq_len(ncol(m)))
  rownames(m) <- paste0("row", seq_len(nrow(m)))
  m
}
