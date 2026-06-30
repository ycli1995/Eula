
ITHScore <- function(cell.embeddings, ...) {



  mu <- Matrix::rowMeans(cell.embeddings)
  sigma <- MatrixGenerics::rowSds(cell.embeddings)

}

.filter_cells_for_ith <- function(cell.embeddings) {
  filtered <- list()
  for (i in seq_len(3)) {
    
  }
  filtered
}
