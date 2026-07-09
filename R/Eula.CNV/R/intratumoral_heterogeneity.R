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

#' @export
runStemness <- function(X, stem.sig = NULL, species = c("human", "mouse")) {
  if (is.null(stem.sig)) {
    species <- match.arg(species)
    stem.sig.file <- system.file(
      "txt", paste0("pcbc-stemsig-", species, ".tsv"),
      package = packageName()
    )
    stem.sig <- read.delim(stem.sig.file, header = FALSE, row.names = 1)
  }

  common.genes <- intersect(rownames(stem.sig), rownames(X))
  X <- X[common.genes, , drop = FALSE]
  stem.sig <- stem.sig[common.genes, , drop = FALSE]

  s <- apply(X, 2, function(z) {
    cor(z, stem.sig, method = "sp", use = "complete.obs")
  })
  names(s) <- colnames(X)

  s <- s - min(s)
  s <- s / max(s)
  s
}

#' @export
getGenesets <- function(
  geneset = c("hallmark-pathways", "CancerSEA"),
  species = c("human", "mouse")
) {
  species <- match.arg(species)
  geneset <- match.arg(geneset)
  file.name <- paste0(geneset, "-", species, ".txt")
  file.name <- system.file("txt", file.name, package = packageName())
  genesets <- read.table(
    file.name,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  genesets <- split(genesets[[2]], f = genesets[[1]])
  genesets <- lapply(genesets, function(x) unique(unlist(strsplit(x, ", "))))
  genesets
}
