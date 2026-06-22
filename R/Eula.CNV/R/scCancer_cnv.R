
.filter_gene_cnv <- function(expr, cutoff = 0.1, minCell = 3) {
  gene.mean <- Matrix::rowMeans(expr)
  gene.sum <- Matrix::rowSums(expr > 0)
  genes.sel <- rownames(expr)[gene.mean >= cutoff & gene.sum >= minCell]

  expr <- expr[genes.sel, , drop = FALSE]
  expr
}

#' @importFrom Eula.matrix relativeCounts
.normalize_cnv <- function(expr) {
    cs <- Matrix::colSums(expr)
    normalize.factor <- 10^round(log10(mean(cs)))
    relativeCounts(expr, scale.factor = normalize.factor)
}
