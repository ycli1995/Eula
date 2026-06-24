
#' @export
runCNV <- function(
  expr,
  chr.data,
  ref.cells,
  cut.off = 0.1,
  min.cell = 3,
  window.size = 101,
  center.method = "median",
  sd.amplifier = 1.0
) {
  expr <- .filter_gene_cnv(expr, cut.off = cut.off, min.cell = min.cell)
  expr <- .normalize_cnv(expr)
  expr <- .anscombe_transform_cnv(expr)
  expr <- .log2p1_cnv(expr)
  expr <- .bound_expr_cnv(expr, thres = .get_avg_bounds(expr))
  expr <- .smooth_by_chr_cnv(expr, chr.data, window.size = window.size)
  expr <- .center_expr_cnv(expr, center.method = center.method)
  expr <- .subtract_ref_expr_cnv(expr, ref.cells = ref.cells, inv.log = TRUE)
  expr <- .denoise_expr_cnv(expr, sd.amplifier = sd.amplifier)
  expr <- .remove_outliers_cnv(expr)
  expr
}


.filter_gene_cnv <- function(expr, cut.off = 0.1, min.cell = 3) {
  gene.mean <- Matrix::rowMeans(expr)
  gene.sum <- Matrix::rowSums(expr > 0)
  genes.sel <- rownames(expr)[gene.mean >= cut.off & gene.sum >= min.cell]

  expr <- expr[genes.sel, , drop = FALSE]
  expr
}

#' @importFrom Eula.matrix relativeCounts
.normalize_cnv <- function(expr) {
  cs <- Matrix::colSums(expr)
  normalize.factor <- 10^round(log10(mean(cs)))
  relativeCounts(expr, scale.factor = normalize.factor)
}

.anscombe_transform_cnv <- function(expr) {
  2 * sqrt(as.matrix(expr) + 0.375)
}

.log2p1_cnv <- function(expr) {
  log2(expr + 1)
}

#' @importFrom MatrixGenerics colMaxs colMins
.get_avg_bounds <- function(expr) {
  lower.bound <- mean(MatrixGenerics::colMins(expr))
  upper.bound <- mean(MatrixGenerics::colMaxs(expr))
  mean(abs(c(lower.bound, upper.bound)))
}

.bound_expr_cnv <- function(expr, thres) {
  expr[expr > thres] <- thres
  expr[expr < -thres] <- -thres
  expr
}

.smooth_one_cnv <- function(ori.data, window.size = 101) {
  half.window <- (window.size - 1) / 2

  pad.data <- c(
    rep.int(0, half.window),
    ori.data,
    rep.int(0, half.window)
  )
  bool.data <- c(
    rep.int(0, half.window),
    rep.int(1, length(ori.data)),
    rep.int(0, half.window)
  )

  kernel.vec <- c(1:half.window, half.window + 1, half.window:1)

  sum.data <- filter(pad.data, kernel.vec, sides = 2)
  num.data <- filter(bool.data, kernel.vec, sides = 2)
  sum.data <- sum.data[!is.na(sum.data)]
  num.data <- num.data[!is.na(num.data)]

  smo.data <- sum.data / num.data
  smo.data
}

.smooth_by_chr_cnv <- function(expr, chr.data, window.size = 101) {
  chr.data <- chr.data[rownames(expr), , drop = FALSE]
  chrList <- chr.data$CHR
  chrs <- as.character(unique(chrList))

  if (window.size < 2) {
    warning("Window length < 2, returning original data.", immediate. = TRUE)
    return(expr)
  }

  if (window.size %% 2 == 0) {
    window.size <- window.size + 1
    warning("'window.size' is even, adding 1.", immediate. = TRUE)
  }

  for (chr in chrs) {
    cur.genes.ix <- which(chrList == chr)
    if (length(cur.genes.ix) < 2) {
      next
    }

    cur.data <- expr[cur.genes.ix, , drop = FALSE]
    for (i in seq_len(ncol(cur.data))) {
      cur.data[, i] <- .smooth_one_cnv(cur.data[, i], window.size = window.size)
    }
    expr.data[cur.genes.ix, ] <- cur.data
  }
  expr.data
}

#' @importFrom MatrixGenerics colMeans colMedians
.center_across_chr_cnv <- function(expr, method = "median") {
  if (method == "median") {
    centers <- MatrixGenerics::colMedians(expr, na.rm = TRUE)
  } else {
    centers <- MatrixGenerics::colMeans(expr, na.rm = TRUE)
  }
  for (i in seq_len(nrow(expr))) {
    expr[i, ] <- expr[i, ] - centers
  }
  expr
}

#' @importFrom Matrix rowMeans
.subtract_ref_expr_cnv <- function(expr, ref.cells, inv.log = TRUE) {
  if (inv.log) {
    ref.means <- log2(
      Matrix::rowMeans(2^expr[, ref.cells, drop = FALSE] - 1) + 1
    )
  } else {
    ref.means <- Matrix::rowMeans(expr[, ref.cells, drop = FALSE])
  }
  expr <- expr - ref.means
  if (inv.log) {
    expr <- 2^expr
  }
  expr
}

#' @importFrom MatrixGenerics colSds
.denoise_by_ref_mean_sd <- function(expr, ref.cells, sd.amplifier = 1.5) {
  mean.ref.vals <- mean(expr[, ref.cells])
  mean.ref.sd <- mean(MatrixGenerics::colSds(expr[, ref.cells, na.rm = TRUE]))
  mean.ref.sd <- mean.ref.sd * sd.amplifier

  up.bound <- mean.ref.vals + mean.ref.sd
  low.bound <- mean.ref.vals - mean.ref.sd

  expr[expr > low.bound & expr.data < up.bound] <- mean.ref.vals
  expr
}

.remove_outliers_cnv <- function(expr) {
  lower.bound <- mean(MatrixGenerics::colMins(expr))
  upper.bound <- mean(MatrixGenerics::colMaxs(expr))

  expr[expr < lower.bound] <- lower.bound
  expr[expr > upper.bound] <- upper.bound
  expr
}
