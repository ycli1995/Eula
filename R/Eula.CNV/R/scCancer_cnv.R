#' @importFrom collapse fquantile
#' @importFrom Eula.utils verboseMsg
#' @export
runMalignancy <- function(
  expr,
  celltypes,
  gene.chr = NULL,
  cut.off = 0.1,
  min.cell = 3,
  window.size = 101,
  center.method = "median",
  sd.amplifier = 1,
  adj.mat = NULL,
  ref.data = NULL,
  ref.adj.mat = NULL,
  species = "human",
  genome = "hg38",
  verbose = TRUE
) {
  if (is.null(ref.data)) {
    if (is.null(adj.mat)) {
      stop("`adj.mat` must be provided if `ref.data` is NULL")
    }
    verboseMsg("Get reference adjacency matrix")
    ref.adj.mat <- .get_adj_mat(species)
  }
  verboseMsg("Running `prepareCNV`")
  cnv.list <- prepareCNV(
    expr = expr,
    celltypes = celltypes,
    gene.chr = gene.chr,
    ref.data = ref.data,
    species = species,
    genome = genome
  )
  obs.anno <- cnv.list$cell.anno[colnames(expr), , drop = FALSE]
  ref.anno <- cnv.list$cell.anno[cnv.list$ref.cells, , drop = FALSE]
  verboseMsg("Running `runCNV`")
  cnv.expr <- runCNV(
    expr = cnv.list$expr,
    chr.data = cnv.list$gene.chr,
    ref.cells = cnv.list$ref.cells,
    cut.off = cut.off,
    min.cell = min.cell,
    window.size = window.size,
    center.method = center.method,
    sd.amplifier = sd.amplifier
  )

  verboseMsg("Running `getMalignScore` for reference cells")
  referScore.smooth <- getMalignScore(
    expr = cnv.expr,
    cells = cnv.list$ref.cells,
    method = "smooth",
    adjMat = ref.adj.mat
  )
  verboseMsg("Running `getMalignScore` for observation cells")
  obserScore.smooth <- getMalignScore(
    expr = cnv.expr,
    cells = setdiff(colnames(cnv.expr), cnv.list$ref.cells),
    method = "smooth",
    adjMat = adj.mat
  )
  up.refer <- fquantile(referScore.smooth, 0.995)
  low.refer <- fquantile(referScore.smooth, 0.005)
  referScore.smooth <- (referScore.smooth - low.refer) / (up.refer - low.refer)
  obserScore.smooth <- (obserScore.smooth - low.refer) / (up.refer - low.refer)

  verboseMsg("Calculate threshold values")
  all.thres <- .get_bimodal_thres(c(referScore.smooth, obserScore.smooth))
  malign.thres <- .get_bimodal_thres(obserScore.smooth)

  ## malignancy type
  verboseMsg("Assign malignancy type")
  if (!is.null(all.thres)) {
    malign.type <- rep("Malignant", length(obserScore.smooth))
    if (!is.null(malign.thres)) {
      malign.type[obserScore.smooth < malign.thres] <- "Non-malignant"
    }
  } else {
    malign.type <- rep("Non-malignant", length(obserScore.smooth))
    if (!is.null(malign.thres)) {
      malign.type[obserScore.smooth >= malign.thres] <- "Malignant"
    }
  }
  names(malign.type) <- names(obserScore.smooth)

  ## add score and type to cell.annotation
  verboseMsg("Add score and type to cell.annotation")
  obs.anno$Malign.score <- obserScore.smooth[rownames(obs.anno)]
  obs.anno$Malign.type <- malign.type[rownames(obs.anno)]
  ref.anno$Malign.score <- referScore.smooth[rownames(ref.anno)]
  ref.anno$Malign.type <- "Non-malignant"

  results <- list(
    cnv.expr = cnv.expr,
    obs.anno = obs.anno,
    ref.anno = ref.anno,
    gene.chr = gene.chr,
    malign.thres = malign.thres,
    all.thres = all.thres
  )
  results
}

#' @importFrom stats density
.get_bimodal_thres <- function(scores) {
  x.density <- density(scores)
  d.x.density <- diff(x.density$y)
  d.sign <- (d.x.density > 0) + 0

  ext.pos <- which(d.sign[-1] - d.sign[-length(d.sign)] != 0)
  ext.density <- x.density$y[ext.pos]
  y.max <- max(ext.density)
  if (length(ext.pos) >= 3) {
    del.ix <- c()
    for (ei in seq_along(ext.density)[-1]) {
      if (abs(ext.density[ei] - ext.density[ei - 1]) < y.max * 0.001) {
        del.ix <- c(del.ix, ei - 1, ei)
      }
    }
    sel.ix <- !(seq_along(ext.density) %in% unique(del.ix))
    ext.density <- ext.density[sel.ix]
    ext.pos <- ext.pos[sel.ix]
  }
  if (length(ext.pos) < 3) {
    return(NULL)
  }

  t.ext.density <- c(0, ext.density, 0)
  ext.height <- numeric(length(ext.pos))
  for (i in seq_along(ext.height)) {
    ext.height[i] <- min(
      abs(t.ext.density[i + 1] - t.ext.density[i]),
      abs(t.ext.density[i + 1] - t.ext.density[i + 2])
    )
  }
  ext <- data.frame(x = ext.pos, y = ext.density, height = ext.height)
  max.ix <- order(ext.density, decreasing = TRUE)
  threshold <- NULL
  if (ext.height[max.ix[2]] / ext.height[max.ix[1]] > 0.01) {
    cut.df <- ext[c(max.ix[2]:max.ix[1]), ]
    threshold <- x.density$x[cut.df[which.min(cut.df$y), ]$x]
  }
  threshold
}

#' @importFrom dplyr %>% arrange
#' @importFrom utils packageName
prepareCNV <- function(
  expr,
  celltypes,
  gene.chr = NULL,
  ref.data = NULL,
  species = "human",
  genome = "hg38"
) {
  ## gene.chr
  if (is.null(gene.chr)) {
    gene.chr <- .get_chr_data(expr, species, genome)
  }
  rownames(gene.chr) <- gene.chr$EnsemblID
  gene.chr$CHR <- as.factor(gene.chr$CHR)
  gene.chr$C_START <- as.numeric(gene.chr$C_START)
  gene.chr$C_STOP <- as.numeric(gene.chr$C_STOP)
  gene.chr <- gene.chr %>% dplyr::arrange(CHR, C_START, C_STOP)

  common.genes <- intersect(gene.chr$EnsemblID, rownames(expr))
  if (length(common.genes) == 0) {
    stop("No common genes between expr and gene.chr.\n")
  }
  expr <- expr[common.genes, , drop = FALSE]
  gene.chr <- gene.chr[common.genes, , drop = FALSE]

  ## cell.anno
  cell.anno <- data.frame(
    cellName = colnames(expr),
    cellAnno = as.factor(celltypes),
    stringsAsFactors = FALSE
  )

  ## reference.data
  if (is.null(ref.data)) {
    ref.data.file <- if (species == "human") {
      "cnvRef_Data-HM.RDS"
    } else if (species == "mouse") {
      "cnvRef_Data-boneMarrow-MS.RDS"
    } else {
      stop("species must be 'human' or 'mouse'.")
    }
    ref.data <- readRDS(system.file(
      "rds", ref.data.file,
      package = packageName()
    ))
  }
  ref.anno <- data.frame(
    cellName = colnames(ref.data),
    cellAnno = "Reference",
    stringsAsFactors = FALSE
  )
  ref.cells <- colnames(ref.data)

  ## combine data
  com.genes <- intersect(rownames(expr), rownames(ref.data))
  if (length(com.genes) == 0) {
    stop("No common genes between expr and ref.data.\n")
  }
  ref.data <- ref.data[com.genes, , drop = FALSE]
  expr <- expr[com.genes, , drop = FALSE]
  expr <- cbind(expr, ref.data)
  gene.chr <- gene.chr[com.genes, , drop = FALSE]

  cell.anno <- rbind(cell.anno, ref.anno)
  rownames(cell.anno) <- cell.anno$cellName

  list(
    expr = expr,
    gene.chr = gene.chr,
    cell.anno = cell.anno,
    ref.cells = ref.cells
  )
}

#' @importFrom utils packageName
.get_chr_data <- function(
  expr,
  species = c("human", "mouse"),
  genome = c("hg38", "hg19", "mm10")
) {
  species <- match.arg(species)
  genome <- match.arg(genome)
  if (species == "human") {
    if (genome == "hg38") {
      gene.chr.file <- system.file(
        "txt", "gene-chr-hg38.txt",
        package = packageName()
      )
    } else if (genome == "hg19") {
      gene.chr.file <- system.file(
        "txt", "gene-chr-hg19.txt",
        package = packageName()
      )
    } else {
      stop(genome, " is not allowed for 'human'.")
    }
  }
  if (species == "mouse") {
    if (genome == "mm10") {
      gene.chr.file <- system.file(
        "txt", "gene-chr-mm10.txt",
        package = "scCancer"
      )
    } else {
      stop(genome, " is not allowed for 'mouse'.\n")
    }
  }
  gene.chr <- read.table(
    gene.chr.file,
    sep = "\t",
    header = FALSE,
    col.names = c("EnsemblID", "CHR", "C_START", "C_STOP"),
    stringsAsFactors = FALSE
  )
  gene.chr
}

#' @importFrom utils packageName
.get_adj_mat <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  if (species == "human") {
    adj.mat.file <- system.file(
      "rds", "cnvRef_SNN-HM.RDS",
      package = packageName()
    )
  } else if (species == "mouse") {
    adj.mat.file <- system.file(
      "rds", "cnvRef_SNN-boneMarrow-MS.RDS",
      package = packageName()
    )
  } else {
    stop("species must be 'human' or 'mouse'.")
  }
  adj.mat <- readRDS(adj.mat.file)
  adj.mat
}

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
  expr <- .smooth_by_chr_cnv(
    expr,
    chr.data = chr.data,
    window.size = window.size
  )
  expr <- .center_across_chr_cnv(expr, method = center.method)
  expr <- .subtract_ref_expr_cnv(expr, ref.cells = ref.cells, inv.log = TRUE)
  expr <- .denoise_by_ref_mean_sd(
    expr,
    ref.cells = ref.cells,
    sd.amplifier = sd.amplifier
  )
  expr <- .remove_outliers_cnv(expr)
  expr
}

#' @importFrom collapse fquantile
#' @importFrom Matrix nnzero
#' @importFrom methods as
#' @export
getMalignScore <- function(
  expr,
  cells,
  method = "smooth",
  adjMat = NULL
) {
  cur.data <- expr[, cells, drop = FALSE]

  if (is.null(adjMat) & method == "smooth") {
    warning(
      "`adjMat` is not provided, use 'direct' method instead.",
      call. = FALSE, immediate. = TRUE
    )
    method <- "direct"
  }
  if (method == "direct") {
    malignScore <- colSums((cur.data - 1)^2) / nrow(cur.data)
    return(malignScore)
  }

  if (!is.null(adjMat)) {
    adjMat <- as(adjMat, "CsparseMatrix")
    adjMat <- adjMat[cells, cells, drop = FALSE]
  }
  thres <- collapse::fquantile(
    adjMat@x,
    1 - (nrow(adjMat) * 10 / Matrix::nnzero(adjMat))
  )

  indexes <- as.matrix(adjMat > thres)
  tt <- 0.5 / (rowSums(indexes) - 1)
  tt[is.infinite(tt)] <- 0

  indexes <- indexes * tt
  indexes <- indexes * (1 - diag(rep(1, nrow(indexes))))
  diagValue <- rep(0.5, nrow(indexes))
  diagValue[tt == 0] <- 1

  indexes <- t(indexes + diag(diagValue))

  new.cur.data <- tcrossprod(cur.data, indexes + diag(diagValue))
  colnames(new.cur.data) <- colnames(cur.data)
  malignScore <- colSums((new.cur.data - 1)^2) / nrow(new.cur.data)
  names(malignScore) <- colnames(new.cur.data)
  malignScore
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

#' @importFrom stats filter
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

  sum.data <- stats::filter(pad.data, kernel.vec, sides = 2)
  num.data <- stats::filter(bool.data, kernel.vec, sides = 2)
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
    expr[cur.genes.ix, ] <- cur.data
  }
  expr
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
  mean.ref.sd <- mean(
    MatrixGenerics::colSds(expr[, ref.cells, drop = FALSE], na.rm = TRUE)
  )
  mean.ref.sd <- mean.ref.sd * sd.amplifier

  up.bound <- mean.ref.vals + mean.ref.sd
  low.bound <- mean.ref.vals - mean.ref.sd

  expr[expr > low.bound & expr < up.bound] <- mean.ref.vals
  expr
}

.remove_outliers_cnv <- function(expr) {
  lower.bound <- mean(MatrixGenerics::colMins(expr))
  upper.bound <- mean(MatrixGenerics::colMaxs(expr))

  expr[expr < lower.bound] <- lower.bound
  expr[expr > upper.bound] <- upper.bound
  expr
}
