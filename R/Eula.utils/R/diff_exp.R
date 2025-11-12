
#' @export
#' @concept data
DE.METHODS <- list(
  counts = c(),
  nocorrect = c("roc"),
  noprefilter = c("DESeq2"),
  latent = c("MAST")
)

#' @export differExp
differExp <- function(object, ...) {
  UseMethod("differExp", object)
}

#' @export
#' @method differExp CsparseMatrix
differExp.CsparseMatrix <- function(
    object,
    cells.1,
    cells.2,
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    max.cells.per.ident = Inf,
    min.cells.group = 3,
    logfc.threshold = 0.1,
    min.mean.exp = 0,
    min.pct = 0.01,
    min.diff.pct = -Inf,
    min.cells.feature = 3,
    only.pos = FALSE,
    mean.fxn = rowExpMean,
    min.exp = 0,
    pseudocount.use = 1,
    base = 2,
    seed = 42,
    latent.vars = NULL,
    densify = FALSE,
    ...
) {
  .validate_cell_groups(
    data.use = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    min.cells.group = min.cells.group
  )
  if (test.use %in% DE.METHODS$noprefilter) {
    min.diff.pct <- -Inf
    logfc.threshold <- 0
  }
  fc.results <- foldChange(
    object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    mean.fxn = mean.fxn,
    min.exp = min.exp,
    pseudocount.use = pseudocount.use,
    base = base
  )
  fc.results <- selectDE(
    fc.results = fc.results,
    logfc.threshold = logfc.threshold,
    min.mean.exp = min.mean.exp,
    min.pct = min.pct,
    min.diff.pct = min.diff.pct,
    min.cells.feature = min.cells.feature,
    only.pos = only.pos
  )

  # Actually perform the DE test
  testDE(
    object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    fc.results = fc.results,
    test.use = test.use,
    p.adjust.method = p.adjust.method,
    min.cells.group = min.cells.group,
    max.cells.per.ident = max.cells.per.ident,
    latent.vars = latent.vars,
    densify = densify,
    ...
  )
}

#' @importFrom stats p.adjust
#' @export
#' @method differExp matrix
differExp.matrix <- function(
    object,
    cells.1,
    cells.2,
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    max.cells.per.ident = Inf,
    min.cells.group = 3,
    logfc.threshold = 0.1,
    min.mean.exp = 0,
    min.pct = 0.01,
    min.diff.pct = -Inf,
    min.cells.feature = 3,
    only.pos = FALSE,
    mean.fxn = rowExpMean,
    min.exp = 0,
    pseudocount.use = 1,
    base = 2,
    seed = 42,
    latent.vars = NULL,
    densify = FALSE,
    ...
) {
  differExp.CsparseMatrix(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    test.use = test.use,
    p.adjust.method = p.adjust.method,
    max.cells.per.ident = max.cells.per.ident,
    min.cells.group = min.cells.group,
    logfc.threshold = logfc.threshold,
    min.mean.exp = min.mean.exp,
    min.pct = min.pct,
    min.diff.pct = min.diff.pct,
    min.cells.feature = min.cells.feature,
    only.pos = only.pos,
    mean.fxn = mean.fxn,
    min.exp = min.exp,
    pseudocount.use = pseudocount.use,
    base = base,
    seed = seed,
    latent.vars = latent.vars,
    densify = FALSE,
    ...
  )
}

#'@export
selectDE <- function(
    fc.results,
    logfc.threshold = 0.1,
    min.mean.exp = 0,
    min.pct = 0.01,
    min.diff.pct = -Inf,
    min.cells.feature = 3,
    only.pos = FALSE,
    verbose = TRUE
) {
  if (logfc.threshold < 0) {
    stop("'logfc.threshold' must be >= 0.")
  }
  if (only.pos) {
    selected <- fc.results[, 1] >= logfc.threshold
  } else {
    selected <- abs(fc.results[, 1]) >= logfc.threshold
  }
  ncells.min <- pmax(fc.results$ncells.1, fc.results$ncells.2)
  pct.max <- pmax(fc.results$pct.1, fc.results$pct.2)
  pct.diff <- pct.max - pmin(fc.results$pct.1, fc.results$pct.2)
  mean.max <- pmax(fc.results$mean.1, fc.results$mean.2)
  selected <- selected &
    (ncells.min >= min.cells.feature) &
    (pct.max >= min.pct) &
    (pct.diff >= min.diff.pct) &
    (mean.max >= min.mean.exp)

  features <- rownames(fc.results)[selected]
  if (length(features) == 0) {
    fastWarning(
      "No feature pass thresholds:\n ",
      "min.cells.feature = ", min.cells.feature, "\n ",
      "min.pct = ", min.pct, "\n ",
      "min.diff.pct = ", min.diff.pct, "\n ",
      "logfc.threshold = ", logfc.threshold, "\n ",
      "only.pos = ", only.pos
    )
  }
  if (verbose) {
    message("Select ", length(features), " features")
  }
  fc.results[features, , drop = FALSE]
}

#' @export testDE
testDE <- function(object, ...) {
  UseMethod("testDE", object)
}

#' @export
#' @method testDE CsparseMatrix
testDE.CsparseMatrix <- function(
    object,
    cells.1,
    cells.2,
    features = NULL,
    fc.results = NULL,
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    latent.vars = NULL,
    min.cells.group = 3,
    max.cells.per.ident = Inf,
    seed = 42,
    densify = FALSE,
    ...
) {
  if (!(test.use %in% DE.METHODS$latent) && !is.null(latent.vars)) {
    fastWarning(
      "'latent.vars' is only used for the following tests: ",
      paste(DE.METHODS$latent, collapse = ", "),
    )
  }
  .validate_cell_groups(
    data.use = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    min.cells.group = min.cells.group
  )
  total.features <- nrow(object)
  if (!is.null(fc.results)) {
    checkColumns(fc.results, c("pct.1", "pct.2"))
    features <- features %||% rownames(fc.results)
    fc.results <- fc.results[features, , drop = FALSE]
  }
  features <- features %||% rownames(object)
  if (length(features) == 0) {
    if (test.use == "roc") {
      empty.df <- as.data.frame(matrix(nrow = 0, ncol = 2))
      colnames(empty.df) <- c("AUC", "power")
    } else {
      empty.df <- as.data.frame(matrix(nrow = 0, ncol = 1))
      colnames(empty.df) <- "p_val"
    }
    if (!is.null(fc.results)) {
      return(cbind(fc.results[features, , drop = FALSE], empty.df))
    }
    return(empty.df)
  }

  # subsample cell groups if they are too large
  if (max.cells.per.ident < Inf) {
    if (length(cells.1) > max.cells.per.ident) {
      set.seed(seed)
      cells.1 <- sample(cells.1, size = max.cells.per.ident)
    }
    if (length(cells.2) > max.cells.per.ident) {
      set.seed(seed)
      cells.2 <- sample(cells.2, size = max.cells.per.ident)
    }
  }
  if (!is.null(latent.vars)) {
    latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
  }
  object <- object[features, c(cells.1, cells.2), drop = FALSE]
  if (densify) {
    object <- as.matrix(object)
  }
  de.results <- switch(
    EXPR = test.use,
    'wilcox' = differWilcox(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      ...
    ),
    'wilcox_limma' = differWilcox(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      limma = TRUE,
      ...
    ),
    'roc' = differROC(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2
    ),
    't' = differTTest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2
    ),
    'MAST' = differMAST(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      ...
    ),
    "DESeq2" = differDESeq2(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      ...
    ),
    stop("Unknown test: ", test.use)
  )
  if (!test.use %in% DE.METHODS$nocorrect) {
    de.results$p_val_adj <- p.adjust(
      p = de.results$p_val,
      method = p.adjust.method,
      n = total.features
    )
  }
  if (is.null(fc.results)) {
    if (test.use == "roc") {
      return(de.results[order(-de.results$AUC), , drop = FALSE])
    }
    return(de.results[order(de.results$p_val), , drop = FALSE])
  }
  de.results <- cbind(
    fc.results[rownames(de.results), , drop = FALSE],
    de.results
  )
  p.col <- "p_val"
  if (test.use == "roc") {
    p.col <- "AUC"
  }
  idx <- order(
    de.results[[p.col]],
    -de.results[, 1],
    -abs(de.results$pct.1 - de.results$pct.2)
  )
  de.results[idx, , drop = FALSE]
}

#' @export
#' @method testDE matrix
testDE.matrix <- function(
    object,
    cells.1,
    cells.2,
    features = NULL,
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    latent.vars = NULL,
    max.cells.per.ident = Inf,
    seed = 42,
    ...
) {
  testDE.CsparseMatrix(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    test.use = test.use,
    p.adjust.method = p.adjust.method,
    latent.vars = latent.vars,
    max.cells.per.ident = max.cells.per.ident,
    seed = seed,
    densify = FALSE,
    ...
  )
}

#' @export foldChange
foldChange <- function(object, ...) {
  UseMethod("foldChange", object)
}

#' @importFrom Matrix rowSums
#' @export
#' @method foldChange CsparseMatrix
foldChange.CsparseMatrix <- function(
    object,
    cells.1,
    cells.2,
    mean.fxn = rowExpMean,
    min.exp = 0,
    pseudocount.use = 1,
    base = 2,
    ...
) {
  group.by <- c(
    rep.int("cells.1", length(cells.1)),
    rep.int("cells.2", length(cells.2))
  )
  group.by <- factor(group.by, c("cells.1", "cells.2"))
  object <- object[, c(cells.1, cells.2), drop = FALSE]
  out <- rowMeanPct(
    object = object,
    group.by = group.by,
    mean.fxn = mean.fxn,
    min.exp = min.exp,
    ...
  )
  colnames(out$n.cells) <- c("ncells.1", "ncells.2")
  colnames(out$avg.pct) <- c("pct.1", "pct.2")
  colnames(out$avg.exp) <- c("mean.1", "mean.2")
  out <- out[c("avg.exp", "n.cells", "avg.pct")]
  out <- as.data.frame(Reduce(cbind, out))

  p1 <- pseudocount.use / length(cells.1)
  p2 <- pseudocount.use / length(cells.2)
  fc.results <- data.frame(
    fc = log(out$mean.1 + p1, base = base) - log(out$mean.2 + p2, base = base),
    row.names = rownames(object)
  )
  if (base == exp(1)) {
    base <- ""
  }
  fc.name <- paste0("avg_log", base, "FC")
  colnames(fc.results)[1] <- fc.name
  cbind(fc.results, out)
}

#' @export
#' @method foldChange matrix
foldChange.matrix <- function(
    object,
    cells.1,
    cells.2,
    mean.fxn,
    min.exp = 0,
    pseudocount.use = 1,
    base = 2,
    ...
) {
  foldChange.CsparseMatrix(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    mean.fxn = mean.fxn,
    min.exp = min.exp,
    pseudocount.use = pseudocount.use,
    base = base,
    ...
  )
}

#' @export differWilcox
differWilcox <- function(object, ...) {
  UseMethod("differWilcox", object)
}

#' @importFrom stats wilcox.test
#' @export
#' @method differWilcox CsparseMatrix
differWilcox.CsparseMatrix <- function(
    object,
    cells.1,
    cells.2,
    limma = FALSE,
    ...
) {
  presto.check <- requireNamespace("presto", quietly = TRUE)
  limma.check <- requireNamespace("limma", quietly = TRUE)
  my.sapply <- .get_sapply()
  if (!setequal(c(cells.1, cells.2), colnames(object))) {
    object <- object[, c(cells.1, cells.2), drop = FALSE]
  }
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(group.info[, "group"])
  group.info <- group.info[colnames(object), , drop = FALSE]

  if (presto.check[1] && (!limma)) {
    res <- presto::wilcoxauc(X = object, y = group.info[, "group"])
    res <- res[1:(nrow(res)/2), ]
    p_val <- res$pval
    return(data.frame(p_val, row.names = rownames(object)))
  }
  if (limma.check) {
    limma <- TRUE
  }
  if (!limma) {
    p_val <- my.sapply(seq_len(nrow(object)), function(i) {
      wilcox.test(object[i, ] ~ group.info[, "group"], ...)$p.value
    })
    return(data.frame(p_val, row.names = rownames(object)))
  }
  if (!limma.check) {
    stop(
      "To use the limma implementation of the Wilcoxon Rank Sum Test,
       please install the limma package:
       --------------------------------------------
       install.packages('BiocManager')
       BiocManager::install('limma')
       --------------------------------------------"
    )
  }
  tmp <- object[1, ]
  overflow <- is.na(suppressWarnings(length(tmp) * length(tmp)))
  if (overflow) {
    stop("Use 'limma' with ", ncol(object), " cells may cause overflow.")
  }
  j <- seq_along(cells.1)
  p_val <- my.sapply(seq_len(nrow(object)), function(i) {
    min(1, 2 * min(limma::rankSumTestWithCorrelation(j, object[i, ])))
  })
  data.frame(p_val, row.names = rownames(object))
}

#' @export
#' @method differWilcox matrix
differWilcox.matrix <- function(
    object,
    cells.1,
    cells.2,
    limma = FALSE,
    ...
) {
  differWilcox.CsparseMatrix(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    limma = limma,
    ...
  )
}

#' @export differTTest
differTTest <- function(object, ...) {
  UseMethod("differTTest", object)
}

#' @importFrom stats t.test
#' @export
#' @method differWilcox CsparseMatrix
differTTest.CsparseMatrix <- function(object, cells.1, cells.2, ...) {
  object.1 <- object[, cells.1, drop = FALSE]
  object.2 <- object[, cells.2, drop = FALSE]
  my.sapply <- .get_sapply()
  p_val <- my.sapply(
    seq_len(nrow(object)),
    function(i) t.test(object.1[i, ], object.2[i, ], ...)$p.value
  )
  data.frame(unlist(p_val), row.names = rownames(object))
}

#' @export
#' @method differWilcox matrix
differTTest.matrix <- function(object, cells.1, cells.2, ...) {
  differTTest.CsparseMatrix(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    ...
  )
}

#' @export differMAST
differMAST <- function(object, ...) {
  UseMethod("differMAST", object)
}

#' @method differMAST matrix
differMAST.matrix <- function(
    object,
    cells.1,
    cells.2,
    latent.vars = NULL,
    ...
) {
  if (!requireNamespace("MAST", quietly = TRUE)) {
    stop("Please install MAST - learn more at https://github.com/RGLab/MAST")
  }
  if (!setequal(c(cells.1, cells.2), colnames(object))) {
    object <- object[, c(cells.1, cells.2), drop = FALSE]
  }

  group.info <- data.frame(row.names = c(cells.1, cells.2))
  latent.vars <- latent.vars %||% group.info
  group.info[cells.1, "condition"] <- "Group1"
  group.info[cells.2, "condition"] <- "Group2"
  group.info[, "condition"] <- factor(group.info[, "condition"])

  latent.vars.names <- c("condition", colnames(latent.vars))
  latent.vars <- cbind(latent.vars, group.info)
  latent.vars$wellKey <- rownames(latent.vars)
  fdat <- data.frame(primerid = rownames(object))
  rownames(fdat) <- fdat[, 1]

  sca <- MAST::FromMatrix(
    exprsArray = object,
    check_sanity = FALSE,
    cData = latent.vars[colnames(object), , drop = FALSE],
    fData = fdat
  )
  SummarizedExperiment::colData(sca)$condition <- relevel(
    factor(SummarizedExperiment::colData(sca)$condition),
    ref = "Group1"
  )

  fmla <- as.formula(paste0(" ~ ", paste(latent.vars.names, collapse = "+")))

  zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
  summaryDt <- MAST::summary(zlmCond, doLRT = 'conditionGroup2')$datatable
  summaryDt <- summaryDt[summaryDt[, "component"] == "H", , drop = FALSE]
  data.frame(p_val = summaryDt[, 4], row.names = summaryDt[, 1])
}

#' @export
#' @method differMAST CsparseMatrix
differMAST.CsparseMatrix <- function(
    object,
    cells.1,
    cells.2,
    latent.vars = NULL,
    ...
) {
  if (!setequal(c(cells.1, cells.2), colnames(object))) {
    object <- object[, c(cells.1, cells.2), drop = FALSE]
  }
  differMAST(
    object = as.matrix(object),
    cells.1 = cells.1,
    cells.2 = cells.2,
    latent.vars = latent.vars,
    ...
  )
}

#' @export differROC
differROC <- function(object, ...) {
  UseMethod("differROC", object)
}

#' @importFrom ROCR performance prediction
#' @export
#' @method differROC CsparseMatrix
differROC.CsparseMatrix <- function(object, cells.1, cells.2, ...) {
  if (!setequal(c(cells.1, cells.2), colnames(object))) {
    object <- object[, c(cells.1, cells.2), drop = FALSE]
  }
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- 1
  group.info[cells.2, "group"] <- 0
  group.info <- group.info[colnames(object), , drop = FALSE]

  my.lapply <- .get_lapply()
  myAUC <- unlist(my.lapply(rownames(object), function(g) {
    prediction.use <- prediction(
      predictions = object[g, ],
      labels = group.info[, "group"],
      label.ordering = 0:1
    )
    perf.use <- performance(prediction.obj = prediction.use, measure = "auc")
    round(perf.use@y.values[[1]], digits = 3)
  }))
  myAUC[is.na(myAUC)] <- 0
  res <- data.frame(AUC = myAUC, row.names = rownames(object))
  res$power <- abs(res$AUC - 0.5) * 2
  res
}

#' @export
#' @method differROC matrix
differROC.matrix <- function(object, cells.1, cells.2, ...) {
  differROC.CsparseMatrix(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    ...
  )
}

#' @export differDESeq2
differDESeq2 <- function(object, ...) {
  UseMethod("differDESeq2", object)
}

#' @export
#' @method differDESeq2 matrix
differDESeq2.matrix <- function(object, cells.1, cells.2, ...) {
  if (!requireNamespace("DESeq2", quietly=TRUE)) {
    stop(
      "Please install DESeq2 - learn more at ",
      "https://bioconductor.org/packages/release/bioc/html/DESeq2.html"
    )
  }
  if (!setequal(c(cells.1, cells.2), colnames(object))) {
    object <- object[, c(cells.1, cells.2), drop = FALSE]
  }
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(group.info[, "group"])
  group.info$wellKey <- rownames(group.info)
  group.info <- group.info[colnames(object), , drop = FALSE]

  dds1 <- DESeq2::DESeqDataSetFromMatrix(
    countData = object,
    colData = group.info,
    design = ~ group
  )
  dds1 <- DESeq2::estimateSizeFactors(dds1)
  dds1 <- DESeq2::estimateDispersions(dds1, fitType = "local")
  dds1 <- DESeq2::nbinomWaldTest(dds1)
  res <- DESeq2::results(
    object = dds1,
    contrast = c("group", "Group1", "Group2"),
    alpha = 0.05,
    ...
  )
  data.frame(p_val = res$pvalue, row.names = rownames(res))
}

#' @export
#' @method differDESeq2 CsparseMatrix
differDESeq2.CsparseMatrix <- function(object, cells.1, cells.2, ...) {
  if (!setequal(c(cells.1, cells.2), colnames(object))) {
    object <- object[, c(cells.1, cells.2), drop = FALSE]
  }
  differDESeq2.matrix(
    object = as.matrix(object),
    cells.1 = cells.1,
    cells.2 = cells.2,
    ...
  )
}

.get_lapply <- function() {
  if (requireNamespace("future.apply", quietly = TRUE)) {
    if (future::nbrOfWorkers() > 1) {
      return(future.apply::future_lapply)
    }
  }
  return(lapply)
}

.get_sapply <- function() {
  if (requireNamespace("future.apply", quietly = TRUE)) {
    if (future::nbrOfWorkers() > 1) {
      return(future.apply::future_sapply)
    }
  }
  return(sapply)
}

.validate_cell_groups <- function(data.use, cells.1, cells.2, min.cells.group) {
  if (length(cells.1) == 0) {
    stop("Cell group 1 is empty - no cells with identity class")
  }
  if (length(cells.2) == 0) {
    stop("Cell group 2 is empty - no cells with identity class")
  }
  if (length(cells.1) < min.cells.group) {
    stop("Cell group 1 has fewer than ", min.cells.group, " cells")
  }
  if (length(cells.2) < min.cells.group) {
    stop("Cell group 2 has fewer than ", min.cells.group, " cells")
  }
  all.cells <- colnames(data.use)
  if (any(!cells.1 %in% all.cells)) {
    bad.cells <- all.cells[which(!as.character(cells.1) %in% all.cells)]
    stop(
      "The following cell names provided to cells.1 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  }
  if (any(!cells.2 %in% all.cells)) {
    bad.cells <- all.cells[which(!as.character(cells.2) %in% all.cells)]
    stop(
      "The following cell names provided to cells.2 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  }
  invisible(NULL)
}


