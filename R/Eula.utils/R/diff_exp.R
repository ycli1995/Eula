
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
    cells.1 = NULL,
    cells.2 = NULL,
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    max.cells.per.ident = Inf,
    min.cells.group = 3,
    logfc.threshold = 0.1,
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
  fc.results <- selectDE(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    mean.fxn = mean.fxn,
    logfc.threshold = logfc.threshold,
    min.pct = min.pct,
    min.diff.pct = min.diff.pct,
    min.cells.feature = min.cells.feature,
    only.pos = only.pos,
    min.exp = min.exp,
    pseudocount.use = pseudocount.use,
    base = base
  )
  if (nrow(fc.results) == 0) {
    if (test.use == "roc") {
      empty.df <- as.data.frame(matrix(nrow = 0, ncol = 2))
      colnames(empty.df) <- c("AUC", "power")
      return(cbind(fc.results, empty.df))
    }
    empty.df <- as.data.frame(matrix(nrow = 0, ncol = 1))
    colnames(empty.df) <- "p_val"
    return(cbind(fc.results, empty.df))
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
  # Actually perform the DE test
  de.results <- testDE(
    object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = rownames(fc.results),
    test.use = test.use,
    p.adjust.method = p.adjust.method,
    latent.vars = latent.vars,
    densify = densify,
    ...
  )
  cbind(fc.results, de.results[rownames(fc.results), , drop = FALSE])
}

#' @importFrom stats p.adjust
#' @export
#' @method differExp matrix
differExp.matrix <- function(
    object,
    cells.1 = NULL,
    cells.2 = NULL,
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    max.cells.per.ident = Inf,
    min.cells.group = 3,
    logfc.threshold = 0.1,
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
    object,
    cells.1,
    cells.2,
    mean.fxn = rowExpMean,
    logfc.threshold = 0.1,
    min.pct = 0.01,
    min.diff.pct = -Inf,
    min.cells.feature = 3,
    only.pos = FALSE,
    min.exp = 0,
    pseudocount.use = 1,
    base = 2,
    ...
) {
  fc.results <- foldChange(
    object,
    cells.1,
    cells.2,
    mean.fxn = mean.fxn,
    min.exp = min.exp,
    pseudocount.use = pseudocount.use,
    base = base,
    ...
  )
  if (only.pos) {
    selected <- fc.results[, 1] > logfc.threshold
  } else {
    selected <- abs(fc.results[, 1]) > logfc.threshold
  }
  ncells.min <- pmax(fc.results$ncells.1, fc.results$ncells.2)
  alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
  alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
  selected <- selected &
    (ncells.min > min.cells.feature) &
    (alpha.min > min.pct) &
    (alpha.diff > min.diff.pct)

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
    cells.1 = NULL,
    cells.2 = NULL,
    features = NULL,
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    latent.vars = NULL,
    densify = FALSE,
    ...
) {
  if (!(test.use %in% DE.METHODS$latent) && !is.null(latent.vars)) {
    fastWarning(
      "'latent.vars' is only used for the following tests: ",
      paste(DE.METHODS$latent, collapse = ", "),
    )
  }
  total.features <- nrow(object)
  features <- features %||% rownames(object)
  object <- object[features, c(cells.1, cells.2), drop = FALSE]
  if (!is.null(latent.vars)) {
    latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
  }
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
  if (test.use %in% DE.METHODS$nocorrect) {
    return(de.results)
  }
  de.results$p_val_adj <- p.adjust(
    p = de.results$p_val,
    method = p.adjust.method,
    n = total.features
  )
  de.results
}

#' @export
#' @method testDE matrix
testDE.matrix <- function(
    object,
    cells.1 = NULL,
    cells.2 = NULL,
    features = NULL,
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    latent.vars = NULL,
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
  if (!is.function(mean.fxn)) {
    stop("'mean.fxn' must be a function.")
  }
  ncells.1 <- rowSums(object[, cells.1, drop = FALSE] > min.exp)
  ncells.2 <- rowSums(object[, cells.2, drop = FALSE] > min.exp)
  pct.1 <- ncells.1 / length(cells.1)
  pct.2 <- ncells.2 / length(cells.2)
  mean.1 <- mean.fxn(object[, cells.1, drop = FALSE])
  mean.2 <- mean.fxn(object[, cells.2, drop = FALSE])

  p1 <- pseudocount.use / length(cells.1)
  p2 <- pseudocount.use / length(cells.2)
  fold.change <- log(mean.1 + p1, base = base) - log(mean.2 + p2, base = base)
  fc.results <- as.data.frame(
    fold_change = fold.change,
    mean.1 = mean.1,
    mean.2 = mean.2,
    ncells.1 = ncells.1,
    ncells.2 = ncells.2,
    pct.1 = pct.1,
    pct.2 = pct.2,
    row.names = rownames(object)
  )
  if (base == exp(1)) {
    base <- ""
  }
  fc.name <- paste0("avg_log", base, "FC")
  colnames(fc.results)[1] <- fc.name
  fc.results
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
  my.sapply <- .get_sapply()
  presto.check <- requireNamespace("presto", quietly = TRUE)
  limma.check <- requireNamespace("limma", quietly = TRUE)

  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(group.info[, "group"])
  object <- object[, rownames(group.info), drop = FALSE]

  if (presto.check[1] && (!limma)) {
    res <- presto::wilcoxauc(X = object, y = group.info[, "group"])
    res <- res[1:(nrow(res)/2), ]
    p_val <- res$pval
    return(data.frame(p_val, row.names = rownames(object)))
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
  my.sapply <- .get_sapply()
  p_val <- my.sapply(
    seq_len(nrow(object)),
    function(i) t.test(object[i, cells.1], object[i, cells.2], ...)$p.value
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
    cData = latent.vars,
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
  my.lapply <- .get_lapply()

  object <- object[, c(cells.1, cells.2), drop = FALSE]
  myAUC <- unlist(my.lapply(rownames(object), function(g) {
    prediction.use <- prediction(
      predictions = object[g, ],
      labels = c(rep(1, length(cells.1)), rep(0, length(cells.2))),
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
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(group.info[, "group"])
  group.info$wellKey <- rownames(group.info)
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


