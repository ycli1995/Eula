#' @export
#' @concept data
DE.METHODS <- list(
  counts = c(),
  nocorrect = c("roc"),
  noprefilter = c("DESeq2"),
  latent = c("MAST")
)

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
  test.use = "wilcox",
  p.adjust.method = "bonferroni",
  p.thresh = 0.05,
  use.adjust = TRUE,
  latent.vars = NULL,
  min.cells.group = 3,
  max.cells.per.ident = Inf,
  seed = 42,
  densify = FALSE,
  verbose = TRUE,
  ...
) {
  if (!(test.use %in% DE.METHODS$latent) && !is.null(latent.vars)) {
    fastWarn(
      "'latent.vars' is only used for the following tests: ",
      paste(DE.METHODS$latent, collapse = ", "),
    )
  }
  validCellGroups(
    data.use = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    min.cells.group = min.cells.group
  )
  de.results0 <- data.frame(row.names = rownames(object))
  if (test.use == "roc") {
    de.results0$AUC <- 0.5
    de.results0$power <- 0
  } else {
    de.results0$p_val <- 1
    if (!test.use %in% DE.METHODS$nocorrect) {
      de.results0$p_val_adj <- 1
    }
  }
  features <- features %||% rownames(object)
  verboseMsg(
    "Test ", length(features), " of ", nrow(object), " features ",
    "with method '", test.use, "'"
  )
  if (length(features) == 0) {
    fastWarn("No feature tested. All features will be not significant.")
    de.results0$significance <- "nosig"
    return(de.results0)
  }

  # subsample cell groups if they are too large
  if (max.cells.per.ident < Inf) {
    verboseMsg("downsample cells to ", max.cells.per.ident)
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
    "wilcox" = differWilcox(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      ...
    ),
    "wilcox_limma" = differWilcox(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      limma = TRUE,
      ...
    ),
    "roc" = differROC(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2
    ),
    "t" = differTTest(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2
    ),
    "MAST" = differMAST(
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
  de.results0[rownames(de.results), colnames(de.results)] <- de.results
  if (test.use == "roc") {
    de.results0$significance <- ifelse(
      de.results0$AUC > 0.7 | de.results0$AUC < 0.3,
      yes = "sig",
      no = "nosig"
    )
    return(de.results0)
  }
  p.col <- "p_val"
  if (!test.use %in% DE.METHODS$nocorrect) {
    de.results0$p_val_adj <- p.adjust(
      p = de.results0$p_val,
      method = p.adjust.method
    )
    if (use.adjust) {
      p.col <- "p_val_adj"
    }
  }
  de.results0$significance <- ifelse(
    test = de.results0[[p.col]] < p.thresh,
    yes = "sig",
    no = "nosig"
  )
  de.results0
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
  p.thresh = 0.05,
  use.adjust = TRUE,
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
    p.thresh = p.thresh,
    use.adjust = use.adjust,
    latent.vars = latent.vars,
    max.cells.per.ident = max.cells.per.ident,
    seed = seed,
    ...
  )
}

#' @export differWilcox
differWilcox <- function(object, ...) {
  UseMethod("differWilcox", object)
}

#' @importFrom stats wilcox.test
#' @importFrom Eula.utils fastWarn
#'
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
    res <- res[seq_len(nrow(res) / 2), , drop = FALSE]
    p_val <- res$pval
    return(data.frame(p_val, row.names = rownames(object)))
  }
  if (limma.check) {
    limma <- TRUE
  }
  if (limma) {
    tmp <- object[1, ]
    overflow <- is.na(suppressWarnings(length(tmp) * length(tmp)))
    if (overflow) {
      fastWarn(
        "Use 'limma' with ", ncol(object), " cells may cause overflow. ",
        "Set 'limma = FALSE'."
      )
      limma <- FALSE
    }
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
  j <- seq_along(cells.1)
  p_val <- my.sapply(seq_len(nrow(object)), function(i) {
    min(1, 2 * min(limma::rankSumTestWithCorrelation(j, object[i, ])))
  })
  out <- data.frame(p_val, row.names = rownames(object))
  out
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
  out <- data.frame(p_val = unlist(p_val), row.names = rownames(object))
  out
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

#' @export
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
  object <- .check_subset_cells(mat = object, cells = c(cells.1, cells.2))
  latent.vars <- latent.vars %||% data.frame(row.names = colnames(object))
  group.info <- .deseq2_group_info(
    mat = object,
    cells.1 = cells.1,
    cells.2 = cells.2
  )
  colnames(group.info) <- c("condition")
  latent.vars.names <- c("condition", colnames(latent.vars))
  latent.vars <- cbind(latent.vars, group.info)
  fdat <- data.frame(primerid = rownames(object))
  rownames(fdat) <- fdat[, 1]
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(object),
    check_sanity = FALSE,
    cData = latent.vars,
    fData = fdat
  )
  cond <- factor(SummarizedExperiment::colData(sca)$condition)
  cond <- relevel(cond, ref = "Group1")
  SummarizedExperiment::colData(sca)$condition <- cond

  fmla <- as.formula(paste0(" ~ ", paste(latent.vars.names, collapse = "+")))
  zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
  summaryDt <- MAST::summary(zlmCond, doLRT = "conditionGroup2")$datatable
  summaryDt <- as.data.frame(summaryDt)
  summaryDt <- summaryDt[summaryDt[, "component"] == "H", , drop = FALSE]
  res <- data.frame(p_val = summaryDt[, 4], row.names = summaryDt[, 1])
  res
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
  differMAST.matrix(
    object = object,
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

#' @export
#' @method differROC CsparseMatrix
differROC.CsparseMatrix <- function(object, cells.1, cells.2, ...) {
  if (!requireNamespace("ROCR", quietly = TRUE)) {
    stop("Please install ROCR - `install.packages('ROCR')`")
  }
  my.lapply <- .get_lapply()

  object <- object[, c(cells.1, cells.2), drop = FALSE]
  myAUC <- unlist(my.lapply(rownames(object), function(g) {
    prediction.use <- ROCR::prediction(
      predictions = object[g, ],
      labels = c(rep(1, length(cells.1)), rep(0, length(cells.2))),
      label.ordering = 0:1
    )
    perf.use <- ROCR::performance(
      prediction.obj = prediction.use, 
      measure = "auc"
    )
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
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop(
      "Please install DESeq2 - learn more at ",
      "https://bioconductor.org/packages/release/bioc/html/DESeq2.html"
    )
  }
  res0 <- data.frame(row.names = rownames(object))
  res0$p_val <- 1
  countData <- .check_subset_cells(object, cells = c(cells.1, cells.2))
  group.info <- .deseq2_group_info(
    mat = countData, 
    cells.1 = cells.1, 
    cells.2 = cells.2
  )
  features <- rownames(object)[Matrix::rowSums(object > 0) == ncol(object)]
  countData <- .check_subset_features(countData, features = features)
  res <- .deseq2_matrix(
    countData = as.matrix(countData),
    group.info = group.info,
    ...
  )
  res0[rownames(res), "p_val"] <- res$p_val
  res0
}

#' @export
#' @method differDESeq2 CsparseMatrix
differDESeq2.CsparseMatrix <- function(object, cells.1, cells.2, ...) {
  differDESeq2.matrix(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    ...
  )
}

.check_subset_cells <- function(mat, cells) {
  if (setequal(colnames(mat), cells)) {
    return(mat)
  }
  mat <- mat[, cells, drop = FALSE]
  mat
}

.check_subset_features <- function(mat, features) {
  if (setequal(colnames(mat), features)) {
    return(mat)
  }
  mat <- mat[features, , drop = FALSE]
  mat
}

.deseq2_group_info <- function(mat, cells.1, cells.2) {
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(group.info[, "group"])
  group.info$wellKey <- rownames(group.info)
  group.info <- group.info[colnames(mat), , drop = FALSE]
  group.info
}

.deseq2_matrix <- function(countData, group.info, ...) {
  countData <- as.matrix(countData)
  dds1 <- DESeq2::DESeqDataSetFromMatrix(
    countData = countData,
    colData = group.info,
    design = ~group
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
  res <- data.frame(p_val = res$pvalue, row.names = rownames(res))
  res
}

.get_lapply <- function() {
  if (requireNamespace("future.apply", quietly = TRUE)) {
    if (future::nbrOfWorkers() > 1) {
      return(future.apply::future_lapply)
    }
  }
  lapply
}

.get_sapply <- function() {
  if (requireNamespace("future.apply", quietly = TRUE)) {
    if (future::nbrOfWorkers() > 1) {
      return(future.apply::future_sapply)
    }
  }
  sapply
}

.get_vapply <- function() {
  if (requireNamespace("future.apply", quietly = TRUE)) {
    if (future::nbrOfWorkers() > 1) {
      return(future.apply::future_vapply)
    }
  }
  vapply
}

#' @export
validCellGroups <- function(data.use, cells.1, cells.2, min.cells.group) {
  .valid_cell_groups(
    data.use = data.use,
    cell.groups = list(cells.1, cells.2),
    min.cells = min.cells.group
  )
}

.valid_cell_groups <- function(data.use, cell.groups, min.cells) {
  all.cells <- colnames(data.use)
  for (i in seq_along(cell.groups)) {
    if (length(cell.groups[[i]]) == 0) {
      stop("Cell group ", i, " is empty - no cells with identity class")
    }
    if (length(cell.groups[[i]]) < min.cells) {
      stop("Cell group ", i, " has fewer than ", min.cells, " cells")
    }
    if (any(!cell.groups[[i]] %in% all.cells)) {
      bad.cells <- head(setdiff(cell.groups[[i]], all.cells))
      stop(
        "The following cell names provided to cell group ", i,
        " are not present: ",
        paste(bad.cells, collapse = ", "),
        "..."
      )
    }
  }
  invisible(NULL)
}
