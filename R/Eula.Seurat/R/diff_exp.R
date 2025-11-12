#' @importFrom Eula.utils DE.METHODS
NULL

#' @export
#' @method foldChange Assay
foldChange.Assay <- function(
    object,
    cells.1,
    cells.2,
    slot = "data",
    features = NULL,
    mean.fxn = NULL,
    min.exp = 0,
    pseudocount.use = 1,
    base = 2,
    ...
) {
  mean.fxn <- .select_mean_fxn(slot = slot, mean.fxn = mean.fxn)
  object <- .check_join_layers(object, slot = slot)

  data <- GetAssayData(object, slot)
  if (!is.null(features)) {
    data <- data[features, , drop = FALSE]
  }
  foldChange(
    object = data,
    cells.1 = cells.1,
    cells.2 = cells.2,
    mean.fxn = mean.fxn,
    min.exp = min.exp,
    pseudocount.use = pseudocount.use,
    base = base,
    ...
  )
}

#' @export
#' @method foldChange StdAssay
foldChange.StdAssay <- foldChange.Assay

#' @importClassesFrom SeuratObject DimReduc
#' @export
#' @method foldChange DimReduc
foldChange.DimReduc <- function(
    object,
    cells.1,
    cells.2,
    slot = NULL,
    features = NULL,
    mean.fxn = NULL,
    min.exp = 0,
    pseudocount.use = 1,
    base = 2,
    ...
) {
  mean.fxn <- mean.fxn %||% Matrix::rowMeans
  data <- t(Embeddings(object))
  if (!is.null(features)) {
    data <- data[features, , drop = FALSE]
  }
  foldChange(
    object = data,
    cells.1 = cells.1,
    cells.2 = cells.2,
    mean.fxn = mean.fxn,
    min.exp = min.exp,
    pseudocount.use = pseudocount.use,
    base = base,
    ...
  )
}

#' @export
#' @method foldChange Seurat
foldChange.Seurat <- function(
    object,
    ident.1 = NULL,
    ident.2 = NULL,
    group.by = NULL,
    slot = "data",
    features = NULL,
    assay = NULL,
    reduction = NULL,
    mean.fxn = NULL,
    min.exp = 0,
    pseudocount.use = 1,
    base = 2,
    ...
) {
  if (!is.null(assay) && !is.null(reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  if (!is.null(group.by)) {
    Idents(object) <- group.by
  }
  # select which data to use
  if (is.null(reduction)) {
    assay <- assay %||% DefaultAssay(object)
    data.use <- object[[assay]]
    cellnames.use <-  colnames(data.use)
  } else {
    data.use <- object[[reduction]]
    cellnames.use <- rownames(data.use)
  }
  cells <- .ident_to_cells_diff(
    object = object,
    ident.1 = ident.1,
    ident.2 = ident.2,
    cells = cellnames.use
  )
  foldChange(
    object = data.use,
    cells.1 = cells$cells.1,
    cells.2 = cells$cells.2,
    slot = slot,
    features = features,
    mean.fxn = mean.fxn,
    min.exp = min.exp,
    pseudocount.use = pseudocount.use,
    base = base,
    ...
  )
}

#' @importFrom SeuratObject WhichCells
.ident_to_cells <- function(object, ident, cells) {
  if (is.null(ident)) {
    stop("Please provide 'ident'.")
  }
  if (length(as.vector(ident)) > 1 && any(as.character(ident) %in% cells)) {
    bad.cells <- setdiff(cells, ident)
    if (length(bad.cells) > 0) {
      stop(
        "The following cells provided to 'ident' are not present:\n ",
        paste(bad.cells, collapse = ", ")
      )
    }
    return(ident)
  }
  WhichCells(object, idents = ident)
}

.ident_to_cells_diff <- function(object, ident.1, ident.2, cells) {
  cells.1 <- .ident_to_cells(object = object, ident = ident.1, cells = cells)
  if (is.null(ident.2)) {
    cells.2 <- setdiff(cells, cells.1)
  } else {
    cells.2 <- .ident_to_cells(object = object, ident = ident.2, cells = cells)
  }
  list(cells.1 = cells.1, cells.2 = cells.2)
}


#' @importFrom Eula.utils testDE
#' @export
#' @method testDE Assay
testDE.Assay <- function(
    object,
    cells.1,
    cells.2,
    slot = "data",
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
  data.slot <- ifelse(
    test = test.use %in% DE.METHODS$counts,
    yes = "counts",
    no = slot
  )
  object <- .check_join_layers(object, slot = data.slot)
  testDE(
    object = GetAssayData(object, slot),
    cells.1 = cells.1,
    cells.2 = cells.2,
    test.use = test.use,
    features = features,
    fc.results = fc.results,
    p.adjust.method = p.adjust.method,
    min.cells.group = min.cells.group,
    max.cells.per.ident = max.cells.per.ident,
    seed = seed,
    latent.vars = latent.vars,
    densify = densify,
    ...
  )
}

#' @export
#' @method testDE DimReduc
testDE.DimReduc <- function(
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
    ...
) {
  testDE(
    object = t(Embeddings(object)),
    cells.1 = cells.1,
    cells.2 = cells.2,
    test.use = test.use,
    features = features,
    fc.results = fc.results,
    p.adjust.method = p.adjust.method,
    min.cells.group = min.cells.group,
    max.cells.per.ident = max.cells.per.ident,
    seed = seed,
    latent.vars = latent.vars,
    densify = densify,
    ...
  )
}

#' @export
#' @method testDE Seurat
testDE.Seurat <- function(
    object,
    ident.1 = NULL,
    ident.2 = NULL,
    group.by = NULL,
    slot = "data",
    features = NULL,
    fc.results = NULL,
    assay = NULL,
    reduction = NULL,
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    latent.vars = NULL,
    min.cells.group = 3,
    max.cells.per.ident = Inf,
    seed = 42,
    densify = FALSE,
    ...
) {
  if (!is.null(assay) && !is.null(reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  if (!is.null(group.by)) {
    Idents(object) <- group.by
  }
  # select which data to use
  if (is.null(reduction)) {
    assay <- assay %||% DefaultAssay(object)
    data.use <- object[[assay]]
    cellnames.use <-  colnames(data.use)
  } else {
    data.use <- object[[reduction]]
    cellnames.use <- rownames(data.use)
  }
  cells <- .ident_to_cells_diff(
    object = object,
    ident.1 = ident.1,
    ident.2 = ident.2,
    cells = cellnames.use
  )
  # fetch latent.vars
  if (!is.null(latent.vars)) {
    latent.vars <- FetchData(
      object = object,
      vars = latent.vars,
      cells = c(cells$cells.1, cells$cells.2)
    )
  }
  testDE(
    object = data.use,
    cells.1 = cells$cells.1,
    cells.2 = cells$cells.2,
    test.use = test.use,
    features = features,
    fc.results = fc.results,
    p.adjust.method = p.adjust.method,
    min.cells.group = min.cells.group,
    max.cells.per.ident = max.cells.per.ident,
    seed = seed,
    latent.vars = latent.vars,
    densify = densify,
    ...
  )
}

.select_mean_fxn <- function(slot, mean.fxn = NULL) {
  new.mean.fxn <- switch(
    EXPR = slot,
    'data' = rowExpMean,
    'scale.data' = Matrix::rowMeans,
    'counts' = Matrix::rowMeans,
    rowExpMean
  )
  mean.fxn <- mean.fxn %||% new.mean.fxn
  mean.fxn
}

#' @importFrom Eula.utils selectDE
#' @export
#' @method testDE StdAssay
testDE.StdAssay <- testDE.Assay

#' @importFrom Eula.utils differExp
#' @export
#' @method differExp Assay
differExp.Assay <- function(
    object,
    cells.1,
    cells.2,
    slot = "data",
    fc.slot = "data",
    features = NULL,
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    min.cells.group = 3,
    max.cells.per.ident = Inf,
    logfc.threshold = 0.1,
    min.mean.exp = 0,
    min.pct = 0.01,
    min.diff.pct = -Inf,
    min.cells.feature = 3,
    only.pos = FALSE,
    mean.fxn = NULL,
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
    slot = fc.slot,
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
  testDE(
    object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    slot = slot,
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

#' @export
#' @method differExp DimReduc
differExp.DimReduc <- function(
    object,
    cells.1,
    cells.2,
    features = NULL,
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    min.cells.group = 3,
    max.cells.per.ident = Inf,
    logfc.threshold = 0.1,
    min.mean.exp = 0,
    min.pct = 0.01,
    min.diff.pct = -Inf,
    min.cells.feature = 3,
    only.pos = FALSE,
    mean.fxn = NULL,
    min.exp = 0,
    pseudocount.use = 1,
    base = 2,
    seed = 42,
    latent.vars = NULL,
    ...
) {
  if (test.use %in% DE.METHODS$counts) {
    stop(
      "The following tests cannot be used for differential expression on ",
      "a reduction as they assume a count model: ",
      paste(DE.METHODS$counts, collapse = ", ")
    )
  }
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
  mean.fxn <- mean.fxn %||% Matrix::rowMeans
  data <- t(Embeddings(object))
  if (!is.null(features)) {
    data <- data[features, , drop = FALSE]
  }
  fc.results <- foldChange(
    data,
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
    data,
    cells.1 = cells.1,
    cells.2 = cells.2,
    fc.results = fc.results,
    test.use = test.use,
    p.adjust.method = p.adjust.method,
    min.cells.group = min.cells.group,
    max.cells.per.ident = max.cells.per.ident,
    latent.vars = latent.vars,
    ...
  )
}

#' @export
#' @method differExp Seurat
differExp.Seurat <- function(
    object,
    ident.1 = NULL,
    ident.2 = NULL,
    group.by = NULL,
    slot = "data",
    features = NULL,
    assay = NULL,
    reduction = NULL,
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    latent.vars = NULL,
    min.cells.group = 3,
    max.cells.per.ident = Inf,
    logfc.threshold = 0.1,
    min.mean.exp = 0,
    min.pct = 0.01,
    min.diff.pct = -Inf,
    min.cells.feature = 3,
    only.pos = FALSE,
    mean.fxn = NULL,
    min.exp = 0,
    pseudocount.use = 1,
    base = 2,
    seed = 42,
    densify = FALSE,
    ...
) {
  if (!is.null(assay) && !is.null(reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  if (!is.null(group.by)) {
    Idents(object) <- group.by
  }
  # select which data to use
  if (is.null(reduction)) {
    assay <- assay %||% DefaultAssay(object)
    data.use <- object[[assay]]
    cellnames.use <-  colnames(data.use)
  } else {
    data.use <- object[[reduction]]
    cellnames.use <- rownames(data.use)
  }
  cells <- .ident_to_cells_diff(
    object = object,
    ident.1 = ident.1,
    ident.2 = ident.2,
    cells = cellnames.use
  )
  # fetch latent.vars
  if (!is.null(latent.vars)) {
    latent.vars <- FetchData(
      object = object,
      vars = latent.vars,
      cells = c(cells$cells.1, cells$cells.2)
    )
  }
  .validate_cell_groups(
    data.use = object,
    cells.1 = cells$cells.1,
    cells.2 = cells$cells.2,
    min.cells.group = min.cells.group
  )
  fc.results <- foldChange(
    object = data.use,
    cells.1 = cells$cells.1,
    cells.2 = cells$cells.2,
    slot = slot,
    features = features,
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

  # fetch latent.vars
  if (!is.null(latent.vars)) {
    latent.vars <- FetchData(
      object = object,
      vars = latent.vars,
      cells = c(cells$cells.1, cells$cells.2)
    )
  }
  testDE(
    object = data.use,
    cells.1 = cells$cells.1,
    cells.2 = cells$cells.2,
    test.use = test.use,
    fc.results = fc.results,
    p.adjust.method = p.adjust.method,
    min.cells.group = min.cells.group,
    max.cells.per.ident = max.cells.per.ident,
    seed = seed,
    latent.vars = latent.vars,
    densify = densify,
    ...
  )
}

#' @export
findAllMarkers <- function(
    object,
    assay = NULL,
    features = NULL,
    group.by = NULL,
    slot = 'data',

    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    latent.vars = NULL,
    min.cells.group = 3,
    max.cells.per.ident = Inf,
    return.thresh = 1e-2,
    use.adjust = TRUE,

    logfc.threshold = 0.1,
    min.mean.exp = 0,
    min.pct = 0.01,
    min.diff.pct = -Inf,
    min.cells.feature = 3,
    only.pos = FALSE,

    mean.fxn = NULL,
    min.exp = 0,
    pseudocount.use = 1,
    base = 2,
    seed = 42,
    densify = FALSE,
    verbose = TRUE,
    ...
) {
  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh <- 0.7
  }
  if (!is.null(group.by) && !identical(group.by, "ident")) {
    if (length(group.by) == 1 && ! group.by %in% colnames(object@meta.data)) {
      stop("'", group.by, "' not found in object metadata")
    }
    Idents(object) <- group.by
  }
  idents.all <- sort(unique(Idents(object)))

  genes.de <- list()
  messages <- list()
  for (i in 1:length(x = idents.all)) {
    if (verbose) {
      message("Calculating cluster ", idents.all[i])
    }
    genes.de[[i]] <- tryCatch(
      expr = {
        differExp(
          object = object,
          ident.1 = idents.all[i],
          ident.2 = NULL,
          slot = slot,
          features = features,
          assay = assay,
          test.use = test.use,
          p.adjust.method = p.adjust.method,
          latent.vars = latent.vars,
          min.cells.group = min.cells.group,
          max.cells.per.ident = max.cells.per.ident,
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
          densify = densify,
          ...
        )
      },
      error = function(cond) return(cond$message)
    )
    if (is.character(genes.de[[i]])) {
      messages[[i]] <- genes.de[[i]]
      genes.de[[i]] <- NULL
    }
  }
  gde.all <- data.frame()
  p.col <- "p_val"
  if (use.adjust) {
    p.col <- "p_val_adj"
  }
  if (test.use == "roc") {
    p.col <- "AUC"
  }
  for (i in 1:length(idents.all)) {
    if (is.null(unlist(genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(gde) > 0) {
      if (test.use == "roc") {
        idx <- (gde$AUC > return.thresh | gde$AUC < (1 - return.thresh))
      } else {
        idx <- gde[[p.col]] < return.thresh
      }
      gde <- gde[idx, , drop = FALSE]
      gde <- gde[order(gde[[p.col]], -abs(gde$pct.1-gde$pct.2)), , drop = FALSE]

      if (nrow(gde) > 0) {
        gde$cluster <- idents.all[i]
        gde$gene <- rownames(gde)
        gde.all <- rbind(gde.all, gde)
      }
    }
  }
  if (nrow(gde.all) == 0) {
    fastWarning("No DE genes identified")
  }
  rownames(gde.all) <- make.unique(as.character(gde.all$gene))
  if (length(messages) == 0) {
    return(gde.all)
  }
  fastWarning("The following tests were not performed: ")
  for (i in seq_along(messages)) {
    if (!is.null(messages[[i]])) {
      fastWarning("When testing ", idents.all[i], ":\n\t", messages[[i]])
    }
  }
  gde.all
}

#' @importFrom dplyr arrange desc slice_head starts_with
#' @export
getTopMarkers <- function(markers, top.n = 5, group.by = "cluster", ...) {
  p.cols <- intersect(colnames(markers), c("p_val", "p_val_adj", "AUC"))[1]
  logfc.col <- grep("^avg_log", colnames(markers), value = TRUE)[1]
  top <- markers %>%
    group_by(across(all_of(group.by))) %>%
    dplyr::arrange(desc(!!logfc.col), !!p.cols) %>%
    dplyr::slice_head(n = top.n)
  top
}

.check_join_layers <- function(object, slot) {
  if (packageVersion("Seurat") >= package_version("5.0.0")) {
    if (length(SeuratObject::Layers(object, search = slot)) > 1) {
      fastWarning(slot, " layers are not joined. Running `JoinLayers`.")
      object <- JoinLayers(object)
    }
  }
  object
}
