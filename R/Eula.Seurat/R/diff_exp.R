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
    base = -1,
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
  object <- .check_group_by(object, group.by = group.by)
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
    bad.cells <- setdiff(ident, cells)
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
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    p.thresh = 0.05,
    use.adjust = TRUE,
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
    p.adjust.method = p.adjust.method,
    p.thresh = p.thresh,
    use.adjust = use.adjust,
    min.cells.group = min.cells.group,
    max.cells.per.ident = max.cells.per.ident,
    seed = seed,
    latent.vars = latent.vars,
    densify = densify,
    ...
  )
}

#' @importFrom Eula.utils selectDE
#' @export
#' @method testDE StdAssay
testDE.StdAssay <- testDE.Assay

#' @export
#' @method testDE DimReduc
testDE.DimReduc <- function(
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
    ...
) {
  testDE(
    object = t(Embeddings(object)),
    cells.1 = cells.1,
    cells.2 = cells.2,
    test.use = test.use,
    features = features,
    p.adjust.method = p.adjust.method,
    p.thresh = p.thresh,
    use.adjust = use.adjust,
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
    assay = NULL,
    reduction = NULL,
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    p.thresh = 0.05,
    use.adjust = TRUE,
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
    p.adjust.method = p.adjust.method,
    p.thresh = p.thresh,
    use.adjust = use.adjust,
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
    p.thresh = 0.05,
    use.adjust = TRUE,
    filter.nosig = TRUE,
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
  if (test.use %in% DE.METHODS$counts) {
    slot <- "counts"
  }
  mean.fxn <- .select_mean_fxn(slot = fc.slot, mean.fxn = mean.fxn)
  data <- GetAssayData(object, slot)
  fc.data <- GetAssayData(object, fc.slot)
  features <- features %||% rownames(data)
  if (!setequal(features, rownames(data))) {
    data <- data[features, , drop = FALSE]
  }
  if (!setequal(features, rownames(fc.data))) {
    fc.data <- fc.data[features, , drop = FALSE]
  }
  differExp(
    object = data,
    cells.1 = cells.1,
    cells.2 = cells.2,
    fc.data = fc.data,
    test.use = test.use,
    p.adjust.method = p.adjust.method,
    p.thresh = p.thresh,
    use.adjust = use.adjust,
    filter.nosig = filter.nosig,
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
    latent.vars = latent.vars,
    densify = densify,
    ...
  )
}

#' @export
#' @method differExp StdAssay
differExp.StdAssay <- differExp.Assay

#' @export
#' @method differExp DimReduc
differExp.DimReduc <- function(
    object,
    cells.1,
    cells.2,
    features = NULL,
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    p.thresh = 0.05,
    use.adjust = TRUE,
    filter.nosig = TRUE,
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
    base = -1,
    seed = 42,
    latent.vars = NULL,
    densify = FALSE,
    ...
) {
  if (test.use %in% DE.METHODS$counts) {
    stop(
      "The following tests cannot be used for differential expression on ",
      "a reduction as they assume a count model: ",
      paste(DE.METHODS$counts, collapse = ", ")
    )
  }
  data <- t(Embeddings(object))
  .validate_cell_groups(
    data.use = data,
    cells.1 = cells.1,
    cells.2 = cells.2,
    min.cells.group = min.cells.group
  )
  if (test.use %in% DE.METHODS$noprefilter) {
    min.diff.pct <- -Inf
    logfc.threshold <- 0
  }
  mean.fxn <- mean.fxn %||% Matrix::rowMeans
  features <- features %||% rownames(data)
  if (!setequal(features, rownames(data))) {
    data <- data[features, , drop = FALSE]
  }
  de.results <- differExp(
    object = data,
    cells.1 = cells.1,
    cells.2 = cells.2,
    fc.data = data,
    test.use = test.use,
    p.adjust.method = p.adjust.method,
    p.thresh = p.thresh,
    use.adjust = use.adjust,
    filter.nosig = filter.nosig,
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
    latent.vars = latent.vars,
    densify = FALSE,
    ...
  )
  if (test.use == "roc") {
    return(de.results[order(-de.results$power, -de.results[[1]]), ])
  }
  de.results[order(-de.results$p_val, -de.results[[1]]), ]
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
    p.thresh = 0.05,
    use.adjust = TRUE,
    filter.nosig = TRUE,
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
  object <- .check_group_by(object, group.by = group.by)
  # select which data to use
  if (is.null(reduction)) {
    assay <- assay %||% DefaultAssay(object)
    data.use <- object[[assay]]
    cellnames.use <-  colnames(data.use)
  } else {
    data.use <- object[[reduction]]
    cellnames.use <- rownames(data.use)
    base <- -1
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
  differExp(
    object = data.use,
    cells.1 = cells$cells.1,
    cells.2 = cells$cells.2,
    slot = slot,
    fc.slot = slot,
    features = features,
    test.use = test.use,
    p.adjust.method = p.adjust.method,
    p.thresh = p.thresh,
    use.adjust = use.adjust,
    filter.nosig = filter.nosig,
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
    p.thresh = 1e-2,
    use.adjust = TRUE,
    filter.nosig = TRUE,

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
  object <- .check_group_by(object, group.by = group.by)
  idents.all <- sort(unique(Idents(object)))

  gde.all <- list()
  for (i in seq_along(idents.all)) {
    ident <- idents.all[i]
    if (verbose) {
      message("Calculating cluster ", ident)
    }
    gde <- tryCatch(
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
          p.thresh = p.thresh,
          use.adjust = use.adjust,
          filter.nosig = filter.nosig,
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
    if (is.character(gde)) {
      fastWarning("Testing ", ident, " failed:\n\t", gde)
      next
    }
    if (nrow(gde) == 0) {
      next
    }
    gde$cluster <- ident
    gde$gene <- rownames(gde)
    rownames(gde) <- NULL
    gde.all[[ident]] <- gde
  }
  gde.all <- Reduce(rbind, gde.all)
  gde.all
}

#' @export
findGroupDiffer <- function(
    object,
    group.data,
    differs,

    group.by = NULL,
    assay = NULL,
    reduction = NULL,
    features = NULL,
    slot = 'data',

    diff.type = c("all", "only_bulk", "no_bulk"),

    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    p.thresh = 1e-2,
    use.adjust = TRUE,
    filter.nosig = FALSE,
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
    verbose = TRUE,
    ...
) {
  if (!is.list(differs)) {
    stop("'differs' must be a list containing comparison information.")
  }

  diff.type <- match.arg(diff.type)
  is.bulk <- diff.type %in% c("all", "only_bulk")
  is.cluster <- diff.type %in% c("all", "no_bulk")

  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh <- 0.7
  }
  p.col <- "p_val"
  if (use.adjust) {
    p.col <- "p_val_adj"
  }
  if (test.use == "roc") {
    p.col <- "AUC"
  }

  if (!is.null(assay) && !is.null(reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  if (is.null(reduction)) {
    assay <- assay %||% DefaultAssay(object)
    DefaultAssay(object) <- assay
    data.use <- object[[assay]]
    cellnames.use <-  colnames(data.use)
  } else {
    data.use <- object[[reduction]]
    cellnames.use <- rownames(data.use)
    base <- -1
  }

  object <- .check_group_by(object, group.by = group.by)
  cells.cluster <- list()
  if (is.bulk) {
    cells.cluster[["BULK"]] <- cellnames.use
  }
  if (is.cluster) {
    cells.cluster <- c(cells.cluster, split(cellnames.use, f = Idents(object)))
  }

  markers <- list()
  for (i in seq_along(differs)) {
    diff <- differs[[i]]
    if (any(!diff %in% colnames(group.data))) {
      bad.groups <- paste(setdiff(diff, colnames(group.data)), collapse = ", ")
      stop("The following groups are not found:\n ", bad.groups)
    }
    gde.all <- list()
    for (cluster in names(cells.cluster)) {
      cells <- cells.cluster[[cluster]]
      cells.1 <- intersect(cells, rownames(group.data)[group.data[, diff[[1]]]])
      cells.2 <- intersect(cells, rownames(group.data)[group.data[, diff[[2]]]])

      if (verbose) {
        message("Calculating cluster ", cluster)
      }
      gde <- tryCatch(
        expr = {
          differExp(
            object = data.use,
            cells.1 = cells.1,
            cells.2 = cells.2,
            slot = slot,
            features = features,
            assay = assay,
            test.use = test.use,
            p.adjust.method = p.adjust.method,
            p.thresh = p.thresh,
            use.adjust = use.adjust,
            filter.nosig = filter.nosig,
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
      if (is.character(gde)) {
        fastWarning("Testing ", cluster, " failed:\n\t", gde)
        next
      }
      if (nrow(gde) == 0) {
        next
      }
      gde$cluster <- cluster
      gde$gene <- rownames(gde)
      rownames(gde) <- NULL
      gde.all[[cluster]] <- gde
    }
    gde.all <- Reduce(rbind, gde.all)
    diff <- paste(diff, collapse = "-vs-")
    markers[[diff]] <- gde.all
  }
  markers
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

.check_group_by <- function(object, group.by = NULL) {
  if (!is.null(group.by) && !identical(group.by, "ident")) {
    group.by <- group.by[1]
    if (length(group.by) == 1 && ! group.by %in% colnames(object@meta.data)) {
      stop("'", group.by, "' not found in object metadata")
    }
    Idents(object) <- group.by
  }
  object
}


