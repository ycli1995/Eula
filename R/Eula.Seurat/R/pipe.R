#' @importFrom Eula.utils captureMsg getArgList getDefaultArgs pipeMsg
NULL

#' @export
pipe_RenameObject <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "RenameObject")
  pipeMsg("Start")

  if (length(params[["RenameObject_yaml"]]) > 0) {
    file <- normalizePath(params[["RenameObject_yaml"]], mustWork = TRUE)
    pipeMsg("Read 'RenameObject_yaml': ", file)
    params <- readYAML(file2)
  }

  obj <- pipe_AddMetaData(
    obj,
    params = params[['AddMetaData']],
    pipe.name = pipe.name
  )

  pipeMsg('Rename columns in meta.data', pipe.name = pipe.name)
  obj <- RenameSeuratColumns(obj, params[['rename_colnames']])

  obj <- pipe_RenameSeurat(obj, params[['rename']], pipe.name = pipe.name)

  obj <- pipe_UpdateSeurat(obj, params[['UpdateSeurat']], pipe.name = pipe.name)

  pipeMsg('Check Seurat extra information', pipe.name = pipe.name)
  obj <- CheckMySeuratObj(obj, column.map = params[['default_colnames']])

  obj <- pipe_SubsetObject(obj, params[['subset']], pipe.name = pipe.name)

  obj <- pipe_ShrinkSeuratObject(obj, params, pipe.name = pipe.name)

  pipeMsg('Finally check Seurat object')
  captureMsg(print(obj))
  captureMsg(str(obj@meta.data))
  captureMsg(print(obj@misc$colors))

  obj
}

#' @export
pipe_RenameSeuratMetaData <- function(
    obj,
    params = list(),
    pipe.name = NULL,
    ...
) {
  pipe.name <- c(pipe.name, "RenameSeuratMetaData")
  pipeMsg("Start")

  COLOR.NAMES <- list(
    "orig.ident" = "color.sample",
    "Groups" = "color.group",
    "seurat_clusters" = "color.cluster",
    "Samples" = "color.sample",
    "Cluster" = "color.cluster",
    "Clusters" = "color.cluster"
  )

  params <- fetch_rename_table(params)
  captureMsg(str(params))

  column <- intersect(params[["column"]], colnames(obj[[]]))
  if (length(column) == 0) {
    fastWarning("Cannot found 'column' name. No rename will do.")
    return(obj)
  }

  set.ident <- params[["set.ident"]] %||% FALSE
  keep.orders <- params[["keep.orders"]] %||% FALSE

  color.name <- params[["color_name"]] %||% COLOR.NAMES[[column]]
  color.name <- color.name %||% column

  params[["rename_map"]] <- params[["rename_map"]] %||%
    if (is.factor(obj@meta.data[, column])) {
      obj@meta.data[, column] %>%
        levels() %>%
        sapply(function(x) x, simplify = FALSE)
    } else {
      obj@meta.data[, column] %>%
        factor() %>%
        levels() %>%
        sapply(function(x) x, simplify = FALSE)
    }

  colors <- unlist(params[["colors"]])
  obj <- RenameSeurat(
    object = obj,
    map = params[["rename_map"]],
    column = column,
    misc.name = color.name,
    colors = colors,
    set.ident = set.ident,
    keep.orders = keep.orders
  )

  obj
}

#' @export
pipe_RenameSeurat <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "RenameSeurat")
  pipeMsg("Start")

  captureMsg(str(obj[[]]))
  captureMsg(str(obj@misc))

  METADATA.COLUMNS <- list(
    "sample" = "orig.ident",
    "group" = "Groups",
    "cluster" = "seurat_clusters"
  )

  for (i in names(params)) {
    pipeMsg("Rename '", i, "': ")
    params[[i]][["column"]] <- params[[i]][["column"]] %||%
      params[[i]]
    obj <- pipe_RenameSeuratMetaData(obj, params[[i]], pipe.name = pipe.name)
  }

  obj
}

#' @importFrom SeuratObject AddMetaData
#' @export
pipe_AddMetaData <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "AddMetaData")
  pipeMsg("Start")

  defaults <- list(infile = NULL)
  params <- getDefaultArgs(defaults, params)
  captureMsg(str(params))
  list2env(params, envir = environment())

  for (i in seq_along(infile)) {
    df <- readTable(infile[[i]], row.names = 1, keep.row.names = FALSE)
    captureMsg(str(df))
    obj <- AddMetaData(obj, df)
  }

  obj
}

#' @export
pipe_ShrinkSeuratObject <- function(
    obj,
    params = list(),
    pipe.name = NULL,
    ...
) {
  pipe.name <- c(pipe.name, "ShrinkSeuratObject")
  pipeMsg("Start")

  defaults <- list(
    scale.data = FALSE,
    misc.counts = TRUE,
    assays = Assays(obj)
  )
  params <- getDefaultArgs(defaults, params)
  params$assays <- intersect(params$assays, Assays(obj))
  params <- params[names(defaults)]

  captureMsg(str(params))
  list2env(params, envir = environment())

  obj <- ShrinkSeuratObject(
    object = obj,
    assays = assays,
    scale.data = scale.data,
    misc.counts = misc.counts
  )

  obj
}

#' @export
pipe_SubsetObject <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "SubsetObject")
  pipeMsg('Start')

  params$cells <- params$cells %||% params$cells_use
  params$features <- params$features %||% params$features_use
  captureMsg(str(params))

  cells_features <- get_cells_and_features(
    object = obj,
    cells = params$cells,
    features = params$features,
    cells_exclude = params$cells_exclude,
    features_exclude = params$features_exclude
  )

  params <- params[intersect(names(params), colnames(obj[[]]))]
  obj <- SubsetObject(
    object = obj,
    column_map = params,
    cells = cells_features[['cells']],
    features = cells_features[['features']]
  )

  obj
}

#' @export
pipe_UpdateSeurat <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "UpdateSeurat")
  pipeMsg('Start')

  params$Assay5 <- params$Assay5 %||% FALSE
  captureMsg(str(params))

  obj <- UpdateSeuratAll(object = obj, Assay5 = params$Assay5)

  obj
}

#' @export
pipe_DoubletFinder <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "DoubletFinder")
  pipeMsg('Start')
  checkPackages("DoubletFinder")

  defaults <- list(
    outdir = getwd(),
    pN = 0.25,
    dims = 1:50,
    rate = NULL
  )
  params <- getDefaultArgs(defaults, params)
  captureMsg(str(params))
  list2env(params, envir = environment())

  obj <- pipe_NormalizeData(obj, params = list(), pipe.name = pipe.name)
  obj <- pipe_FindVariableFeatures(obj, params = list(), pipe.name = pipe.name)
  obj <- pipe_ScaleData(obj, params = list(), pipe.name = pipe.name)
  obj <- pipe_RunPCA(
    obj,
    params = list(npcs = max(dims)),
    pipe.name = pipe.name
  )
  obj <- pipe_RunTSNE(obj, params = list(dims = dims), pipe.name = pipe.name)
  obj <- pipe_RunUMAP(obj, params = list(dims = dims), pipe.name = pipe.name)

  sweep.res.list <- DoubletFinder::paramSweep(obj, PCs = dims, sct = FALSE)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list)

  pdf("find_pK.pdf", width = 5, height = 4)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  dev.off()
  pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)[1]]))

  rate <- rate %||% 7.6 * 10^-6 * ncol(obj) + 5.27 * 10^-4

  nExp_poi <- round(as.numeric(rate) * ncol(obj))

  if (exists("seurat_clusters", obj@meta.data)) {
    annotations <- obj@meta.data$seurat_clusters
    homotypic.prop <- DoubletFinder::modelHomotypic(annotations)
    nExp_poi <- round(nExp_poi * (1 - homotypic.prop))
  }

  obj <- DoubletFinder::doubletFinder(
    obj,
    PCs = dims,
    pN = pN,
    pK = pK,
    nExp = nExp_poi,
    reuse.pANN = NULL,
    sct = FALSE
  )
  colnames(obj@meta.data)[grep('pANN', colnames(obj@meta.data))] <- "pANN"
  colnames(obj@meta.data)[grep('DF.classifications', colnames(obj@meta.data))] <- "classifications"

  data <- FetchSeuratData(obj, c("pANN", "classifications"))
  writeTable(data, "DF.classify.xls")

  colors <- c("Singlet" = "black", "Doublet" = "red")

  pipeMsg('Plotting dim_plot')
  group.by <- c("classifications")
  reductions <- c("tsne", "umap")
  corner.axis <- TRUE
  save_dim_plots(
    obj,
    outdir = outdir,
    reductions = reductions,
    group.by = group.by,
    corner.axis = corner.axis,
    colors = colors
  )
  pipeMsg('Plotting feature_dim_plot')
  save_feature_dim_plots(
    obj,
    features = "pANN",
    outdir = outdir,
    reductions = reductions,
    combine = FALSE,
    corner.axis = corner.axis
  )

  obj
}

#' @export
pipe_StatFeaturePercentage <- function(
    obj,
    params = list(),
    pipe.name = NULL,
    ...
) {
  pipe.name <- c(pipe.name, "UpdateSeurat")
  pipeMsg('Start')

  stat.pct <- params[['stat.pct']] %||% TRUE
  captureMsg(str(params))

  params[['stat.pct']] <- NULL
  for (i in names(params)) {
    i2 <- gsub("\\/| ", "_", i)
    obj <- StatFeatures(
      object = obj,
      features = params[[i]],
      col.name = paste0("percent.", i2),
      stat.pct = stat.pct
    )
  }

  obj
}

#' @export
pipe_StatCellQuality <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "StatCellQuality")
  pipeMsg("Start")

  defaults <- list(
    outdir = getwd(),
    sample.by = "orig.ident",
    group.by = "Groups",
    features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.rrna")
  )
  params <- getDefaultArgs(defaults, params)
  captureMsg(str(params))
  list2env(params, envir = environment())

  features.use <- intersect(features, colnames(obj[[]]))
  if (length(features.use) == 0) {
    fastWarning(
      "The following features are not found in obj@meta.data:\n ",
      paste(features, collapse = ", ")
    )
    return(invisible(NULL))
  }
  sample.by <- sample.by[1]
  if (!sample.by %in% colnames(obj[[]])) {
    stop("'", sample.by, "' not found in obj@meta.data.")
  }
  group.by <- group.by[1]
  if (!group.by %in% colnames(obj[[]])) {
    group.by <- sample.by
  }
  pipeMsg("Ploting violin plots")
  mkdir(outdir)
  save_violin_plots(
    obj = obj,
    features = features.use,
    outdir = outdir,
    basic.size = 4,
    group.by = sample.by,
    combine = TRUE,
    pt.size = 0.1,
    raster = TRUE,
    raster.dpi = 300,
    ncol = length(features.use),
    split.by = NULL
  )
  save_violin_plots(
    obj = obj,
    features = features.use,
    outdir = outdir,
    basic.size = 4,
    group.by = group.by,
    combine = TRUE,
    pt.size = 0.1,
    raster = TRUE,
    raster.dpi = 300,
    ncol = length(features.use),
    split.by = NULL
  )
  invisible(NULL)
}

#' @export
pipe_MakeSeuratObj <- function(params, pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "MakeSeuratObject")
  pipeMsg("Start")
  captureMsg(print(params))

  if (length(params) == 0) {
    stop("'MakeSeuratObj' params is empty.")
  }
  data.paths <- params[["paths"]]
  if (length(data.paths) == 0) {
    stop("No 'paths' for matrices.")
  }
  data.names <- params[["names"]]
  if (length(data.names) != length(data.paths)) {
    stop("'names' must has the same length as 'paths'")
  }
  name.list <- params[['name_list']]
  if (length(name.list) == 0) {
    fastWarning("No 'name_list' for MakeSeuratObj")
  }
  assay <- params[['assay']] %||% "RNA"

  obj <- MakeSeuratObj(
    data.paths = data.paths,
    data.names = data.names,
    assay = assay,
    name.list = name.list
  )

  pipeMsg("Checking object created", pipe.name = pipe.name)
  captureMsg(print(obj))

  pipeMsg("Checking misc data", pipe.name = pipe.name)
  captureMsg(str(obj@misc))

  obj
}

#' @importFrom Seurat FindVariableFeatures
#' @importFrom SeuratObject VariableFeatures VariableFeatures<-
#' @export
pipe_FindVariableFeatures <- function(
    obj,
    params = list(),
    pipe.name = NULL,
    ...
) {
  pipe.name <- c(pipe.name, "FindVariableFeatures")
  pipeMsg("Start")

  defaults <- list(
    assay = DefaultAssay(obj),
    selection.method = "vst",
    nfeatures = 2000
  )
  params <- getDefaultArgs(defaults, params)
  params[['vfeatures']] <- getArgList(params[['vfeatures']])
  params[['vfeature.must']] <- getArgList(params[['vfeature.must']])
  params[['vfeature.remove']] <- getArgList(params[['vfeature.remove']])
  captureMsg(str(params))

  if (length(params[['vfeatures']]) > 0) {
    pipeMsg("Use preset variable features")
    list2env(params, envir = environment())

    vfeatures <- getFeaturesID(obj, features = vfeatures)
    VariableFeatures(obj[[assay]]) <- vfeatures
    return(obj)
  }

  pipeMsg("Find variable features")
  list2env(params, envir = environment())
  obj <- FindVariableFeatures(
    object = obj,
    assay = assay,
    selection.method = selection.method,
    nfeatures = nfeatures,
    ...
  )

  if (length(params[['vfeature.must']]) > 0) {
    pipeMsg("Check vfeatures required")
    vfeature.must <- getFeaturesID(obj, features = vfeature.must)
    VariableFeatures(obj) <- union(VariableFeatures(obj), vfeature.must)
  }

  if (length(params[['vfeature.remove']]) > 0) {
    pipeMsg("Check vfeatures excluded")
    vfeature.remove <- getFeaturesID(obj, features = vfeature.remove)
    VariableFeatures(obj) <- setdiff(VariableFeatures(obj), vfeature.remove)
  }


  obj
}

#' @export
pipe_CellCycleScoring <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "CellCycleScoring")
  pipeMsg("Start")

  defaults <- list(
    s.features = Seurat::cc.genes$s.genes,
    g2m.features = Seurat::cc.genes$g2m.genes,
    regress = "diff",
    assay = DefaultAssay(obj),
    slot = "data",
    set.ident = FALSE
  )

  params[['s.features']] <- getArgList(params[['s.features']])
  params[['g2m.features']] <- getArgList(params[['g2m.features']])
  params <- getDefaultArgs(defaults, params)
  captureMsg(str(params))

  if (any(lengths(params[c('s.features', 'g2m.features')]) < 2)) {
    print(str(params))
    fastWarning("No enough cell cycle gene (n < 2). Skip cell cycle scoring.")
    return(obj)
  }

  params[['s.features']] <- getFeaturesID(obj, params[['s.features']])
  params[['g2m.features']] <- getFeaturesID(obj, params[['g2m.features']])
  list2env(params, envir = environment())
  ctrl <- params[['ctrl']]
  regress <- match.arg(regress, c("none", "diff", "all"))

  obj <- CellCycle(
    object = obj,
    s.features = s.features,
    g2m.features = g2m.features,
    set.ident = set.ident,
    slot = slot,
    assay = assay,
    ctrl = ctrl,
    ...
  )
  if (regress == "none") {
    return(obj)
  }
  if (regress == "diff") {
    obj@misc$vars.to.regress <- c("CC.Difference", obj@misc$vars.to.regress)
    return(obj)
  }
  obj@misc$vars.to.regress <- c(
    "S.Score", "G2M.Score",
    obj@misc$vars.to.regress
  )
  obj
}

#' @importFrom Seurat NormalizeData
#' @export
pipe_NormalizeData <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "NormalizeData")
  pipeMsg("Start")

  defaults <- list(
    method = "LogNormalize",
    scale.factor = 10000,
    assay = DefaultAssay(obj)
  )
  params <- getDefaultArgs(defaults, params)
  captureMsg(str(params))

  list2env(params, envir = environment())
  obj <- NormalizeData(
    object = obj,
    assay = assay,
    normalization.method = method,
    scale.factor = scale.factor,
    verbose = FALSE,
    ...
  )
  obj
}

#' @importFrom Seurat ScaleData
#' @export
pipe_ScaleData <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "ScaleData")
  pipeMsg("Start")

  defaults <- list(
    assay = DefaultAssay(obj),
    scale.max = 10,
    only.var.features = TRUE
  )
  params <- getDefaultArgs(defaults, params)
  params[['vars.to.regress']] <- params[['vars.to.regress']] %||%
    c(paste0("nCount_", params[['assay']]), "percent.mito")
  params[['vars.to.regress']] <- unique(c(
    params[['vars.to.regress']],
    obj@misc$vars.to.regress
  ))
  captureMsg(str(params))
  list2env(params, envir = environment())

  vars.to.regress <- intersect(vars.to.regress, colnames(obj[[]]))
  features <- if (only.var.features) NULL else rownames(obj)

  obj <- ScaleData(
    object = obj,
    assay = assay,
    vars.to.regress = vars.to.regress,
    features = features,
    scale.max = scale.max,
    verbose = TRUE,
    ...
  )
  obj
}

#' @importFrom Seurat RunPCA
#' @export
pipe_RunPCA <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "RunPCA")
  pipeMsg("Start")

  defaults <- list(
    assay = DefaultAssay(obj),
    npcs = 50,
    rev.pca = FALSE,
    weight.by.var = TRUE,
    seed = 42,
    approx = TRUE,
    reduction.key = "PC_",
    reduction.name = "pca"
  )
  params <- getDefaultArgs(defaults, params)
  params[['reduction.surfix']] <- params[['reduction.surfix']] %||%
    params[['assay']]
  captureMsg(str(params))
  list2env(params, envir = environment())

  obj <- RunPCA(
    object = obj,
    assay = assay,
    npcs = npcs,
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    verbose = FALSE,
    reduction.key = reduction.key,
    reduction.name = reduction.name,
    seed.use = seed,
    approx = approx,
    ...
  )
  obj[[paste0(reduction.name, "_", reduction.surfix)]] <- obj[[reduction.name]]
  obj@misc$reduction <- list(
    RunUMAP = reduction.name,
    RunTSNE = reduction.name,
    FindNeighbors = reduction.name
  )
  obj
}

#' @importFrom SeuratObject Embeddings
#' @importFrom Seurat RunTSNE
#' @export
pipe_RunTSNE <- function(obj, params = list(), pipe.name = NULL, ...) {
  run <- params[['run']] %||% TRUE
  if (!run) {
    fastWarning("'run' is FALSE. Skip running tSNE.")
    return(obj)
  }
  pipe.name <- c(pipe.name, "RunTSNE")
  pipeMsg("Start")
  defaults <- list(
    reduction = obj@misc$reduction[['RunTSNE']] %||% "pca",
    reduction.key = "tSNE_",
    seed = 42,
    dim.embed = 2,
    num.threads = 1,
    perplexity = 30,
    reduction.name = "tsne"
  )
  params <- getDefaultArgs(defaults, params)
  params[['dims']] <- params[['dims']] %0%
    seq_len(ncol(Embeddings(obj, reduction = params[['reduction']])))
  params[['reduction.surfix']] <- params[['reduction.surfix']] %||%
    params[['reduction']]

  captureMsg(str(params))
  list2env(params, envir = environment())
  perplexity <- min(perplexity, floor((ncol(obj) - 1)/3))

  obj <- RunTSNE(
    object = obj,
    reduction = reduction,
    dims = dims,
    seed.use = seed,
    tsne.method = "Rtsne",
    dim.embed = dim.embed,
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    num_threads = num.threads,
    perplexity = perplexity,
    ...
  )
  obj[[paste0(reduction.name, "_", reduction.surfix)]] <- obj[[reduction.name]]
  obj
}

#' @importFrom Seurat RunUMAP
#' @importFrom SeuratObject Embeddings
#' @export
pipe_RunUMAP <- function(obj, params = list(), pipe.name = NULL, ...) {
  run <- params[['run']] %||% TRUE
  if (!run) {
    fastWarning("'run' is FALSE. Skip running UMAP.")
    return(obj)
  }

  pipe.name <- c(pipe.name, "RunUMAP")
  pipeMsg("Start")
  defaults <- list(
    reduction = obj@misc$reduction[['RunUMAP']] %||% "pca",
    assay = DefaultAssay(obj),
    reduction.key = "UMAP_",
    seed = 42,
    umap.method = "uwot",
    dim.embed = 2,
    n.neighbors = 30,
    metric = "cosine",
    learning.rate = 1,
    local.connectivity = 1,
    min.dist = 0.3,
    spread = 1,
    set.op.mix.ratio = 1,
    reduction.name = "umap"
  )
  params <- getDefaultArgs(defaults, params)
  params[['reduction.surfix']] <- params[['reduction.surfix']] %||%
    params[['reduction']]
  params[['dims']] <- params[['dims']] %0%
    seq_len(ncol(Embeddings(obj, reduction = params[['reduction']])))
  captureMsg(str(params))

  list2env(params, envir = environment())
  graph <- params[["graph"]]
  nn.name <- params[['nn.name']]
  reduction.model <- params[['reduction.model']]
  n.epochs <- params[['n.epochs']]
  if (!is.null(nn.name)) {
    dims <- NULL
  }

  obj <- RunUMAP(
    object = obj,
    dims = dims,
    reduction = reduction,
    graph = graph,
    assay = assay,
    nn.name = nn.name,
    umap.method = umap.method,
    reduction.model = reduction.model,
    return.model = TRUE,
    n.neighbors = as.integer(n.neighbors),
    n.components = as.integer(dim.embed),
    metric = metric,
    n.epochs = n.epochs,
    learning.rate = learning.rate,
    min.dist = min.dist,
    spread = spread,
    set.op.mix.ratio = set.op.mix.ratio,
    local.connectivity = local.connectivity,
    seed.use = as.integer(seed),
    verbose = TRUE,
    reduction.name = reduction.name,
    reduction.key = reduction.key,
    ...
  )
  obj[[paste0(reduction.name, "_", reduction.surfix)]] <- obj[[reduction.name]]
  obj
}

#' @importFrom Seurat FindClusters FindNeighbors
#' @export
pipe_FindClusters <- function(obj, params = list(), pipe.name = NULL, ...) {
  run <- params[['run']] %||% TRUE
  if (!run) {
    fastWarning("'run' is FALSE. Skip running FindClusters.")
    return(obj)
  }

  pipe.name <- c(pipe.name, "FindClusters")
  pipeMsg("Start")
  defaults <- list(
    reduction = obj@misc$reduction[['FindNeighbors']] %||% "pca",
    assay = DefaultAssay(obj),
    k.param = 20,
    prune.SNN = 1/15,
    n.trees = 50,
    l2.norm = FALSE,
    cluster.name = "seurat_clusters",
    resolution = 0.5,
    algorithm = "leiden",
    seed = 42,
    group.singletons = TRUE
  )
  params <- getDefaultArgs(defaults, params)
  params[['dims']] <- params[['dims']] %0%
    seq_len(ncol(Embeddings(obj, reduction = params[['reduction']])))
  captureMsg(str(params))

  list2env(params, envir = environment())
  graph.name <- params[['graph.name']]

  obj <- FindNeighbors(
    object = obj,
    assay = assay,
    reduction = reduction,
    dims = dims,
    k.param = k.param,
    prune.SNN = prune.SNN,
    n.trees = n.trees,
    verbose = FALSE,
    graph.name = graph.name,
    l2.norm = l2.norm,
    ...
  )
  if (seed < 1) {
    seed <- 1
  }

  obj <- FindClusters(
    object = obj,
    graph.name = graph.name,
    cluster.name = cluster.name,
    resolution = resolution,
    algorithm = algorithm,
    random.seed = as.integer(seed),
    group.singletons = group.singletons,
    verbose = TRUE,
    ...
  )

  pipeMsg("Check clustering results")
  captureMsg(table(obj[[cluster.name]], useNA = 'ifany'))

  obj@misc$colors[[cluster.name]] <- setColors(
    obj[[cluster.name, drop = TRUE]],
    type = "tsne",
    tag = "set1"
  )

  obj[['Clusters']] <- obj[[cluster.name]]
  obj@misc$colors[['Clusters']] <- obj@misc$colors[[cluster.name]]
  obj
}

#' @importFrom Seurat CCAIntegration IntegrateLayers RPCAIntegration
#' @importFrom harmony RunHarmony
#' @export
pipe_IntegrateData <- function(obj, params = list(), pipe.name = NULL, ...) {
  run <- params[['run']] %||% TRUE
  if (!run) {
    fastWarning("'run' is FALSE. Skip running IntegrateData.")
    return(obj)
  }

  pipe.name <- c(pipe.name, "IntegrateData")
  pipeMsg("Start")
  defaults <- list(
    method = "harmony",
    reduction = "pca",
    assay = DefaultAssay(obj),
    split.by = "orig.ident",
    theta = 2
  )
  params <- getDefaultArgs(defaults, params)
  params[['method']] <- match.arg(
    params[['method']],
    choices = c("harmony", "CCA", "RPCA")
  )
  default.new.reduc <- switch(
    params[['method']],
    harmony = "harmony",
    CCA = "cca.integrated",
    RPCA = "rpca.integrated"
  )
  params[['new.reduction']] <- params[['new.reduction']] %||% default.new.reduc

  captureMsg(str(params))
  list2env(params, envir = environment())

  if (length(unique(obj[[split.by, drop = TRUE]])) == 1) {
    fastWarning("Only 1 batch. No need for integration.")
    return(obj)
  }
  obj@misc$reduction <- list(
    RunUMAP = new.reduction,
    RunTSNE = new.reduction,
    FindNeighbors = new.reduction
  )
  if (method == "harmony") {
    obj <- RunHarmony(
      object = obj,
      group.by.vars = split.by,
      project.dim = TRUE,
      reduction.use = params[['reduction']],
      assay.use = assay,
      reduction.save = new.reduction,
      ...
    )
    return(obj)
  }
  method.f <- switch(
    method,
    CCA = CCAIntegration,
    RPCA = RPCAIntegration
  )
  features <- VariableFeatures(obj, assay = assay)
  obj[[assay]] <- split(
    obj[[assay]],
    f = obj[[split.by, drop = TRUE]],
    layers = c("counts", "data", "scale.data")
  )
  obj <- IntegrateLayers(
    obj,
    method = method.f,
    orig.reduction = reduction,
    new.reduction = new.reduction,
    verbose = TRUE,
    features = features,
    ...
  )
  obj <- JoinLayers(obj)
  obj
}

#' @export
pipe_StatFilterCells <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "StatFilterCells")
  pipeMsg("Start")

  defaults <- list(
    outdir = getwd(),
    group.by = c("Samples", "Groups")
  )
  params <- getDefaultArgs(defaults, params)
  captureMsg(str(params))
  list2env(params, envir = environment())

  group.by <- intersect(group.by, colnames(obj[[]]))

  if (length(group.by) == 0) {
    stop("No valid 'group.by'.")
  }
  if (!is.data.frame(obj@misc$colData)) {
    stop("No valid 'colData' in object@misc")
  }
  for (g in group.by) {
    raw <- obj@misc$colData
    filtered <- obj[[]]
    common.cols <- intersect(colnames(raw), colnames(filtered))
    df <- StatFilterCells(
      raw = raw[, common.cols, drop = FALSE],
      filtered = filtered[, common.cols, drop = FALSE],
      group.by = g
    )
    outfile <- paste(c("FilterCells.Stat", g, "xls"), collapse = ".")
    writeTable(df, file.path(outdir, outfile))
  }

  invisible(NULL)
}

#' @export
pipe_GetSeuratObj <- function(params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "GetSeuratObj")
  pipeMsg("Start")

  obj <- params[['obj']]
  if (!is(obj, "Seurat")) {
    obj_file <- normalizePath(params[['obj_file']], mustWork = TRUE)
    obj <- Load(obj_file)
  }

  obj <- pipe_RenameObject(obj, params, pipe.name = pipe.name)

  obj
}

#' @export
pipe_PostClustering <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "PostClustering")
  pipeMsg("Start")

  defaults <- list(
    outdir = getwd(),
    sample.by = "orig.ident",
    group.by = "Groups",
    cluster.by = "seurat_clusters",
    reductions = c("tsne", "umap"),
    corner.axis = TRUE
  )
  params <- getDefaultArgs(defaults, params)
  captureMsg(str(params))
  list2env(params, envir = environment())

  if (!group.by %in% colnames(obj[[]])) {
    obj[[group.by]] <- obj[[sample.by]]
    obj@misc$colors[[group.by]] <- obj@misc$colors[[sample.by]]
  }

  pipeMsg('Plotting dim_plot')
  save_dim_plots(
    obj,
    outdir = outdir,
    reductions = reductions,
    group.by = c(sample.by, cluster.by),
    split.by = sample.by,
    corner.axis = corner.axis,
    ...
  )
  save_dim_plots(
    obj,
    outdir = outdir,
    reductions = reductions,
    group.by = c(group.by, cluster.by),
    split.by = group.by,
    corner.axis = corner.axis,
    ...
  )

  pipeMsg('Get average expression')
  out.names <- c(Samples = sample.by, Groups = group.by, Clusters = cluster.by)
  for (i in names(out.names)) {
    outfile <- file.path(outdir, paste0("AllGene.", i, ".avg_exp.xls"))
    CalAvgExp(obj, group.by = out.names[i], outfile = outfile)

    outfile <- file.path(outdir, paste0("AllGene.", i, ".avg_pct.xls"))
    CalAvgPct(obj, group.by = out.names[i], outfile = outfile)
  }

  pipeMsg('Get clustering results')
  df <- FetchSeuratData(obj, out.names)
  writeTable(df, file.path(outdir, "Cells.cluster.list.xls"))

  invisible(NULL)
}

#' @export
pipe_FindAllMarkers <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "FindAllMarkers")
  pipeMsg("Start")

  defaults <- list(
    outdir = getwd(),
    assay = DefaultAssay(obj),
    group.by = "seurat_clusters",
    slot = 'data',
    test.use = "wilcox",
    p.adjust.method = "bonferroni",
    min.cells.group = 3,
    max.cells.per.ident = Inf,
    return.thresh = 1e-2,
    use.adjust = TRUE,
    logfc.threshold = 0.25,
    min.mean.exp = 0,
    min.pct = 0.01,
    min.diff.pct = -Inf,
    min.cells.feature = 3,
    only.pos = TRUE,
    min.exp = 0,
    pseudocount.use = 0.000001,
    base = 2,
    seed = 42
  )
  params <- getDefaultArgs(defaults, params)
  params[['features']] <- getArgList(params[['features']])
  params[['features']] <- getFeaturesID(obj, params[['features']])

  captureMsg(str(params))
  list2env(params, envir = environment())
  features <- params[['features']]
  latent.vars <- params[['latent.vars']]

  markers <- findAllMarkers(
    object = obj,
    assay = assay,
    features = features,
    group.by = group.by,
    slot = slot,

    test.use = test.use,
    p.adjust.method = p.adjust.method,
    latent.vars = latent.vars,
    min.cells.group = min.cells.group,
    max.cells.per.ident = max.cells.per.ident,
    return.thresh = return.thresh,
    use.adjust = use.adjust,

    logfc.threshold = logfc.threshold,
    min.mean.exp = min.mean.exp,
    min.pct = min.pct,
    min.diff.pct = min.diff.pct,
    min.cells.feature = min.cells.feature,
    only.pos = only.pos,

    min.exp = min.exp,
    pseudocount.use = pseudocount.use,
    base = base,
    seed = seed,
    ...
  )

  annot <- GetRowAnnot(obj, markers$gene)
  markers <- cbind(annot, as.data.frame(markers))
  markers <- markers %>%
    dplyr::mutate(Clusters = cluster, .before = 1) %>%
    select(!c("gene", "cluster")) %>%
    as.data.frame()
  save(markers, file = file.path(outdir, "markers.Rda"))

  pipeMsg("Output marker table: ", file.path(outdir, "DeGene.list.xls"))
  writeTable(markers, file.path(outdir, "DeGene.list.xls"))

  pipeMsg("Stat marker genes")
  StatMarker(markers, colors = obj@misc$colors[[group.by]])

  markers
}

#' @export
pipe_PostFindAllMarkers <- function(
    obj,
    markers,
    params = list(),
    pipe.name = NULL,
    ...
) {
  pipe.name <- c(pipe.name, "PostFindAllMarkers")
  pipeMsg("Start")

  defaults <- list(
    outdir = getwd(),
    top.n = 5,
    group.by = "seurat_clusters",
    coord.flip = TRUE,
    corner.axis = TRUE,
    reductions = "umap"
  )
  params <- getDefaultArgs(defaults, params)
  captureMsg(str(params))
  list2env(params, envir = environment())

  top <- getTopMarkers(markers, top.n = top.n, group.by = "Clusters")
  writeTable(top, file.path(outdir, "Top.list.xls"))

  features <- unique(top$GeneID)

  pipeMsg('Plotting dot_plot')
  save_dot_plots(
    obj,
    features = features,
    outdir = outdir,
    group.by = group.by,
    coord.flip = coord.flip
  )

  pipeMsg('Plotting violin_plot')
  mkdir(file.path(outdir, "violin_plots"))
  save_violin_plots(
    obj,
    features = features,
    outdir = file.path(outdir, "violin_plots"),
    combine = FALSE,
    group.by = group.by
  )

  pipeMsg('Plotting feature_dim_plot')
  mkdir(file.path(outdir, "feature_plots"))
  save_feature_dim_plots(
    obj,
    features = features,
    outdir = file.path(outdir, "feature_plots"),
    reductions = reductions,
    combine = FALSE,
    corner.axis = corner.axis
  )

  invisible(markers)
}

#' @export
pipe_MultiCore <- function(params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "MultiCore")
  pipeMsg("Start")

  n.cores <- params[['n.cores']] %||% 1

  require("future", quietly = TRUE)
  options(future.globals.maxSize = 100 * 1024 * 1024^2)
  if (n.cores > 1) {
    message("Use 'multicore', workers = ", n.cores)
    plan("multicore", workers = n.cores)
  }
  invisible(NULL)
}

#' @export
pipe_GroupCells <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "GroupCells")
  pipeMsg("Start")

  captureMsg(str(params))

  out0 <- data.frame(row.names = Cells(obj))
  if (length(params[['group.by']]) > 0) {
    group.by <- params[['group.by']]
    group.by <- Reduce(pasteFactors, as.list(obj[[group.by, drop = FALSE]]))
    cell.groups <- split(Cells(obj), f = group.by)
    params[['cell.groups']] <- c(params[['cell.groups']], cell.groups)
  }
  if (length(params[['cell.groups']]) > 0) {
    out <- group_cells(
      cell.groups = params[['cell.groups']],
      cells = rownames(out0)
    )
    out0 <- cbind(out0, out[rownames(out0), , drop = FALSE])
  }
  if (length(params[['meta.data']]) > 0) {
    out <- group_cells_by_meta(obj[[]], groups = params[['meta.data']])
    if (any(colnames(out) %in% colnames(out0))) {
      bad.names <- intersect(colnames(out), colnames(out0))
      stop("Duplicated group names:\n ", paste(bad.names, collapse = ", "))
    }
    out0 <- cbind(out0, out[rownames(out0), , drop = FALSE])
  }
  if (length(params[['features']]) > 0) {
    for (i in seq_along(params[['features']])) {
      genes <- names(params[['features']][[i]])
      genes <- getFeaturesID(obj, genes, uniq = TRUE)
      if (length(genes) == 0) {
        stop("No valid features found.")
      }
      if (length(genes) != length(params[['features']][[i]])) {
        stop("Invalid features found.")
      }
      if (any(!genes %in% rownames(obj))) {
        bad.genes <- setdiff(genes, rownames(obj))
        stop("Invalid features:\n ", paste(bad.genes, collapse = ", "))
      }
      names(params[['features']][[i]]) <- genes
    }
    out <- group_cells_by_features(
      GetAssayData(obj),
      groups = params[['features']]
    )
    if (any(colnames(out) %in% colnames(out0))) {
      bad.names <- intersect(colnames(out), colnames(out0))
      stop("Duplicated group names:\n ", paste(bad.names, collapse = ", "))
    }
    out0 <- cbind(out0, out[rownames(out0), , drop = FALSE])
  }
  n.cells <- colSums(out0)
  bad.groups <- names(n.cells)[n.cells == 0]
  if (length(bad.groups) > 0) {
    stop(
      "The following groups contain no cells:\n ",
      paste(bad.groups, collapse = ", ")
    )
  }


  out0
}

#' @export
pipe_GroupDiffer <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "GroupDiffer")
  pipeMsg("Start")

  outdir <- params[['outdir']] %||% getwd()

  defaults <- list(
    group.by = NULL,
    assay = NULL,
    reduction = NULL,
    features = NULL,
    slot = 'data',
    diff.type = "all",
    test.use = "MAST",
    p.adjust.method = "bonferroni",
    p.thresh = 0.05,
    use.adjust = TRUE,
    latent.vars = NULL,
    min.cells.group = 3,
    max.cells.per.ident = Inf,
    logfc.threshold = 0.1,
    min.mean.exp = 0,
    min.pct = 0.01,
    min.diff.pct = -Inf,
    min.cells.feature = 3,
    only.pos = FALSE,
    min.exp = 0,
    pseudocount.use = 0.000001,
    base = 2,
    seed = 42
  )
  params[['findGroupDiffer']] <- getDefaultArgs(
    defaults,
    params[['findGroupDiffer']]
  )
  captureMsg(str(params))

  differs <- params[['differs']]
  if (length(differs) == 0) {
    stop("No 'differs' found.")
  }
  group.data <- pipe_GroupCells(
    obj,
    params[['GroupCells']],
    pipe.name = pipe.name
  )

  params <- params[['findGroupDiffer']]
  params[['features']] <- getArgList(params[['features']])
  params[['features']] <- getFeaturesID(obj, params[['features']])
  list2env(params, envir = environment())

  markers <- findGroupDiffer(
    object = obj,
    group.data = group.data,
    differs = differs,
    assay = assay,
    group.by = group.by,
    features = features,
    slot = slot,
    diff.type = diff.type,
    test.use = test.use,
    p.adjust.method = p.adjust.method,
    p.thresh = p.thresh,
    use.adjust = use.adjust,
    filter.nosig = FALSE,
    latent.vars = latent.vars,
    min.cells.group = min.cells.group,
    max.cells.per.ident = max.cells.per.ident,
    logfc.threshold = logfc.threshold,
    min.mean.exp = min.mean.exp,
    min.pct = min.pct,
    min.diff.pct = min.diff.pct,
    min.cells.feature = min.cells.feature,
    only.pos = only.pos,
    mean.fxn = NULL,
    min.exp = min.exp,
    pseudocount.use = pseudocount.use,
    base = base,
    seed = seed,
    densify = FALSE,
    verbose = TRUE,
    ...
  )

  pipeMsg("Writing results for all genes")
  for (i in seq_along(markers)) {
    diff <- names(markers)[i]
    annot <- GetRowAnnot(obj, markers[[i]]$gene)
    markers[[i]] <- cbind(annot, as.data.frame(markers[[i]])) %>%
      dplyr::mutate(Clusters = cluster, .before = 1) %>%
      select(!c("gene", "cluster")) %>%
      as.data.frame()
    outfile <- paste0("GroupDiffer.", diff, ".all.xls")
    writeTable(markers[[i]], file.path(outdir, outfile))
  }
  save(markers, file = file.path(outdir, "markers.Rda"))

  pipeMsg("Writing results for significant genes")
  for (i in seq_along(markers)) {
    diff <- names(markers)[i]
    markers[[i]] <- markers[[i]] %>%
      dplyr::filter(significance != "nosig")
    outfile <- paste0("GroupDiffer.", diff, ".filtered.xls")
    writeTable(markers[[i]], file.path(outdir, outfile))
  }


  markers
}
