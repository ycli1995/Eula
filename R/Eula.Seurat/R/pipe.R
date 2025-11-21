#' @importFrom Eula.utils capture.msg fetch_default_params
NULL

#' @export
pipe_RenameObject <- function(obj, params = list()) {
  if (length(params[["RenameObject_yaml"]]) > 0) {
    file <- normalizePath(params[["RenameObject_yaml"]], mustWork = TRUE)
    Message("Read 'RenameObject_yaml': ", file)
    params <- readYAML(file2)
  }

  obj <- RenameSeuratColumns(obj, params[['rename_colnames']])

  obj <- RenameSeuratWrapper(obj, params[['rename']])

  obj <- pipe_UpdateSeurat(obj, params[['UpdateSeurat']])

  obj <- CheckMySeuratObj(obj, column.map = params$default_colnames)

  obj <- pipe_SubsetObject(obj, params[['subset']])

  obj <- pipe_ShrinkSeuratObject(obj, params)

  Message('>>>>> Finally check...')
  print(obj)
  print(str(obj@meta.data))
  print(obj@misc$colors)

  obj
}

#' @export
pipe_ShrinkSeuratObject <- function(obj, params = list(), ...) {
  Message('>>>>> Shrinking Seurat object...')

  defaults <- list(
    scale.data = FALSE,
    misc.counts = TRUE,
    assays = Assays(obj)
  )
  params <- fetch_default_params(defaults, params)
  params$assays <- intersect(params$assays, Assays(obj))
  params <- params[names(defaults)]

  capture.msg(str(params))
  list2env(params, envir = environment())

  obj <- ShrinkSeuratObject(
    object = obj,
    assays = params$assays,
    scale.data = params$scale.data,
    misc.counts = params$misc.counts
  )
  obj
}

#' @export
pipe_SubsetObject <- function(obj, params = list(), ...) {
  Message('>>>>> Subset Seurat object...')

  params$cells <- params$cells %||% params$cells_use
  params$features <- params$features %||% params$features_use
  capture.msg(str(params))

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
pipe_UpdateSeurat <- function(obj, params = list(), ...) {
  Message('>>>>> Check Seurat version...')
  params$Assay5 <- params$Assay5 %||% FALSE
  capture.msg(str(params))

  obj <- UpdateSeuratAll(object = obj, Assay5 = params$Assay5)
  obj
}

#' @export
pipe_StatFeaturePercentage <- function(obj, params = list(), ...) {
  Message('>>>>> Stat features')
  stat.pct <- params[['stat.pct']] %||% TRUE
  capture.msg(str(params))

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
pipe_MakeSeuratObj <- function(params) {
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
    stop("No 'name_list' for MakeSeuratObj")
  }
  assay <- params[['assay']] %||% "RNA"

  Message('>>>>> Making Seurat object')
  obj <- MakeSeuratObj(
    data.paths = data.paths,
    data.names = data.names,
    assay = assay,
    name.list = name.list
  )
  obj
}

#' @importFrom Seurat FindVariableFeatures
#' @importFrom SeuratObject VariableFeatures VariableFeatures<-
#' @export
pipe_FindVariableFeatures <- function(obj, params = list(), ...) {
  Message("Runing FindVariableFeatures")

  defaults <- list(
    assay = DefaultAssay(obj),
    selection.method = "vst",
    nfeatures = 2000
  )
  params <- fetch_default_params(defaults, params)
  params[['vfeatures']] <- norm_list_param(params[['vfeatures']])
  params[['vfeature.must']] <- norm_list_param(params[['vfeature.must']])
  params[['vfeature.remove']] <- norm_list_param(params[['vfeature.remove']])
  capture.msg(str(params))

  if (length(params[['vfeatures']]) > 0) {
    Message("Use preset variable features")
    list2env(params, envir = environment())

    vfeatures <- getFeaturesID(obj, features = vfeatures)
    VariableFeatures(obj[[assay]]) <- vfeatures
    return(obj)
  }

  Message("Find variable features")
  list2env(params, envir = environment())
  obj <- FindVariableFeatures(
    object = obj,
    assay = assay,
    selection.method = selection.method,
    nfeatures = nfeatures,
    ...
  )

  if (length(params[['vfeature.must']]) > 0) {
    Message("Check vfeatures required")
    vfeature.must <- getFeaturesID(obj, features = vfeature.must)
    VariableFeatures(obj) <- union(VariableFeatures(obj), vfeature.must)
  }

  if (length(params[['vfeature.remove']]) > 0) {
    Message("Check vfeatures excluded")
    vfeature.remove <- getFeaturesID(obj, features = vfeature.remove)
    VariableFeatures(obj) <- setdiff(VariableFeatures(obj), vfeature.remove)
  }
  obj
}

#' @export
pipe_CellCycleScoring <- function(obj, params = list(), ...) {
  Message("Runing cell cycle scoring.")

  defaults <- list(
    s.features = Seurat::cc.genes$s.genes,
    g2m.features = Seurat::cc.genes$g2m.genes,
    regress = "diff",
    assay = DefaultAssay(obj),
    slot = "data",
    set.ident = FALSE
  )

  params[['s.features']] <- norm_list_param(params[['s.features']])
  params[['g2m.features']] <- norm_list_param(params[['g2m.features']])
  params <- fetch_default_params(defaults, params)
  capture.msg(str(params))

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
  obj@misc$vars.to.regress <- c("S.Score", "G2M.Score", obj@misc$vars.to.regress)
  obj
}

#' @importFrom Seurat NormalizeData
#' @export
pipe_NormalizeData <- function(obj, params = list(), ...) {
  Message("Runing NormalizeData")

  defaults <- list(
    method = "LogNormalize",
    scale.factor = 10000,
    assay = DefaultAssay(obj)
  )
  params <- fetch_default_params(defaults, params)
  capture.msg(str(params))

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
pipe_ScaleData <- function(obj, params = list(), ...) {
  Message("Runing ScaleData.")

  defaults <- list(
    assay = DefaultAssay(obj),
    scale.max = 10,
    only.var.features = TRUE
  )
  params <- fetch_default_params(defaults, params)
  params[['vars.to.regress']] <- params[['vars.to.regress']] %||%
    c(paste0("nCount_", params[['assay']]), "percent.mito")
  params[['vars.to.regress']] <- unique(c(
    params[['vars.to.regress']],
    obj@misc$vars.to.regress
  ))
  capture.msg(str(params))
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
pipe_RunPCA <- function(obj, params = list(), ...) {
  Message("Runing PCA.")

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
  params <- fetch_default_params(defaults, params)
  params[['reduction.surfix']] <- params[['reduction.surfix']] %||%
    params[['assay']]
  capture.msg(str(params))
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
pipe_RunTSNE <- function(obj, params = list(), ...) {
  run <- params[['run']] %||% TRUE
  if (!run) {
    fastWarning("'run' is FALSE. Skip running tSNE.")
    return(obj)
  }

  Message("Runing tSNE.")
  defaults <- list(
    reduction = obj@misc$reduction[['RunTSNE']] %||% "pca",
    reduction.key = "tSNE_",
    seed = 42,
    dim.embed = 2,
    num.threads = 1,
    perplexity = 30,
    reduction.name = "tsne"
  )
  params <- fetch_default_params(defaults, params)
  params[['dims']] <- params[['dims']] %0%
    seq_len(ncol(Embeddings(obj, reduction = params[['reduction']])))
  params[['reduction.surfix']] <- params[['reduction.surfix']] %||%
    params[['reduction']]

  capture.msg(str(params))
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
pipe_RunUMAP <- function(obj, params = list(), ...) {
  run <- params[['run']] %||% TRUE
  if (!run) {
    fastWarning("'run' is FALSE. Skip running UMAP.")
    return(obj)
  }

  Message("Runing UMAP.")
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
  params <- fetch_default_params(defaults, params)
  params[['reduction.surfix']] <- params[['reduction.surfix']] %||%
    params[['reduction']]
  params[['dims']] <- params[['dims']] %0%
    seq_len(ncol(Embeddings(obj, reduction = params[['reduction']])))
  capture.msg(str(params))

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
pipe_FindClusters <- function(obj, params = list(), ...) {
  run <- params[['run']] %||% TRUE
  if (!run) {
    fastWarning("'run' is FALSE. Skip running FindClusters.")
    return(obj)
  }

  Message("Runing FindClusters.")
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
  params <- fetch_default_params(defaults, params)
  params[['dims']] <- params[['dims']] %0%
    seq_len(ncol(Embeddings(obj, reduction = params[['reduction']])))
  capture.msg(str(params))

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
  capture.msg(table(obj[[cluster.name]], useNA = 'ifany'))

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
pipe_IntegrateData <- function(obj, params = list(), ...) {
  run <- params[['run']] %||% TRUE
  if (!run) {
    fastWarning("'run' is FALSE. Skip running IntegrateData.")
    return(obj)
  }

  Message('>>>>> Running IntegrateData')
  defaults <- list(
    method = "harmony",
    reduction = "pca",
    assay = DefaultAssay(obj),
    split.by = "orig.ident",
    theta = 2
  )
  params <- fetch_default_params(defaults, params)
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

  capture.msg(str(params))
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
pipe_StatFilterCells <- function(obj, params = list(), ...) {
  outdir <- params[['outdir']] %||% getwd()

  group.by <- params[['group.by']] %||% "Samples"
  group.by <- intersect(group.by, colnames(obj[[]]))

  Message('>>>>> Stat filtered cells')
  if (length(group.by) == 0) {
    stop("No valid 'group.by'.")
  }
  if (!is.data.frame(obj@misc$colData)) {
    stop("No valid 'colData' in object@misc")
  }
  raw <- obj@misc$colData
  filtered <- obj[[]]
  common.cols <- intersect(colnames(raw), colnames(filtered))
  df <- StatFilterCells(
    raw = raw[, common.cols, drop = FALSE],
    filtered = filtered[, common.cols, drop = FALSE],
    group.by = group.by
  )
  outfile <- paste(c("FilterCells.Stat", group.by, "xls"), collapse = ".")
  writeTable(df, file.path(outdir, outfile))
}

#' @export
pipe_GetSeuratObj <- function(params = list(), ...) {
  obj <- params[['obj']]
  if (!is(obj, "Seurat")) {
    obj_file <- normalizePath(params[['obj_file']], mustWork = TRUE)
    obj <- Load(obj_file)
  }

  obj <- pipe_RenameObject(obj, params)
  obj
}

#' @export
pipe_MakeSeuratObj_R <- function(params) {
  outfile <- params[['obj_out']] %||% file.path(getwd(), "obj.Rda")
  outdir <- dirname(outfile)

  checkKeys(params, "MakeSeuratObj")

  obj <- pipe_MakeSeuratObj(params[["MakeSeuratObj"]])

  obj <- pipe_StatFeaturePercentage(obj, params[["StatFeaturePercentage"]])

  obj <- RenameSeuratColumns(obj, params[['rename_colnames']])

  obj <- RenameSeuratWrapper(obj, params[['rename']])

  obj <- pipe_UpdateSeurat(obj, params[['UpdataSeurat']])

  obj <- CheckMySeuratObj(obj, column.map = params[['default_colnames']])

  obj@misc$colData <- obj[[]]

  obj <- pipe_SubsetObject(obj, params[['subset']])

  mkdir(outdir)
  params[['StatFilterCells']][['outdir']] <- outdir
  pipe_StatFilterCells(obj, params[['StatFilterCells']])

  Message('>>>>> Save object to: ', outfile)
  save(obj, file = outfile)

  Message('>>>>> Output some lists: ', outdir)
  WriteRenameConf(obj, outdir)

  obj
}

#' @export
pipe_RenameObject_R <- function(params = list(), ...) {
  outfile <- params[['obj_out']]
  if (length(outfile) == 0) {
    stop("Missing 'obj_out' params.")
  }
  if (is.na(outfile)) {
    stop("Invalid 'obj_out' params.")
  }

  obj <- pipe_GetSeuratObj(params, ...)

  outdir <- dirname(outfile)
  mkdir(outdir)
  Message('>>>>> Save object to: ', outfile)
  save(obj, file = outfile)

  Message('>>>>> Output some lists: ', outdir)
  WriteRenameConf(obj, outdir)

  obj
}

#' @export
pipe_SeuratClustering_R <- function(params = list(), ...) {
  outdir <- params[['outdir']] %||% getwd()
  pipe_MultiCore(params)

  obj <- pipe_GetSeuratObj(params)

  obj <- pipe_NormalizeData(obj, params[['NormalizeData']])

  obj <- pipe_FindVariableFeatures(obj, params[['FindVariableFeatures']])

  obj <- pipe_CellCycleScoring(obj, params[['CellCycleScoring']])

  obj <- pipe_ScaleData(obj, params[['ScaleData']])

  obj <- pipe_RunPCA(obj, params[['RunPCA']])

  obj <- pipe_IntegrateData(obj, params[['IntegrateData']])

  obj <- pipe_RunTSNE(obj, params[['RunTSNE']])

  obj <- pipe_RunUMAP(obj, params[['RunUMAP']])

  obj <- pipe_FindClusters(obj, params[['FindClusters']])

  outfile <- file.path(outdir, "obj.Rda")
  mkdir(outdir, chdir = TRUE)
  Message('>>>>> Save object to: ', outfile)
  save(obj, file = outfile)

  Message('>>>>> Output some lists: ', outdir)
  WriteRenameConf(obj)

  pipe_PostClustering(obj, params[['PostClustering']])

  invisible(obj)
}

#' @export
pipe_PostClustering <- function(obj, params = list(), ...) {
  Message('>>>>> Running PostClustering')
  defaults <- list(
    outdir = getwd(),
    sample.by = "orig.ident",
    group.by = "Groups",
    cluster.by = "seurat_clusters",
    reductions = c("tsne", "umap"),
    corner.axis = TRUE
  )
  params <- fetch_default_params(defaults, params)
  capture.msg(str(params))
  list2env(params, envir = environment())

  if (!group.by %in% colnames(obj[[]])) {
    obj[[group.by]] <- obj[[sample.by]]
    obj@misc$colors[[group.by]] <- obj@misc$colors[[sample.by]]
  }

  Message('Plotting DimPlot')
  save_multi_DimPlot(
    obj,
    outdir = outdir,
    reductions = reductions,
    group.by = c(sample.by, cluster.by),
    split.by = sample.by,
    corner.axis = corner.axis,
    ...
  )
  save_multi_DimPlot(
    obj,
    outdir = outdir,
    reductions = reductions,
    group.by = c(group.by, cluster.by),
    split.by = group.by,
    corner.axis = corner.axis,
    ...
  )

  Message('Get average expression')
  out.names <- c(Samples = sample.by, Groups = group.by, Clusters = cluster.by)
  for (i in names(out.names)) {
    outfile <- file.path(outdir, paste0("AllGene.", i, ".avg_exp.xls"))
    CalAvgExp(obj, group.by = out.names[i], outfile = outfile)

    outfile <- file.path(outdir, paste0("AllGene.", i, ".avg_pct.xls"))
    CalAvgPct(obj, group.by = out.names[i], outfile = outfile)
  }

  Message('Get clustering results')
  df <- FetchSeuratData(obj, out.names)
  writeTable(df, file.path(outdir, "Cells.cluster.list.xls"))

  obj
}

#' @export
pipe_FindAllMarkers <- function(obj, params = list(), ...) {
  Message('>>>>> Running FindAllMarkers')

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
  params <- fetch_default_params(defaults, params)
  params[['features']] <- norm_list_param(params[['features']])
  params[['features']] <- getFeaturesID(obj, params[['features']])

  capture.msg(str(params))
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
  writeTable(markers, file.path(outdir, "DeGene.list.xls"))

  StatMarker(markers, colors = obj@misc$colors[[group.by]])

  markers
}

#' @export
pipe_PostFindAllMarkers <- function(obj, markers, params = list(), ...) {
  Message('>>>>> Running post-FindAllMarkers')
  defaults <- list(
    outdir = getwd(),
    top.n = 5,
    group.by = "seurat_clusters",
    coord.flip = TRUE
  )
  params <- fetch_default_params(defaults, params)
  capture.msg(str(params))
  list2env(params, envir = environment())

  top <- getTopMarkers(markers, top.n = top.n, group.by = "Clusters")
  writeTable(top, file.path(outdir, "Top.list.xls"))

  features <- unique(top$GeneID)

  Message('Plotting DotPlot')
  save_multi_DotPlot(
    obj,
    features = features,
    outdir = outdir,
    group.by = group.by,
    coord.flip = coord.flip,
    ...
  )

  invisible(markers)
}

#' @export
pipe_MultiCore <- function(params = list(), ...) {
  n.cores <- params[['n.cores']] %||% 1

  require("future", quietly = TRUE)
  options(future.globals.maxSize = 100 * 1024 * 1024^2)
  if (n.cores > 1) {
    Message("Set 'multicore': ", n.cores)
    plan("multicore", workers = n.cores)
  }
  invisible(NULL)
}

#' @export
pipe_GroupCells <- function(obj, params = list(), ...) {
  Message("Grouping cells")
  capture.msg(str(params))

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
pipe_GroupDiffer <- function(obj, params = list(), ...) {
  Message("Running GroupDiffer")

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
  params[['findGroupDiffer']] <- fetch_default_params(
    defaults = defaults,
    params = params[['findGroupDiffer']]
  )
  capture.msg(str(params))

  differs <- params[['differs']]
  if (length(differs) == 0) {
    stop("No 'differs' found.")
  }
  group.data <- pipe_GroupCells(obj, params[['GroupCells']])

  params <- params[['findGroupDiffer']]
  params[['features']] <- norm_list_param(params[['features']])
  params[['features']] <- getFeaturesID(obj, params[['features']])

  list2env(params, envir = environment())
  features <- params[['features']]
  assay <- params[['assay']]
  reduction <- params[['reduction']]
  group.by <- params[['group.by']]
  latent.vars <- params[['latent.vars']]

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

  for (i in seq_along(markers)) {
    diff <- names(markers)[i]
    markers[[i]] <- markers[[i]] %>%
      dplyr::filter(significance != "nosig")
    outfile <- paste0("GroupDiffer.", diff, ".filtered.xls")
    writeTable(markers[[i]], file.path(outdir, outfile))
  }

  markers
}

#' @export
pipe_FindAllMarkers_R <- function(params = list(), ...) {
  outdir <- params[['outdir']] %||% getwd()

  pipe_MultiCore(params)

  obj <- pipe_GetSeuratObj(params)

  mkdir(outdir, chdir = TRUE)
  markers <- pipe_FindAllMarkers(obj, params[['FindAllMarkers']])

  markers <- pipe_PostFindAllMarkers(
    obj = obj,
    markers = markers,
    params = params[['PostFindAllMarkers']]
  )

  invisible(markers)
}

#' @export
pipe_GroupDiffer_R <- function(params = list(), ...) {
  outdir <- params[['outdir']] %||% getwd()

  pipe_MultiCore(params)

  obj <- pipe_GetSeuratObj(params)

  mkdir(outdir, chdir = TRUE)
  params[['GroupDiffer']][['outdir']] <- outdir
  markers <- pipe_GroupDiffer(obj, params[['GroupDiffer']])

  invisible(markers)
}
