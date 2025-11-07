
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
  params$scale.data <- params$scale.data %||% FALSE
  params$misc.counts <- params$misc.counts %||% TRUE

  params$assays <- params$assays %||% Assays(obj)
  params$assays <- intersect(params$assays, Assays(obj))

  Message('>>>>> Shrinking Seurat object...')
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
  params$cells <- params$cells %||% params$cells_use
  params$features <- params$features %||% params$features

  Message('>>>>> Subset Seurat object...')
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
  params$Assay5 <- params$Assay5 %||% FALSE

  Message('>>>>> Check Seurat version...')
  obj <- UpdateSeuratAll(object = obj, Assay5 = params$Assay5)
  obj
}

#' @export
pipe_StatFeaturePercentage <- function(obj, params = list(), ...) {
  stat.pct <- params[['stat.pct']] %||% TRUE
  params[['stat.pct']] <- NULL

  Message('>>>>> Stat features')
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

  params[['assay']] <- params[['assay']] %||% DefaultAssay(obj)
  params[['vfeatures']] <- norm_list_param(params[['vfeatures']])
  if (length(params[['vfeatures']]) > 0) {
    Message("Use preset variable features")
    print(str(params))
    list2env(params, envir = environment())

    vfeatures <- getFeaturesID(obj, features = vfeatures)
    print(str(vfeatures))
    VariableFeatures(obj[[assay]]) <- vfeatures
    return(obj)
  }

  params[['selection.method']] <- params[['selection.method']] %||% "vst"
  params[['nfeatures']] <- params[['nfeatures']] %||% 2000
  params[['vfeature.must']] <- norm_list_param(params[['vfeature.must']])
  params[['vfeature.remove']] <- norm_list_param(params[['vfeature.remove']])
  print(str(params))
  list2env(params, envir = environment())

  Message("Find variable features")
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

  params[['s.features']] <- norm_list_param(params[['s.features']]) %0%
    Seurat::cc.genes$s.genes
  params[['g2m.features']] <- norm_list_param(params[['g2m.features']]) %0%
    Seurat::cc.genes$g2m.genes
  params[['s.features']] <- getFeaturesID(obj, params[['s.features']])
  params[['g2m.features']] <- getFeaturesID(obj, params[['g2m.features']])

  if (any(lengths(params[c('s.features', 'g2m.features')]) < 2)) {
    print(str(params))
    fastWarning("No enough cell cycle gene (n < 2). Skip cell cycle scoring.")
    return(obj)
  }

  params[['regress']] <- params[['regress']] %||% "diff"
  params[['assay']] <- params[['assay']] %||% DefaultAssay(obj)
  params[['slot']] <- params[['slot']] %0% "data"
  params[['set.ident']] <- params[['set.ident']] %0% FALSE
  print(str(params))
  list2env(params, envir = environment())
  ctrl <- params[['ctrl']]

  regress <- match.arg(params[['regress']], c("none", "diff", "all"))

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

  params[['method']] <- params[['method']] %||% "LogNormalize"
  params[['scale.factor']] <- params[['scale.factor']] %||% 10000
  params[['assay']] <- params[['assay']] %||% DefaultAssay(obj)
  print(str(params))

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

  params[['assay']] <- params[['assay']] %||% DefaultAssay(obj)
  params[['vars.to.regress']] <- params[['vars.to.regress']] %||%
    c(paste0("nCount_", params[['assay']]), "percent.mito")
  params[['vars.to.regress']] <- unique(c(
    params[['vars.to.regress']],
    obj@misc$vars.to.regress
  ))
  params[['scale.max']] <- params[['scale.max']] %||% 10
  params[['only.var.features']] <- params[['only.var.features']] %||% TRUE

  print(str(params))
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

  params[['assay']] <- params[['assay']] %||% DefaultAssay(obj)
  params[['npcs']] <- params[['npcs']] %||% 50
  params[['rev.pca']] <- params[['rev.pca']] %||% FALSE
  params[['weight.by.var']] <- params[['weight.by.var']] %||% TRUE
  params[['seed']] <- params[['seed']] %||% 42
  params[['approx']] <- params[['approx']] %||% TRUE

  params[['reduction.key']] <- params[['reduction.key']] %||% "PC_"
  params[['reduction.name']] <- params[['reduction.name']] %||% "pca"
  params[['reduction.surfix']] <- params[['reduction.surfix']] %||%
    params[['assay']]

  print(str(params))
  list2env(params, envir = environment())

  obj <- RunPCA(
    object = obj,
    assay = assay,
    npcs = params[['npcs']],
    rev.pca = params[['rev.pca']],
    weight.by.var = params[['weight.by.var']],
    verbose = FALSE,
    reduction.key = reduction.key,
    reduction.name = reduction.name,
    seed.use = seed,
    approx = params[['approx']],
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
  params[['reduction']] <- params[['reduction']] %||%
    obj@misc$reduction[['RunTSNE']]
  params[['reduction']] <- params[['reduction']] %||% "pca"

  params[['reduction.key']] <- params[['reduction.key']] %||% "tSNE_"
  params[['seed']] <- params[['seed']] %||% 42
  params[['dim.embed']] <- params[['dim.embed']] %||% 2
  params[['num.threads']] <- params[['num.threads']] %||% 1
  params[['dims']] <- params[['dims']] %0%
    seq_len(ncol(Embeddings(obj, reduction = params[['reduction']])))

  params[['perplexity']] <- params[['perplexity']] %||% 30
  params[['reduction.name']] <- params[['reduction.name']] %||% "tsne"
  params[['reduction.surfix']] <- params[['reduction.surfix']] %||%
    params[['reduction']]

  print(str(params))
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
  params[['reduction']] <- params[['reduction']] %||%
    obj@misc$reduction[['RunUMAP']]
  params[['reduction']] <- params[['reduction']] %||% "pca"

  params[['reduction.surfix']] <- params[['reduction.surfix']] %||%
    params[['reduction']]
  params[['reduction.name']] <- params[['reduction.name']] %||% "umap"
  params[['reduction.key']] <- params[['reduction.key']] %||% "UMAP_"
  params[['seed']] <- params[['seed']] %||% 42
  params[['dim.embed']] <- params[['dim.embed']] %||% 2
  params[['umap.method']] <- params[['umap.method']] %||% "uwot"
  params[['n.neighbors']] <- params[['n.neighbors']] %||% 30
  params[['metric']] <- params[['metric']] %||% "cosine"
  params[['learning.rate']] <- params[['learning.rate']] %||% 1
  params[['local.connectivity']] <- params[['local.connectivity']] %||% 1
  params[['min.dist']] <- params[['min.dist']] %||% 0.3
  params[["spread"]] <- params[["spread"]] %||% 1
  params[['set.op.mix.ratio']] <- params[['set.op.mix.ratio']] %||% 1
  params[['dims']] <- params[['dims']] %0%
    seq_len(ncol(Embeddings(obj, reduction = params[['reduction']])))

  print(str(params))
  list2env(params, envir = environment())

  graph <- params[["graph"]]
  assay <- params[['assay']]
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

  reduction <- params[['reduction']] %||%
    obj@misc$reduction[['FindNeighbors']]
  reduction <- reduction %||% "pca"

  assay <- params[['assay']] %||% DefaultAssay(obj)
  k.param <- params[['k.param']] %||% 20
  prune.SNN <- params[['prune.SNN']] %||% 1/15
  n.trees <- params[['n.trees']] %||% 50
  graph.name <- params[['graph.name']]
  l2.norm <- params[['l2.norm']] %||% FALSE

  dims <- params[['dims']] %0%
    seq_len(ncol(Embeddings(obj, reduction = reduction)))

  Message("Runing FindNeighbors.")
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

  cluster.name <- params[['cluster.name']] %||% "seurat_clusters"
  resolution <- params[['resolution']] %||% 0.5
  algorithm <- params[['algorithm']] %||% "leiden"
  seed <- params[['seed']] %||% 42
  group.singletons <- params[['group.singletons']] %||% TRUE
  if (seed < 1) {
    seed <- 1
  }

  Message("Runing FindClusters.")
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
  print(table(obj[[cluster.name]], useNA = 'ifany'))

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
  params[['split.by']] <- params[['split.by']] %||% "orig.ident"
  split.by <- params[['split.by']]
  if (length(unique(obj[[split.by, drop = TRUE]])) == 1) {
    fastWarning("Only 1 batch. No need for integration.")
    return(obj)
  }

  params[['method']] <- params[['method']] %||% "harmony"
  params[['reduction']] <- params[['reduction']] %||% "pca"
  params[['assay']] <- params[['assay']] %||% DefaultAssay(obj)
  assay <- params[['assay']]
  reduction <- params[['reduction']]

  method <- match.arg(params[['method']], choices = c("harmony", "CCA", "RPCA"))
  default.new.reduc <- switch(
    method,
    harmony = "harmony",
    CCA = "cca.integrated",
    RPCA = "rpca.integrated"
  )

  params[['new.reduction']] <- params[['new.reduction']] %||% default.new.reduc
  new.reduction <- params[['new.reduction']]

  obj@misc$reduction <- list(
    RunUMAP = new.reduction,
    RunTSNE = new.reduction,
    FindNeighbors = new.reduction
  )

  if (method == "harmony") {
    params[['theta']] <- params[['theta']] %||% 2
    print(str(params))
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
  print(str(params))
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
  mkdir(outdir)
  Message('>>>>> Save object to: ', outfile)
  save(obj, file = outfile)

  Message('>>>>> Output some lists: ', outdir)
  WriteRenameConf(obj, outdir)

  setwd(outdir)
  pipe_PostClustering(obj, params[['PostClustering']])

  obj
}

#' @export
pipe_PostClustering <- function(obj, params = list(), ...) {
  Message('>>>>> Running PostClustering')

  params[['outdir']] <- params[['outdir']] %||% getwd()
  params[['sample.by']] <- params[['sample.by']] %||% "orig.ident"
  params[['group.by']] <- params[['group.by']] %||% "Groups"
  params[['cluster.by']] <- params[['cluster.by']] %||% "seurat_clusters"
  params[['reductions']] <- params[['reductions']] %||% c("tsne", "umap")
  params[['corner.axis']] <- params[['corner.axis']] %||% TRUE
  print(str(params))
  list2env(params, envir = environment())

  if (!group.by %in% colnames(obj[[]])) {
    obj[[group.by]] <- obj[[sample.by]]
    obj@misc$colors[[group.by]] <- obj@misc$colors[[sample.by]]
  }

  Message('>>>>> Plotting DimPlot')
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

  Message('>>>>> Get average expression')
  out.names <- c(Samples = sample.by, Groups = group.by, Clusters = cluster.by)
  for (i in names(out.names)) {
    outfile <- paste0("AllGene.", i, ".avg_exp.xls")
    CalAvgExp(obj, group.by = out.names[i], outfile = outfile)

    outfile <- paste0("AllGene.", i, ".avg_pct.xls")
    CalAvgPct(obj, group.by = out.names[i], outfile = outfile)
  }

  Message('>>>>> Get clustering results')
  df <- FetchSeuratData(obj, out.names)
  writeTable(df, "Cells.cluster.list.xls")

  obj
}
