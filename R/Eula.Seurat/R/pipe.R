
#' @export
pipe_RenameObject <- function(obj, parameter = list()) {

  print(str(parameter))

  obj <- RenameSeuratColumns(obj, parameter[['rename_colnames']])

  obj <- RenameSeuratWrapper(obj, parameter[['rename']])

  obj <- pipe_UpdateSeurat(obj, parameter[['UpdateSeurat']])

  obj <- CheckMySeuratObj(obj, column.map = parameter$default_colnames)

  obj <- pipe_SubsetObject(obj, parameter[['subset']])

  obj <- pipe_ShrinkSeuratObject(obj, parameter)

  Message('>>>>> Finally check...')
  print(obj)
  print(str(obj@meta.data))
  print(obj@misc$colors)

  obj
}

#' @export
pipe_ShrinkSeuratObject <- function(obj, parameter = list(), ...) {
  parameter$scale.data <- parameter$scale.data %||% FALSE
  parameter$misc.counts <- parameter$misc.counts %||% TRUE

  parameter$assays <- parameter$assays %||% Assays(obj)
  parameter$assays <- intersect(parameter$assays, Assays(obj))

  Message('>>>>> Shrinking Seurat object...')
  obj <- ShrinkSeuratObject(
    object = obj,
    assays = parameter$assays,
    scale.data = parameter$scale.data,
    misc.counts = parameter$misc.counts
  )
  obj
}

#' @export
pipe_SubsetObject <- function(obj, parameter = list(), ...) {
  parameter$cells <- parameter$cells %||% parameter$cells_use
  parameter$features <- parameter$features %||% parameter$features

  Message('>>>>> Subset Seurat object...')
  cells_features <- get_cells_and_features(
    object = obj,
    cells = parameter$cells,
    features = parameter$features,
    cells_exclude = parameter$cells_exclude,
    features_exclude = parameter$features_exclude
  )

  parameter <- parameter[intersect(names(parameter), colnames(obj[[]]))]
  obj <- SubsetObject(
    object = obj,
    column_map = parameter,
    cells = cells_features[['cells']],
    features = cells_features[['features']]
  )
  obj
}

#' @export
pipe_UpdateSeurat <- function(obj, parameter = list(), ...) {
  parameter$Assay5 <- parameter$Assay5 %||% FALSE

  Message('>>>>> Check Seurat version...')
  obj <- UpdateSeuratAll(object = obj, Assay5 = parameter$Assay5)
  obj
}

#' @export
pipe_StatFeaturePercentage <- function(obj, parameter = list(), ...) {
  stat.pct <- parameter[['stat.pct']] %||% TRUE
  parameter[['stat.pct']] <- NULL

  Message('>>>>> Stat features')
  for (i in names(parameter)) {
    i2 <- gsub("\\/| ", "_", i)
    obj <- StatFeatures(
      object = obj,
      features = parameter[[i]],
      col.name = paste0("percent.", i2),
      stat.pct = stat.pct
    )
  }
  obj
}

#' @export
pipe_MakeSeuratObj <- function(parameter) {
  if (length(parameter) == 0) {
    stop("'MakeSeuratObj' parameter is empty.")
  }
  data.paths <- parameter[["paths"]]
  if (length(data.paths) == 0) {
    stop("No 'paths' for matrices.")
  }
  data.names <- parameter[["names"]]
  if (length(data.names) != length(data.paths)) {
    stop("'names' must has the same length as 'paths'")
  }
  name.list <- parameter[['name_list']]
  if (length(name.list) == 0) {
    stop("No 'name_list' for MakeSeuratObj")
  }
  assay <- parameter[['assay']] %||% "RNA"

  Message('>>>>> Making Seurat object')
  obj <- MakeSeuratObj(
    data.paths = data.paths,
    data.names = data.names,
    assay = assay,
    name.list = name.list
  )
  obj
}

pipe_FindVariableFeatures <- function(obj, parameter = list(), ...) {
  assay <- parameter[['assay']] %||% DefaultAssay(obj)
  selection.method <- parameter[['selection.method']] %||% "vst"
  nfeatures <- parameter[['nfeatures']] %||% 2000

  vfeatures <- norm_list_param(parameter[['vfeatures']])
  if (length(vfeatures) > 0) {
    Message("Use preset variable features")
    vfeatures <- getFeaturesID(obj, features = vfeatures)
    print(str(vfeatures))
    VariableFeatures(obj[[assay]]) <- vfeatures
    return(obj)
  }

  Message("Find variable features")
  obj <- FindVariableFeatures(
    object = obj,
    assay = assay,
    selection.method = selection.method,
    nfeatures = nfeatures,
    ...
  )
  vfeature.must <- norm_list_param(parameter[['vfeature.must']])
  if (length(vfeature.must) > 0) {
    Message("Check vfeatures required")
    vfeature.must <- getFeaturesID(obj, features = vfeature.must)
    VariableFeatures(obj) <- union(VariableFeatures(obj), vfeature.must)
  }

  vfeature.remove <- norm_list_param(parameter[['vfeature.remove']])
  if (length(vfeature.remove) > 0) {
    Message("Check vfeatures excluded")
    vfeature.remove <- getFeaturesID(obj, features = vfeature.remove)
    VariableFeatures(obj) <- setdiff(VariableFeatures(obj), vfeature.remove)
  }
  obj
}

#' @export
pipe_CellCycleScoring <- function(obj, parameter = list(), ...) {
  Message("Runing cell cycle scoring.")

  s.features <- norm_list_param(parameter[['s.features']]) %0%
    Seurat::cc.genes$s.genes
  g2m.features <- norm_list_param(parameter[['s.features']]) %0%
    Seurat::cc.genes$g2m.genes
  set.ident <- parameter[['set.ident']] %0% FALSE
  assay <- parameter[['assay']]
  slot <- parameter[['slot']] %0% "data"
  ctrl <- parameter[['ctrl']]

  s.features <- getFeaturesID(obj, s.features)
  g2m.features <- getFeaturesID(obj, g2m.features)

  if (length(s.features) < 2 || length(g2m.features) < 2) {
    fastWarning("No enough cell cycle gene (n < 2). Skip cell cycle scoring.")
    return(obj)
  }
  obj <- CellCycle(
    object = obj,
    s.features = s.features,
    g2m.features = g2m.features,
    set.ident = set.ident,
    set.ident = set.ident,
    slot = slot,
    assay = assay,
    ctrl = ctrl,
    ...
  )
  obj
}

#' @export
pipe_NormalizeData <- function(obj, parameter = list(), ...) {
  method <- parameter[['method']] %||% "LogNormalize"
  scale.factor <- parameter[['scale.factor']] %||% 10000
  assay <- parameter[['assay']] %||% DefaultAssay(obj)

  Message("Runing NormalizeData.")
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

#' @export
pipe_ScaleData <- function(obj, parameter = list(), ...) {
  vars.to.regress <- parameter[['vars.to.regress']]
  assay <- parameter[['assay']] %||% DefaultAssay(obj)
  scale.max <- parameter[['scale.max']] %||% 10
  only.var.features <- parameter[['only.var.features']] %||% TRUE

  features <- if (only.var.features) NULL else rownames(object)

  Message("Runing ScaleData.")
  obj <- ScaleData(
    object = obj,
    assay = assay,
    vars.to.regress = vars.to.regress,
    features = features,
    scale.max = scale.max,
    verbose = FALSE,
    ...
  )
  obj
}

#' @importFrom Seurat RunPCA
#' @export
pipe_RunPCA <- function(obj, parameter = list(), ...) {
  assay <- parameter[['assay']] %||% DefaultAssay(obj)
  reduction.name <- parameter[['reduction.name']] %||% "pca"
  npcs <- parameter[['npcs']] %||% 50
  rev.pca <- parameter[['rev.pca']] %||% FALSE
  weight.by.var <- parameter[['weight.by.var']] %||% TRUE
  reduction.key <- parameter[['reduction.key']] %||% "PC_"
  seed <- parameter[['seed']] %||% 42
  approx <- parameter[['approx']] %||% TRUE
  reduction.surfix <- parameter[['reduction.surfix']] %||% assay

  Message("Runing PCA.")
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
  obj
}

#' @importFrom SeuratObject Embeddings
#' @importFrom Seurat RunTSNE
#' @export
pipe_RunTSNE <- function(obj, parameter = list(), ...) {

  reduction <- parameter[['reduction']] %||% obj@misc$reduction[['RunTSNE']]
  reduction <- reduction %||% "pca"

  reduction.surfix <- parameter[['reduction.surfix']] %||% reduction
  perplexity <- parameter[['perplexity']] %||% 30
  reduction.name <- parameter[['reduction.name']] %||% "tsne"
  reduction.key <- parameter[['reduction.key']] %||% "tSNE_"
  seed <- seed %||% parameter[['seed']] %||% 42
  dim.embed <- dim.embed %||% parameter[['dim.embed']] %||% 2
  num.threads <- parameter[['num.threads']] %||% 1

  dims <- parameter[['dims']] %0%
    seq_len(ncol(Embeddings(obj, reduction = reduction)))

  perplexity <- min(perplexity, floor((ncol(obj) - 1)/3))

  Message("Runing tSNE.")
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
pipe_RunUMAP <- function(obj, parameter = list(), ...) {
  reduction <- parameter[['reduction']] %||% obj@misc$reduction[['RunUMAP']]
  reduction <- reduction %||% "pca"

  reduction.surfix <- parameter[['reduction.surfix']] %||% reduction
  reduction.name <- parameter[['reduction.name']] %||% "umap"
  reduction.key <- parameter[['reduction.key']] %||% "UMAP_"
  seed <- seed %||% parameter[['seed']] %||% 42
  dim.embed <- dim.embed %||% parameter[['dim.embed']] %||% 2

  graph <- parameter[["graph"]]
  assay <- parameter[['assay']]
  nn.name <- parameter[['nn.name']]
  umap.method <- parameter[['umap.method']] %||% "uwot"
  reduction.model <- parameter[['reduction.model']]
  n.neighbors <- parameter[['n.neighbors']] %||% 30
  metric <- parameter[['metric']] %||% "cosine"
  n.epochs <- parameter[['n.epochs']]
  learning.rate <- parameter[['learning.rate']] %||% 1
  min.dist <- parameter[['min.dist']] %||% 0.3
  spread <- parameter[["spread"]] %||% 1
  set.op.mix.ratio <- parameter[['set.op.mix.ratio']] %||% 1
  local.connectivity <- parameter[['local.connectivity']]

  dims <- NULL
  if (is.null(nn.name)) {
    dims <- parameter[['dims']] %0%
      seq_len(ncol(Embeddings(obj, reduction = reduction)))
  }

  Message("Runing UMAP.")
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
pipe_FindClusters <- function(obj, parameter = list(), ...) {
  reduction <- parameter[['reduction']] %||%
    obj@misc$reduction[['FindNeighbors']]
  reduction <- reduction %||% "pca"

  assay <- parameter[['assay']] %||% DefaultAssay(obj)
  k.param <- parameter[['k.param']] %||% 20
  prune.SNN <- parameter[['prune.SNN']] %||% 1/15
  n.trees <- parameter[['n.trees']] %||% 50
  graph.name <- parameter[['graph.name']]
  l2.norm <- parameter[['l2.norm']] %||% FALSE

  dims <- parameter[['dims']] %0%
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

  cluster.name <- parameter[['cluster.name']] %||% "seurat_clusters"
  resolution <- parameter[['resolution']] %||% 0.5
  algorithm <- parameter[['algorithm']] %||% "leiden"
  seed <- parameter[['seed']] %||% 42
  group.singletons <- parameter[['group.singletons']] %||% TRUE
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
  obj
}

#' @importFrom Seurat CCAIntegration IntegrateLayers RPCAIntegration
#' @importFrom harmony RunHarmony
#' @export
pipe_IntegrateData <- function(obj, parameter = list(), ...) {
  method <- parameter[['method']] %||% "harmony"
  split.by <- parameter[['split.by']] %||% "orig.ident"
  reduction <- parameter[['reduction']] %||% "pca"
  assay <- parameter[['assay']] %||% DefaultAssay(obj)

  method <- match.arg(method, choices = c("harmony", "CCA", "RPCA"))

  default.new.reduc <- switch(
    method,
    harmony = "harmony",
    CCA = "cca_integrated",
    RPCA = "rpca_integrated"
  )
  new.reduction <- parameter[['new.reduction']] %||% default.new.reduc
  obj@misc$reduction <- list(
    RunUMAP = new.reduction,
    RunTSNE = new.reduction,
    FindNeighbors = new.reduction
  )

  if (method == "harmony") {
    theta <- parameter[['theta']] %||% 2
    obj <- RunHarmony(
      object = obj,
      group.by.vars = split.by,
      project.dim = TRUE,
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
pipe_StatFilterCells <- function(obj, parameter = list(), ...) {
  group.by <- parameter[['group.by']] %||% "Samples"
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
  group.by <- paste(group.by, collapse = ".")
  writeTable(df, paste0("FilterCells.Stat.", group.by, ".xls"))
}

#' @export
pipe_MakeSeuratObj_R <- function(parameter) {
  checkKeys(parameter, "MakeSeuratObj")

  obj <- pipe_MakeSeuratObj(parameter[["MakeSeuratObj"]])

  obj <- pipe_StatFeaturePercentage(obj, parameter[["StatFeaturePercentage"]])

  obj <- RenameSeuratColumns(obj, parameter[['rename_colnames']])

  obj <- RenameSeuratWrapper(obj, parameter[['rename']])

  obj <- pipe_UpdateSeurat(obj, parameter[['UpdataSeurat']])

  obj <- CheckMySeuratObj(obj, column.map = parameter[['default_colnames']])

  obj@misc$colData <- obj[[]]

  obj <- pipe_SubsetObject(obj, parameter[['subset']])

  pipe_StatFilterCells(obj)

  obj
}


