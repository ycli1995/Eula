#' @include pipe.R
NULL

#' @export
run_seurat_main <- function(method, file) {
  run_func <- get(paste0("pipe_", method, "_R"))

  pipeMsg("##### Found pipeline: ", method)

  pipeMsg("##### Reading parameter: ", file)
  parameter <- readYAML(file)

  pipeMsg("##### Running pipeline")
  run_func(parameter)

  pipeMsg("##### Finish pipeline")
  invisible(NULL)
}

#' @export
run_DoubletFinder <- function(file, indir, sample, outdir) {
  pipeMsg("##### Reading parameter: ", file)
  params <- readYAML(file)

  pipeMsg("##### Format parameters")
  params[['MakeSeuratObj']] <- list(
    names = sample,
    paths = indir
  )
  params[['outdir']] <- outdir

  pipeMsg("##### Running DoubletFinder")
  pipe_DoubletFinder_R(params)

  pipeMsg("##### Finish pipeline")
  invisible(NULL)
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
  pipeMsg('>>>>> Save object to: ', outfile)
  save(obj, file = outfile)

  pipeMsg('>>>>> Output some lists: ', outdir)
  WriteRenameConf(obj, outdir)

  invisible(obj)
}

#' @importFrom Eula.utils checkKeys
#' @export
pipe_MakeSeuratObj_R <- function(params) {

  outfile <- params[['obj_out']] %||% file.path(getwd(), "obj.Rda")
  outdir <- dirname(outfile)

  checkKeys(params, "MakeSeuratObj")

  obj <- pipe_MakeSeuratObj(params[["MakeSeuratObj"]])

  obj <- pipe_StatFeaturePercentage(obj, params[["StatFeaturePercentage"]])

  obj <- pipe_AddMetaData(obj, params = params[['AddMetaData']])

  obj <- RenameSeuratColumns(obj, params[['rename_colnames']])

  obj <- pipe_RenameSeurat(obj, params[['rename']])

  obj <- pipe_UpdateSeurat(obj, params[['UpdataSeurat']])

  obj <- CheckMySeuratObj(obj, column.map = params[['default_colnames']])

  obj@misc$colData <- obj[[]]

  params[["StatCellQuality"]][['outdir']] <- file.path(outdir, "before")
  pipe_StatCellQuality(obj, params[["StatCellQuality"]])

  obj <- pipe_SubsetObject(obj, params[['subset']])

  mkdir(outdir)
  params[['StatFilterCells']][['outdir']] <- outdir
  pipe_StatFilterCells(obj, params[['StatFilterCells']])

  params[["StatCellQuality"]][['outdir']] <- file.path(outdir, "after")
  pipe_StatCellQuality(obj, params[["StatCellQuality"]])

  pipeMsg('Save object to: ', outfile)
  save(obj, file = outfile)

  pipeMsg('Output some lists: ', outdir)
  WriteRenameConf(obj, outdir)

  invisible(obj)
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
  pipeMsg('>>>>> Save object to: ', outfile)
  save(obj, file = outfile)

  pipeMsg('>>>>> Output some lists: ', outdir)
  WriteRenameConf(obj)

  pipe_PostClustering(obj, params[['PostClustering']])

  invisible(obj)
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

#' @export
pipe_DoubletFinder_R <- function(params = list(), ...) {
  outdir <- params[['outdir']] %||% getwd()

  pipe_MultiCore(params)

  checkKeys(params, "MakeSeuratObj")

  obj <- pipe_MakeSeuratObj(params[["MakeSeuratObj"]])

  mkdir(outdir, chdir = TRUE)
  obj <- pipe_DoubletFinder(obj, params[['DoubletFinder']])

  invisible(obj)
}

