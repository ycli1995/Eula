
#' @export
pipe_RenameObject <- function(obj, parameter = list()) {

  parameter$scale.data <- parameter$scale.data %||% TRUE
  parameter$misc.counts <- parameter$misc.counts %||% TRUE
  parameter$Assay5 <- parameter$Assay5 %||% FALSE

  parameter$assays <- parameter$assays %||% Assays(obj)
  parameter$assays <- intersect(parameter$assays, Assays(obj))

  print(parameter)

  Message('>>>>> Rename columns in obj@meta.data...')
  obj <- RenameSeuratColumns(obj, parameter$rename_colnames)

  Message('>>>>> Rename Seurat...')
  obj <- RenameSeuratWrapper(obj, parameter$rename)

  Message('>>>>> Check Seurat version...')
  obj <- UpdateSeuratAll(obj, assay5 = parameter$Assay5)

  Message('>>>>> Check default columns...')
  obj <- CheckMySeuratObj(obj, parameter$default_colnames)

  Message('>>>>> Subset Seurat object...')
  obj <- SubsetObjectWrapper(obj, parameter$subset)

  Message('>>>>> Shrinking Seurat object...')
  obj <- ShrinkSeuratObject(
    object = obj,
    assays = parameter$assays,
    scale.data = parameter$scale.data,
    misc.counts = parameter$misc.counts
  )

  Message('>>>>> Finally check...')
  print(obj)
  print(str(obj@meta.data))
  print(obj@misc$colors)

  obj
}
