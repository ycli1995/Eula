
#' @importFrom Eula.utils pasteFactors pseudobulkMatrix
#' @importFrom SeuratObject DefaultAssay FetchData GetAssayData JoinLayers
#' Version
#' @export
#' @method pseudobulkMatrix Seurat
pseudobulkMatrix.Seurat <- function(
    object,
    group.by = NULL,
    features = NULL,
    assay = NULL,
    slot = "data",
    method = c("aggregate", "average"),
    is.expm1 = NULL,
    reverse = FALSE,
    bulk = FALSE,
    sparse = FALSE,
    ...
) {
  assay <- assay %||% DefaultAssay(object)
  is.expm1 <- is.expm1 %||% slot == "data"
  group.by <- group.by %||% "ident"

  if (Version(object) >= package_version("5.0.0")) {
    object <- JoinLayers(object)
  }
  data <- GetAssayData(object[[assay]], slot)
  if (!is.null(features)) {
    data <- data[features, , drop = FALSE]
  }
  if (is.expm1) {
    data <- expm1(data)
  }
  cell.groups <- FetchData(object, vars = group.by)
  cell.groups <- Reduce(pasteFactors, cell.groups)
  pseudobulkMatrix(
    object = data,
    cell.groups = cell.groups,
    method = method,
    reverse = reverse,
    bulk = bulk,
    sparse = sparse,
    ...
  )
}
