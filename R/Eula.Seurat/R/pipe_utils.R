
#' @export
CalAvgExp <- function(
    object,
    features = NULL,
    group.by = NULL,
    assay = NULL,
    slot = "data",
    method = "average",
    is.expm1 = NULL,
    reverse = FALSE,
    bulk = FALSE,
    outfile = paste0("All", row.type, ".avg_exp.xls"),
    row.type = 'Gene',
    ...
) {
  out <- pseudobulkMatrix(
    object,
    features = features,
    group.by = group.by,
    assay = assay,
    slot = slot,
    method = method,
    is.expm1 = is.expm1,
    reverse = reverse,
    bulk = bulk,
    ...
  )
  out <- AddRowAnnot(object, out, row.type = row.type)
  if (is.null(outfile)) {
    return(out)
  }
  writeTable(out, outfile)
}

#' @importFrom SeuratObject SetAssayData
#' @export
CalAvgPct <- function(
    object,
    features = NULL,
    group.by = NULL,
    assay = NULL,
    slot = "data",
    reverse = FALSE,
    bulk = FALSE,
    outfile = paste0("All", row.type, ".avg_pct.xls"),
    row.type = 'Gene',
    ...
) {
  assay <- assay %||% DefaultAssay(object)
  counts <- GetAssayData(object[[assay]], slot)
  counts[counts > 0] <- 1
  object[[assay]] <- SetAssayData(object[[assay]], "counts", counts)
  out <- pseudobulkMatrix(
    object,
    features = features,
    group.by = group.by,
    assay = assay,
    slot = "counts",
    method = "average",
    is.expm1 = FALSE,
    reverse = reverse,
    bulk = bulk,
    ...
  )
  out <- AddRowAnnot(object, out, row.type = row.type)
  if (is.null(outfile)) {
    return(out)
  }
  writeTable(out, outfile)
}

#' @export
AddRowAnnot <- function(object, mat, row.type = "Gene", ...) {
  gene.ids <- getFeaturesName(object, features = rownames(mat), col = "id")
  gene.names <- getFeaturesName(object, features = gene.ids, col = "name")
  annot <- data.frame(gene.ids, gene.names)
  colnames(annot) <- paste0(row.type, c("ID", "Name"))
  cbind(annot, as.data.frame(mat))
}
