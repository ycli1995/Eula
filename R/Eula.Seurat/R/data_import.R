
#' @importFrom Matrix readMM
.Read10xPeak <- function(
    data_path,
    mtx_file = 'matrix.mtx',
    peaks_file = 'peaks.bed',
    barcodes_file = 'barcodes.tsv'
) {
  barcodes <- readLines(file.path(data_path, barcodes_file))
  peaks <- data.table::fread(
    file.path(data_path, peaks_file),
    sep = "\t",
    header = FALSE,
    data.table = FALSE,
    stringsAsFactors = FALSE
  )
  peaks <- paste0(peaks$V1, ":", peaks$V2, "-", peaks$V3)
  counts <- Matrix::readMM(file.path(data_path, mtx_file))
  colnames(counts) <- barcodes
  rownames(counts) <- peaks
  counts
}

#' @importFrom Seurat Read10X
.Read10xDir <- function(
    data_path,
    use.names = TRUE,
    remove_suffix = FALSE
) {
  if (file.exists(file.path(data_path, "peaks.bed"))) {
    return(.Read10xPeak(data_path, remove_suffix = remove_suffix))
  }
  Seurat::Read10X(
    data_path,
    gene.column = ifelse(use.names, 2, 1),
    strip.suffix = remove_suffix
  )
}

.Read10xH5 <- function(
    data_path,
    use.names = TRUE,
    remove_suffix = FALSE
) {
  mat <- Seurat::Read10X_h5(data_path, use.names = use.names)
  if (!is.list(mat)) {
    mat <- list(mat)
  }
  if (remove_suffix) {
    for (i in seq_along(mat)) {
      colnames(mat[[i]]) <- gsub("-[0-9]+$", "", colnames(mat[[i]]))
    }
  }
  if (length(mat) == 1) {
    return(mat[[1]])
  }
  mat
}

.Read10xMat <- function(
    data_path,
    use.names = TRUE,
    remove_suffix = FALSE,
    assay = NULL
) {
  data_path <- normalizePath(data_path[1], mustWork = TRUE)
  if (dir.exists(data_path)) {
    mat <- .Read10xDir(
      data_path = data_path,
      use.names = use.names,
      remove_suffix = remove_suffix
    )
  } else if (hdf5r::is_hdf5(data_path)) {
    mat <- .Read10xH5(
      data_path = data_path,
      use.names = use.names,
      remove_suffix = remove_suffix
    )
  }
  if (!is.list(mat)) {
    return(mat)
  }
  names(mat)[names(mat) == "Gene Expression"] <- "RNA"
  names(mat)[names(mat) == "Peaks"] <- "ATAC"
  if (!is.null(assay)) {
    mat <- mat[assay]
  }
  if (length(mat) == 1) {
    return(mat[[1]])
  }
  return(mat)
}

.Read10xFeatures <- function(data_path, use.names = TRUE, assay = NULL) {
  data_path <- normalizePath(data_path[1], mustWork = TRUE)
  if (dir.exists(data_path)) {
    features <- data.table::fread(
      file.path(data_path, "features.tsv.gz"),
      sep = "\t",
      header = FALSE,
      data.table = FALSE,
      stringsAsFactors = FALSE
    )[, 1:3]
  } else if (hdf5r::is_hdf5(data_path)) {
    features <- hdf5r.Extra::h5Read(data_path, "matrix/features")
    features <- data.frame(
      V1 = features$id,
      V2 = features$name,
      V3 = features$feature_type
    )
  }
  colnames(features) <- c("id", "name", "type")
  if (use.names) {
    rownames(features) <- make.unique(features$GeneName)
  } else {
    rownames(features) <- features$GeneID
  }
  if (length(unique(features$Type)) < 2) {
    return(features)
  }
  features <- split(features, f = features$Type)
  names(features)[names(features) == "Gene Expression"] <- "RNA"
  names(features)[names(features) == "Peaks"] <- "ATAC"
  if (!is.null(assay)) {
    features <- features[assay]
  }
  if (length(features) == 1) {
    return(features[[1]])
  }
  return(features)
}

#' @importFrom SeuratObject Cells CreateSeuratObject GetAssayData
#' @export
MakeSeuratObj <- function(
    data.paths,
    data.names,
    assay = "RNA",
    name.list = NULL,
    ...
) {
  mats <- list()
  features <- list()
  for (i in seq(data.names)) {
    mats[[i]] <- .Read10xMat(data.paths[i], use.names = FALSE, assay = assay)
    colnames(mats[[i]]) <- paste0(data.names[i], "_", colnames(mats[[i]]))

    features[[i]] <- rownames(mats[[i]])
    mats[[i]] <- CreateSeuratObject(
      mats[[i]],
      assay = assay,
      project = data.names[i],
      min.cells = -1,
      min.features = -1
    )
  }
  # Merge objects, get pdata and matrix
  mats <- MergeObject(mats)
  pdata <- mats[[]]
  mats <- GetAssayData(mats[[assay]], "counts")

  # Match all.features and seurat row names
  all.features <- Reduce(union, features)
  names(all.features) <- gsub("_", "-", all.features)

  # Add fdata
  all.features <- all.features[rownames(mats)]
  fdata <- ReadNameList(all.features, name_list = name_list)
  fdata <- fdata[rownames(mats), ]

  # Change to gene symbols
  rownames(mats) <- rownames(fdata) <- make.unique(fdata$name)
  rownames(mats) <- rownames(fdata) <- gsub("_", "-", rownames(fdata))
  fdata$seurat_id <- rownames(fdata)

  object <- CreateSeuratObject(
    mats,
    assay = assay,
    min.cells = -1,
    min.features = -1
  )
  object$orig.ident <- pdata[Cells(object), "orig.ident"]
  object
}

#' @export
ReadNameList <- function(features, name.list = NULL) {
  if (is.null(name.list)) {
    return(data.frame(
      id = features,
      merge_name = features,
      name = features,
      row.names = features
    ))
  }
  name.list <- normalizePath(name.list, mustWork = TRUE)
  fdata <- readTable(name.list, header = FALSE, quote = '"')
  if (ncol(fdata) > 4) {
    fdata <- fdata[, 1:4]
  }
  rownames(fdata) <- fdata[, 1]
  colnames(fdata) <- c("id", "merge_name", "name", "type")
  fdata <- fdata[features, , drop = FALSE]
  if (length(names(features)) == nrow(fdata)) {
    rownames(fdata) <- names(features)
  }
  fdata
}

#' @export
MergeObject <- function(object.list, ...) {
  if (!is.list(object.list)) {
    return(object.list)
  }
  if (length(object.list) == 1) {
    return(object.list[[1]])
  }
  object <- merge(object.list[[1]], object.list[-1], ...)
  if (Version(object) < 5) {
    return(object)
  }
  object <- SeuratObject::JoinLayers(object)
  object
}

