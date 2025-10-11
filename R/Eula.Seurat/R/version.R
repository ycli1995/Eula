
#' @importFrom SeuratObject Images UpdateSeuratObject Version
#' @export
UpdateReductions5 <- function(object) {
  if (Version(object) >= 5) {
    return(object)
  }
  old.reds <- object@reductions
  keep <- grep("_", names(old.reds), value = TRUE, invert = TRUE)
  old.reds <- old.reds[keep]
  object@reductions <- old.reds
  if ("images" %in% slotNames(object)) {
    old_images <- Images(object)
  }
  object <- UpdateSeuratObject(object)
  if ("images" %in% slotNames(object)) {
    object@images <- object@images[old_images]
  }
  return(object)
}

#' @importFrom SeuratObject Assays JoinLayers VariableFeatures
#' @export
UpdateAssay5 <- function(object) {
  for (i in SeuratObject::Assays(object)) {
    object[[i]] <- as(object[[i]], "Assay5")
    object[[i]] <- JoinLayers(object[[i]])
    print(str(object[[i]]))
    print(str(VariableFeatures(object[[i]])))
  }
  return(object)
}

#' @export
UpdateSeurat5 <- function(object) {
  object <- UpdateAssay5(object)
  object <- UpdateReductions5(object)
  return(object)
}

#' @importFrom SeuratObject Images UpdateSeuratObject Version
#' @export
UpdateSeurat3 <- function(object) {
  if (Version(object) >= package_version("3.2")) {
    return(object)
  }
  old.reds <- object@reductions
  keep <- grep("_", names(old.reds), value = TRUE, invert = TRUE)
  old.reds <- old.reds[keep]
  object@reductions <- old.reds

  for (i in names(object@assays)) {
    if (is.logical(GetAssayData(object@assays[[i]], "scale.data"))) {
      object@assays[[i]] <- SetAssayData(
        object@assays[[i]],
        "scale.data",
        new.data = matrix(0, 0, 0)
      )
    }
  }

  old.images <- tryCatch(
    Images(object),
    error = function(e) {
      warning("No 'images' slot.", call. = FALSE, immediate. = TRUE)
      character()
    }
  )

  object <- UpdateSeuratObject(object)

  if (length(old.images) > 0) {
    object@images <- object@images[old.images]
  }
  object
}

#' @importFrom Seurat CreateAssayObject GetAssayData VariableFeatures<-
#' @export
Assay5to4 <- function(assay, assay.orig = "RNA", ...) {
  if (inherits(assay, "Assay")) {
    return(assay)
  }
  new.assay <- CreateAssayObject(
    counts = GetAssayData(assay, "counts"),
    min.cells = -1,
    min.features = -1
  )
  new.assay@data <- GetAssayData(assay, "data")
  scale.data <- GetAssayData(assay, "scale.data")
  if (!is.null(scale.data)) {
    new.assay@scale.data <- scale.data
  }
  mf <- assay[[]]
  rownames(mf) <- rownames(assay)
  new.assay@meta.features <- mf
  VariableFeatures(new.assay) <- VariableFeatures(assay)
  new.assay@assay.orig <- assay.orig
  print(str(new.assay))
  return(new.assay)
}

#' @importFrom SeuratObject DefaultAssay GetImage GetTissueCoordinates Key
#' @importFrom Seurat scalefactors ScaleFactors
#' @importClassesFrom Seurat VisiumV1
#' @export
VisiumV2toV1 <- function(imageObject) {
  image <- GetImage(imageObject, "raw")
  scale.factors <- ScaleFactors(imageObject)
  unnormalized.radius <- scale.factors$fiducial * scale.factors$lowres
  spot.radius <- unnormalized.radius/max(dim(image))

  coordinates <- GetTissueCoordinates(imageObject, scale = NULL) %>%
    dplyr::rename(imagerow = "x", imagecol = "y") %>%
    dplyr::mutate(tissue = 1L, row = imagerow, col = imagecol, .before = 1) %>%
    dplyr::select(!cell)

  return(new(
    Class = "VisiumV1",
    image = image,
    scale.factors = scalefactors(
      spot = scale.factors$hires,
      fiducial = scale.factors$fiducial,
      hires = scale.factors$hires,
      lowres = scale.factors$lowres
    ),
    coordinates = coordinates,
    spot.radius = spot.radius,
    assay = DefaultAssay(imageObject),
    key = Key(imageObject)
  ))
}

#' @export
ValidGraph4 <- function(graph, assay.used = "RNA", ...) {
  if (length(graph@assay.used) > 0) {
    return(graph)
  }
  graph@assay.used <- assay.used
  return(graph)
}

#' @importFrom SeuratObject Graphs
#' @export
Seurat5to4 <- function(obj, ...) {
  for (i in SeuratObject::Assays(obj)) {
    obj[[i]] <- Assay5to4(obj[[i]], assay.orig = i)
  }
  for (i in SeuratObject::Graphs(obj)) {
    obj[[i]] <- ValidGraph4(obj[[i]], assay.used = DefaultAssay(obj))
  }
  obj@version <- package_version("4.3.0")
  return(obj)
}


