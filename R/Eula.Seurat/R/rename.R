
#' @importFrom Eula.utils recodeFactor
#' @export
#' @method recodeFactor Seurat
recodeFactor.Seurat <- function(
    x, map, column,
    keep.orders = FALSE,
    set.ident = FALSE,
    ...
) {
  column = column[1]
  x@meta.data <- recodeFactor(
    x = x@meta.data,
    map = map,
    column = column,
    keep.orders = keep.orders,
    ...
  )
  if (set.ident & column %in% colnames(x@meta.data)) {
    message("Set new column '", column, "' to Idents(object)")
    Idents(x) <- column
  }
  return(x)
}

#' @importFrom Eula.utils checkColorMap fastWarning
#' @export
#' @method checkColorMap Seurat
checkColorMap.Seurat <- function(
    x,
    column,
    colors = NULL,
    tag = NULL,
    misc.name = NULL,
    ...
) {
  if (!column %in% colnames(x@meta.data)) {
    fastWarning(
      "checkColorMap: '", column, "' is not found in ",
      "the input seuratobj@meta.data"
    )
    return(x)
  }
  if (!is.factor(x@meta.data[, column])) {
    x@meta.data[, column] <- as.factor(x@meta.data[, column])
  }
  orig.colors <- NULL
  if (!is.null(misc.name)) {
    orig.colors <- x@misc[[misc.name]]
  }
  orig.colors <- orig.colors %||% x@misc$colors[[column]]
  old.colors <- checkColorMap(
    x@meta.data,
    column = column,
    colors = orig.colors,
    tag = tag
  )
  new.colors <- old.colors
  ## Replace old colors with new colors
  if (any(names(colors) %in% names(new.colors))) {
    cmm.names <- intersect(names(colors), names(new.colors))
    new.colors[cmm.names] <- colors[cmm.names]
    new.colors <- c(new.colors, colors[setdiff(names(colors), cmm.names)])
  }
  new.colors <- checkColorMap(
    x@meta.data,
    column = column,
    colors = new.colors,
    tag = tag,
    ...
  )
  x@misc$colors[[column]] <- new.colors
  if (is.null(misc.name)) {
    return(x)
  }
  x@misc[[misc.name]] <- new.colors
  print(x@misc[[misc.name]])
  x@misc$color2column[[column]] <- misc.name
  x
}

#' @importFrom SeuratObject Idents<- Idents
#' @export
RenameSeurat <- function(
    object,
    map,
    column,
    colors = NULL,
    misc.name = NULL,
    set.ident = FALSE,
    keep.orders = FALSE,
    ...
) {
  object %>%
    recodeFactor(
      map = map,
      column = column,
      set.ident = set.ident,
      keep.orders = keep.orders
    ) %>%
    checkColorMap(
      column = column,
      misc.name = misc.name,
      colors = colors,
      ...
    )
}

#' @export
RenameSeuratColumns <- function(object, map = list()) {
  Message('>>>>> Rename columns in obj@meta.data...')
  for (i in names(map)) {
    if (map[[i]] %in% colnames(object@meta.data)) {
      message("Reset object@meta.data: '", map[[i]], "' -> '", i, "'")
      object@meta.data[, i] <- object@meta.data[, map[[i]]]
    }
  }
  return(object)
}

.COLOR_NAMES <- list(
  "orig.ident" = "color.sample",
  "Groups" = "color.group",
  "seurat_clusters" = "color.cluster",
  "Samples" = "color.sample",
  "Cluster" = "color.cluster",
  "Clusters" = "color.cluster"
)

#' @export
RenameSeuratMetaData <- function(object, parameter = list()) {
  parameter <- fetch_rename_table(parameter)
  column <- intersect(parameter[["column"]], colnames(object[[]]))
  if (length(column) == 0) {
    fastWarning("Cannot found 'column' name. No rename will do.")
    return(object)
  }

  set.ident <- parameter[["set.ident"]] %||% FALSE
  keep.orders <- parameter[["keep.orders"]] %||% FALSE

  color.name <- parameter[["color_name"]] %||% .COLOR_NAMES[[column]]
  color.name <- color.name %||% column

  parameter[["rename_map"]] <- parameter[["rename_map"]] %||%
    if (is.factor(object@meta.data[, column])) {
      object@meta.data[, column] %>%
        levels() %>%
        sapply(function(x) x, simplify = FALSE)
    } else {
      object@meta.data[, column] %>%
        factor() %>%
        levels() %>%
        sapply(function(x) x, simplify = FALSE)
    }

  colors <- unlist(parameter[["colors"]])
  object <- RenameSeurat(
    object,
    map = parameter[["rename_map"]],
    column = column,
    misc.name = color.name,
    colors = colors,
    set.ident = set.ident,
    keep.orders = keep.orders
  )
  object
}

.METADATA_COLUMNS <- list(
  "sample" = "orig.ident",
  "group" = "Groups",
  "cluster" = "seurat_clusters"
)

#' @export
RenameSeuratWrapper <- function(object, parameter = list()) {
  Message('>>>>> Rename Seurat...')
  print(str(object[[]]))
  print(str(object@misc))
  for (i in names(parameter)) {
    message("Rename for '", i, "': ")
    parameter[[i]][["column"]] <- parameter[[i]][["column"]] %||%
      .METADATA_COLUMNS[[i]]
    object <- RenameSeuratMetaData(object, parameter[[i]])
  }
  object
}

#' @export
WriteRenameConf <- function(
    object,
    outdir = getwd(),
    columns = c("Samples", "Groups", "Clusters")
) {
  columns <- intersect(columns, colnames(object[[]]))
  for (i in columns) {
    labels <- levels(as.factor(object[[]][, i]))
    write(labels, file.path(outdir, paste0(i, ".list")))

    rename.df <- data.frame(V1 = labels, V2 = labels)
    colors <- object@misc$colors[[i]]
    if (length(colors) > 0) {
      rename.df$V3 <- colors[labels]
    }
    writeTable(
      rename.df,
      file.path(outdir, paste0(i, "_conf.xls")),
      col.names = FALSE
    )
  }
  invisible(NULL)
}


