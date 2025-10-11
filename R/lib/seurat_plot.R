
Seurat_DimPlot <- function(
    object,
    dims = 1:2,
    reduction = NULL,
    cells = NULL,
    
    group.by = NULL,
    split.by = NULL,
    shape.by = NULL,
    
    shuffle = FALSE,
    seed = 42,
    
    label = FALSE,
    label.size = 4,
    label.color = "black",
    label.box = FALSE,
    
    colors = NULL,
    pt.size = NULL,
    nrow = NULL,
    combine = TRUE,
    raster = NULL,
    raster.dpi = c(512, 512),
    labs.args = list(),
    theme = NULL,
    ...
) {
  reduction <- reduction %||% DefaultDimReduc(object)
  cells <- cells %||% Cells(object, assay = DefaultAssay(object[[reduction]]))
  dims <- paste0(Key(object[[reduction]]), dims)
  
  orig.groups <- group.by
  group.by <- group.by %||% "ident"
  split.by <- split.by[1]
  shape.by <- shape.by[1]
  
  if (length(colors) > 0) {
    if (length(group.by) > 1) {
      if (!is.list(colors)) {
        stop(
          "To use multiple 'group.by', 'colors' must be a list with \n",
          "corresponeding 'group.by' names."
        )
      }
    } else {
      colors <- list(colors = colors)
      names(colors) <- group.by
    }
  } else {
    colors <- list()
  }
  vars <- c(dims, group.by, split.by, shape.by)
  data <- FetchData(object, vars = vars, cells = cells, clean = "project")
  
  split.by <- intersect(split.by, colnames(data))
  group.by <- intersect(group.by, colnames(data))
  shape.by <- intersect(shape.by, colnames(data))
  
  ncol <- NULL
  if (length(split.by) > 0) {
    data[["split"]] <- as.factor(data[[split.by]])
    if (length(group.by) > 1) {
      nrow <- 1
      ncol <- 1
    }
  }
  
  if (length(shape.by) > 0) {
    data[['shape']] <- data[[shape.by]]
  }
  
  colnames(data)[1:2] <- c("x", "y")
  if (isTRUE(shuffle)) {
    set.seed(seed)
    data <- data[sample(seq_len(nrow(data))), , drop = FALSE]
  }
  
  plots <- list()
  for (i in group.by) {
    use.cols <- c("x", "y", i, "split", "shape")
    use.cols <- intersect(use.cols, colnames(data))
    tmp.data <- data[, use.cols, drop = FALSE]
    
    colnames(tmp.data)[3] <- "color"
    tmp.data[['color']] <- as.factor(tmp.data[['color']])
    
    if (label) {
      label.data <- get_median_position(data, group.by = c("color", "split"))
      label.data$label <- label.data$color
    } else {
      label.data <- NULL
    }
    
    p <- single_dim_plot(
      data = tmp.data,
      colors = colors[[i]],
      label.data = label.data,
      label.repel = label.repel,
      label.size = label.size,
      label.color = label.color,
      label.box = label.box,
      raster = raster,
      raster.dpi = raster.dpi,
      labs.args = labs.args,
      theme = theme,
      nrow = nrow,
      ...
    )
    if ("split" %in% colnames(tmp.data)) {
      p <- p + facet_theme()
    }
    plots[[i]] <- p
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = ncol)
  }
  return(plots)
}

