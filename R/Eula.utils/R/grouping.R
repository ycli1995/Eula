#' @importFrom vctrs %0%
NULL

#' @export
group_cells <- function(cell.groups, cells = NULL) {
  if (!is.list(cell.groups)) {
    stop("'cell.groups' must be a list with grouping cells.")
  }
  cells <- cells %||% unique(unlist(cell.groups))
  out <- data.frame(row.names = cells)
  for (i in seq_along(cell.groups)) {
    group.name <- names(cell.groups)[i]
    if (length(group.name) == 0) {
      stop("Group name is missing for the given cell set.")
    }
    out[[group.name]] <- rownames(out) %in% cell.groups[[i]]
  }
  out
}

#' @export
group_cells_by_meta <- function(metadata, groups) {
  if (!is.list(groups)) {
    stop("'groups' must be a list with grouping info.")
  }
  out <- data.frame(row.names = rownames(metadata))
  for (i in seq_along(groups)) {
    group.name <- names(groups)[i]
    if (length(group.name) == 0) {
      stop("Group name is missing.")
    }
    if (!is.list(groups[[i]])) {
      stop("'groups[[\"", group.name, "\"]]' must be a list.")
    }
    group.info <- groups[[i]]
    tmp <- list()
    for (j in seq_along(group.info)) {
      name <- names(group.info)[j]
      if (length(name) == 0) {
        stop("Name of meta data is missing for group '", group.name, "'.")
      }
      if (!name %in% colnames(metadata)) {
        stop("'", name, "' not found in meta data.")
      }
      tmp[[j]] <- metadata[[name]] %in% groups[[i]][[j]]
    }
    out[[group.name]] <- Reduce("&", tmp)
  }
  out
}

#' @export
group_cells_by_features <- function(data, groups) {
  if (!is.list(groups)) {
    stop("'groups' must be a list with grouping info.")
  }
  out <- data.frame(row.names = colnames(data))
  for (i in seq_along(groups)) {
    group.name <- names(groups)[i]
    if (length(group.name) == 0) {
      stop("Group name is missing.")
    }
    if (!is.list(groups[[i]])) {
      stop("'groups[[\"", group.name, "\"]]' must be a list.")
    }
    group.info <- groups[[i]]
    tmp <- list()
    for (j in seq_along(group.info)) {
      name <- names(group.info)[j]
      if (length(name) == 0) {
        stop("Feature name is missing for group '", group.name, "'.")
      }
      if (!name %in% rownames(data)) {
        stop("Feature '", name, "' not found in the input data.")
      }
      if (length(group.info[[j]]) != 2) {
        stop("Length of grouping info by feature must be 2.")
      }
      tmp[[j]] <- data[name, ] >= group.info[[j]][[1]] &
        data[name, ] <= group.info[[j]][[2]]
    }
    out[[group.name]] <- Reduce("&", tmp)
  }
  out
}
