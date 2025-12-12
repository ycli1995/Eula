
#' @export
splitArgs <- function(arg, split = ",") {
  if (!is.character(arg)) {
    return(arg)
  }
  unique(unlist(strsplit(arg, split = split)))
}

#' @export
filterArgs <- function(args, functions) {
  all.args <- lapply(functions, formalArgs)
  all.args <- unique(unlist(all.args))

  keep.args <- intersect(names(args), all.args)
  args[keep.args]
}

#' Full-fill a list of arguments with defaults
#'
#' @param args A list of arguments.
#' @param defaults A list of default arguments.
#' @param ... `r .dot_param`
#'
#' @details
#' Arguments as `NULL` (no default value) in `defaults` are also kept and will
#' be assigned to `args`.
#'
#' @returns
#' A full-filled `args` containing all key-value present in `defaults`.
#'
#' @export
getDefaultArgs <- function(args = list(), defaults = list(), ...) {
  if (!is.list(defaults)) {
    stop("'defaults' must be a list.")
  }
  null.args <- c()
  for (i in names(defaults)) {
    args[[i]] <- args[[i]] %||% defaults[[i]]
    if (is.null(args[[i]])) {
      null.args <- c(null.args, i)
    }
  }
  null.args <- sapply(null.args, function(x) NULL, simplify = FALSE)
  c(args, null.args)
}

#' @export
getArgList <- function(x, sep = "\t", split = ",", ...) {
  if (length(x) == 0) {
    return(x)
  }
  if (length(x) > 1) {
    return(x)
  }
  if (!is.character(x)) {
    return(x)
  }
  if (file.exists(x)) {
    x <- readTable(x, sep = sep, header = FALSE, ...)[, 1]
    return(x)
  }
  pattern <- paste0("\\S", split, "\\S")
  if (grepl(pattern, x)) {
    x <- splitArgs(x, split = split)
    x <- sapply(X = x, FUN = getArgList, ..., simplify = FALSE)
    return(unlist(x))
  }
  x
}

#' @export
runCmd <- function(..., sep = " ", collapse = NULL) {
  cmd <- paste(..., sep = sep, collapse = collapse)
  cat(cmd, "\n")
  out <- system(cmd)
  if (out == 0) {
    return(invisible(NULL))
  }
  stop("Error when runing shell command: ", cmd)
}


#' @export
fetch_rename_table <- function(parameter) {
  if (length(parameter$rename_file) == 0) {
    return(parameter)
  }
  rename_table <- readTable(
    tools::file_path_as_absolute(parameter$rename_file),
    header = FALSE
  )
  parameter$rename_map <- fetch_rename_dataframe(rename_table)
  if (ncol(rename_table) < 3) {
    return(parameter)
  }
  rename_table <- rename_table[, -2]
  parameter$colors <- fetch_color_dataframe(rename_table)
  parameter
}

#' @export
fetch_color_table <- function(parameter) {
  if (length(parameter$color_file) == 0) {
    return(parameter)
  }
  rename_table <- readTable(
    tools::file_path_as_absolute(parameter$color_file),
    header = FALSE
  )
  parameter$colors <- fetch_color_dataframe(rename_table)
  parameter
}

#' @importFrom tools file_path_as_absolute
#' @export
fetch_rename_table_file <- function(file) {
  rename_table <- readTable(tools::file_path_as_absolute(file), header = FALSE)
  parameter <- list()
  parameter$rename_map <- fetch_rename_dataframe(rename_table)
  if (ncol(rename_table) < 3) {
    return(parameter)
  }
  rename_table <- rename_table[, -2]
  parameter$colors <- fetch_color_dataframe(rename_table)
  parameter
}

#' @export
fetch_rename_dataframe <- function(rename_table) {
  rename_map <- setNames(as.character(rename_table[[2]]), rename_table[[1]])
  rename_map <- as.list(rename_map)
  rename_map <- lapply(rename_map, function(x) unlist(strsplit(x, ",")))
  rename_map
}

#' @importFrom tools file_path_as_absolute
#' @export
fetch_color_table_file <- function(file) {
  rename_table <- readTable(tools::file_path_as_absolute(file), header = FALSE)
  parameter <- list()
  parameter$colors <- fetch_color_dataframe(rename_table)
  parameter
}

#' @export
fetch_color_dataframe <- function(color_table) {
  rename_map <- setNames(as.character(color_table[[2]]), color_table[[1]])
  rename_map <- as.list(rename_map)
  rename_map
}
