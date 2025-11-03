#' @importFrom tools file_path_as_absolute
#' @importFrom rlang %||%
#' @importFrom vctrs %0%
#' @importFrom dplyr %>%
NULL

#' @export recodeFactor
recodeFactor <- function(x, ...) {
  UseMethod("recodeFactor", x)
}

#' @importFrom dplyr recode recode_factor
#'
#' @export
#' @method recodeFactor default
recodeFactor.default <- function(x, map, keep.orders = FALSE, ...) {
  if (keep.orders) {
    recode_func <- dplyr::recode
  } else {
    recode_func <- dplyr::recode_factor
  }
  map <- unlistMap(map = map)
  x <- x %>%
    as.factor() %>%
    recode_func(!!!map) %>%
    as(Class = class(x))
  x
}

#' @export
#' @method recodeFactor data.frame
recodeFactor.data.frame <- function(x, map, column, keep.orders = FALSE, ...) {
  column <- column[1]
  if (!column %in% colnames(x)) {
    return(x)
  }
  message("Replace entries in column '", column, "':")
  x[, paste0(column, "_old")] <- x[, column]
  x[, column] <- recodeFactor(x[, column], map = map, keep.orders = keep.orders)
  print(table(x[, column], x[, paste0(column, "_old")]))
  x
}

#' @export
YAML_HANDLERS <- list(
  "bool#no" = function(x) if (x %in% c("false", "FALSE")) FALSE else x,
  "bool#yes" = function(x) if (x %in% c("true", "TRUE")) TRUE else x
)

#' @export
YAML_WRITE_LOGICAL <- list(
  logical = function(x) {
    result <- ifelse(x, "TRUE", "FALSE")
    class(result) <- "verbatim"
    result
  }
)

#' @export
maxNChar <- function(str) {
  max(nchar(as.character(str)))
}

#' @export
readRDX <- function(file) {
  con <- gzfile(file)
  on.exit(close(con))
  magic <- readChar(con, 5L, useBytes = TRUE)
  if (grepl("RD[ABX][2-9]\n", magic)) {
    return(get(load(file)))
  }
  readRDS(file)
}

#' @export
Load <- function(file) {
  file <- tools::file_path_as_absolute(file)
  Message("Loading: ", file)
  object <- readRDX(file)
  if (!"version" %in% slotNames(object)) {
    return(object)
  }
  if (grepl('^2', object@version)) {
    checkPackages("Seurat")
    return(Seurat::UpdateSeuratObject(object))
  }
  object
}

#' @importFrom data.table fread
#' @export
readTable <- function(
  file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  row.names = NULL,
  keep.row.names = TRUE,
  data.table = FALSE,
  ...
) {
  df <- data.table::fread(
    file,
    sep = sep,
    header = header,
    stringsAsFactors = stringsAsFactors,
    data.table = data.table,
    ...
  )
  if (is.null(row.names)) {
    return(df)
  }
  rownames(df) <- df[, row.names]
  if (keep.row.names) {
    return(df)
  }
  df[[row.names]] <- NULL
  df
}

#' @export
writeTable <- function(
  x,
  file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  ...
) {
  write.table(
    x,
    file = file,
    quote = quote,
    sep = sep,
    row.names = row.names,
    col.names = col.names,
    ...
  )
}

#' @export
splitArgs <- function(arg, split = ",") {
  if (!is.character(arg)) {
    return(arg)
  }
  unique(unlist(strsplit(arg, split = split)))
}

#' @importFrom yaml yaml.load_file
#' @export
readYAML <- function(file, handlers = NULL, ...) {
  if (length(file) == 0) {
    return(list())
  }
  if (is.na(file)) {
    return(list())
  }
  file <- tools::file_path_as_absolute(file)
  handlers <- handlers %||% YAML_HANDLERS
  parameter <- yaml::yaml.load_file(file, handlers = handlers)
  parameter
}

#' @export
Message <- function(..., type = c("Running", "Warning", "Stopped")) {
  type <- match.arg(arg = type)
  message(paste0("\n[", as.character(Sys.time()), "] [", type, "]: ", ...))
  if (type == "Stopped") {
    q(status = 1)
  }
}

#' @export
loadPackages <- function(packages, lib.loc = .libPaths()) {
  Message('>>>> Start Load Rpackages <<<<')
  Message('---> Check R packages <---')
  libs <- list.files(lib.loc)
  error.packages <- packages[! packages %in% libs]
  if (length(error.packages) > 0) {
    error.packages <- paste(error.packages, collapse = ", ")
    Message("Please install packages: [", error.packages, "]", type = "Stopped")
  }
  Message('---> Loading Rpackages <---')
  for (i in packages) {
    message(paste0(" ~~~> ", i, " <~~~"))
    suppressMessages(library(i, lib.loc = lib.loc, character.only = TRUE))
  }
}

#' @export
unlistMap <- function(map) {
  map <- stack(map)
  map <- setNames(as.character(map$ind), nm = map$values)
  map[unique(names(map))]
}

#' @export
pasteFactors <- function(x, y, collapse = "_", rev = FALSE) {
  if (!is.factor(x = x)) {
    x <- as.factor(x = x)
  }
  if (!is.factor(x = y)) {
    y <- as.factor(x = y)
  }
  lv <- sapply(levels(x), paste0, collapse, levels(y), simplify = TRUE)
  if (rev) {
    lv <- t(lv)
  }
  lv <- as.vector(lv)
  str <- factor(paste0(x, collapse, y), lv)
  str
}

#' @export
filterDotArgs <- function(params, functions) {
  all.args <- lapply(functions, formalArgs)
  all.args <- unique(unlist(all.args))

  keep.args <- intersect(names(params), all.args)
  params[keep.args]
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
fetch_color_table_file <- function(file) {
  rename_table <- readTable(tools::file_path_as_absolute(file), header = FALSE)
  parameter <- list()
  parameter$colors <- fetch_color_dataframe(rename_table)
  return(parameter)
}

#' @export
fetch_rename_dataframe <- function(rename_table) {
  rename_map <- setNames(as.character(rename_table[[2]]), rename_table[[1]])
  rename_map <- as.list(rename_map)
  rename_map <- lapply(rename_map, function(x) unlist(strsplit(x, ",")))
  rename_map
}

#' @export
fetch_color_dataframe <- function(color_table) {
  rename_map <- setNames(as.character(color_table[[2]]), color_table[[1]])
  rename_map <- as.list(rename_map)
  rename_map
}

#' @export
mkdir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}

# find_group_index <- function(
#     mat,
#     metadata,
#     group.by = NULL,
#     groups = NULL,
#     by_gene = NULL,
#     ...
# ) {
#   group.by <- group.by %||% "orig.ident"
#   metadata <- metadata %>% select(!!group.by)
#   if (is.null(groups)) {
#     groups <- list()
#   }
#   for (i in group.by) {
#     if (length(groups[[i]]) == 0) {
#       groups[[i]] <- unique(metadata[, i]) %>%
#         as.character() %>%
#         setNames(nm = .) %>%
#         as.list()
#     }
#   }
#   for (i in group.by) {
#     for (j in unique(metadata[, i])) {
#       metadata <- metadata %>%
#         mutate(!!j := .[[i]] %in% j)
#     }
#     for (j in names(groups[[i]])) {
#       metadata <- metadata %>%
#         mutate(!! j := .[[i]] %in% groups[[i]][[j]])
#     }
#     metadata <- metadata %>%
#       select(!all_of(i))
#   }
#   if (missing(mat)) {
#     return(metadata)
#   }
#   for (i in names(by_gene)) {
#     gene_id <- by_gene[[i]][[1]]
#     gene_thres <- by_gene[[i]][[2]]
#     metadata <- metadata %>%
#       mutate(!! i := data[gene_id, ] >= gene_thres[[1]] & mat[gene_id, ] <= gene_thres[[2]])
#   }
#   return(metadata)
# }

#' @export
excu_shell <- function(..., sep = " ", collapse = NULL) {
  cmd <- paste(..., sep = sep, collapse = collapse)
  cat(cmd, "\n")
  out <- system(cmd)
  if (out == 0) {
    return(invisible(NULL))
  }
  stop("Error when runing shell command: ", cmd)
}

#' @export
check_complete <- function(file) {
  if (!file.exists(file)) {
    stop("No such file: ", file, "\nSomething must not complete correctly.")
  }
}

#' @export
norm_list_param <- function(x, ...) {
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
    x <- readTable(x, header = FALSE, ...)[, 1]
    return(x)
  }
  if (grepl("\\S,\\S", x)) {
    x <- split_args(x)
    x <- sapply(X = x, FUN = norm_list_param, ..., simplify = FALSE)
    return(unlist(x))
  }
  x
}

#' @export
norm_fname <- function(x, ...) {
  gsub("\\/| |\\(|\\)", "_", x)
}
