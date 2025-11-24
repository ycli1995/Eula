#' @importFrom tools file_path_as_absolute
NULL

#' @export
normalizeFileName <- function(x, ...) {
  gsub("\\/| |\\(|\\)", "_", x)
}

#' @export
YAML.HANDLERS.READ <- list(
  "bool#no" = function(x) if (x %in% c("false", "FALSE")) FALSE else x,
  "bool#yes" = function(x) if (x %in% c("true", "TRUE")) TRUE else x
)

#' @export
YAML.HANDLERS.WRITE <- list(
  logical = function(x) {
    result <- ifelse(x, "TRUE", "FALSE")
    class(result) <- "verbatim"
    result
  }
)

#' @importFrom yaml yaml.load_file
#' @export
readYAML <- function(file, handlers = NULL, ...) {
  if (length(file) == 0) {
    return(list())
  }
  if (is.na(file)) {
    return(list())
  }
  file <- file_path_as_absolute(file)
  handlers <- handlers %||% YAML.HANDLERS.READ
  yaml::yaml.load_file(file, handlers = handlers, ...)
}

#' @importFrom yaml write_yaml
#' @export
writeYAML <- function(x, file, handlers = NULL, ...) {
  handlers <- handlers %||% YAML.HANDLERS.WRITE
  .check_dir(file)
  yaml::write_yaml(x = x, file = file, handlers = handlers, ...)
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
  .check_dir(file)
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

#' @importFrom tools file_path_as_absolute
#' @export
readRDX <- function(file) {
  file <- file_path_as_absolute(file)
  verboseMsg("Loading: ", file)
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
  readRDX(file)
}

#' @export
mkdir <- function(path, chdir = FALSE) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  if (chdir) {
    setwd(path)
  }
  invisible(NULL)
}

.check_dir <- function(file) {
  dir <- dirname(file)
  if (!dir.exists(dir)) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  }
  invisible(NULL)
}
