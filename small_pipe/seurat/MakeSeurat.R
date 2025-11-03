
me <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
lib.dir <- normalizePath(dirname(me), mustWork = TRUE)

args <- commandArgs(TRUE)

file <- args[1]

if (anyNA(file)) {
  warning(
    "\n  Usage: Rscript ", me, " <parameter.yaml>\n",
    call. = FALSE, immediate. = TRUE
  )
  quit()
}

lib.R <- normalizePath(file.path(lib.dir, "Eula.Seurat.R"), mustWork = TRUE)
source(lib.R, chdir = TRUE)

parameter <- readYAML(file)

obj <- pipe_MakeSeuratObj_R(parameter)
