
me <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
lib.dir <- normalizePath(dirname(me))

args <- commandArgs(TRUE)

method <- args[1]
file <- args[2]

if (anyNA(c(method, file))) {
  warning(
    "\n  Usage: Rscript ", me, "<method> <parameter.yaml>\n",
    call. = FALSE, immediate. = TRUE
  )
  quit()
}

lib.R <- normalizePath(file.path(lib.dir, "Eula.Seurat.R"), mustWork = TRUE)
source(lib.R, chdir = TRUE)

run_func <- get(paste0("pipe_", method, "_R"))

Message("##### Found pipeline: ", method)

Message("##### Reading parameter: ", file)
parameter <- readYAML(file)
print(str(parameter))

Message("##### Running pipeline")
run_func(parameter)

Message("##### Finish pipeline")
