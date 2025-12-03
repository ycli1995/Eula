
me <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
lib.dir <- normalizePath(dirname(me))

args <- commandArgs(TRUE)

file <- args[1]
indir <- args[2]
sample <- args[3]
outdir <- args[4]

if (anyNA(c(file, indir, sample, outdir))) {
  warning(
    "\n  Usage: Rscript ", me, "<parameter.yaml> <indir> <sample> <outdir>\n",
    call. = FALSE, immediate. = TRUE
  )
  quit()
}

lib.R <- normalizePath(file.path(lib.dir, "Eula.Seurat.R"), mustWork = TRUE)
source(lib.R, chdir = TRUE)

run_DoubletFinder(file, indir, sample, outdir)
