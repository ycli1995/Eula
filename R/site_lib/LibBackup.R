me <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
args <- commandArgs(TRUE)

infile <- args[1]
outdir <- args[2]

if (anyNA(infile, outdir)) {
  warning("\n  Usage: Rscript ", me, " <Eula.xxx.R> <outdir>\n", call. = FALSE, immediate. = TRUE)
  quit()
}

infile <- tools::file_path_as_absolute(infile)
setwd(dirname(infile))
source(infile, chdir = TRUE)

name <- basename(infile)
name <- gsub("\\.R$", "", name)

backup_lib(name, lib_conf, outdir)
backup_conf(name, lib_conf, outdir)
