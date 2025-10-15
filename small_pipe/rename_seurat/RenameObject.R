
me <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
lib.dir <- normalizePath(dirname(me))

args <- commandArgs(TRUE)

file <- args[1]
obj_file <- args[2]
obj_out <- args[3]

if (anyNA(c(file, obj_file, obj_out))) {
  warning("\n  Usage: Rscript ", me, " <parameter.yaml> <obj_file> <obj_out>\n", call. = FALSE, immediate. = TRUE)
  quit()
}

source(tools::file_path_as_absolute(file.path(lib.dir, "Eula.Seurat.R")), chdir = TRUE)

parameter <- readYAML(file)
if (length(parameter[["RenameObject_yaml"]]) > 0) {
  file2 <- tools::file_path_as_absolute(parameter[["RenameObject_yaml"]])
  Message("Read 'RenameObject_yaml': ", file2)
  parameter <- readYAML(file2)
}

obj <- Load(obj_file)

obj <- pipe_RenameObject(obj, parameter)

Message('>>>>> Save object to: ', obj_out)
save(obj, file = obj_out)

outdir <- dirname(obj_out)
Message('>>>>> Output some lists: ', outdir)
WriteRenameConf(obj, outdir)
