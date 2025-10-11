
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

file <- tools::file_path_as_absolute(file)
obj_file <- tools::file_path_as_absolute(obj_file)

parameter <- readYAML(file)
if (length(parameter[["RenameObject_yaml"]]) > 0) {
  file2 <- tools::file_path_as_absolute(parameter[["RenameObject_yaml"]])
  Message("Read 'RenameObject_yaml': ", file2)
  parameter <- readYAML(file2)
}

parameter$scale.data <- parameter$scale.data %||% TRUE
parameter$misc.counts <- parameter$misc.counts %||% TRUE
parameter$Assay5 <- parameter$Assay5 %||% FALSE

obj <- Load(obj_file)

parameter$assays <- parameter$assays %||% Assays(obj)
parameter$assays <- intersect(parameter$assays, Assays(obj))

print(parameter)

Message('>>>>> Rename columns in obj@meta.data...')
obj <- RenameSeuratColumns(obj, parameter$rename_colnames)

Message('>>>>> Rename Seurat...')
obj <- RenameSeuratWrapper(obj, parameter$rename)

obj <- UpdateSeurat3(obj)
if (packageVersion("Seurat") >= 5) {
  obj <- UpdateReductions5(obj)
  if (parameter$Assay5) {
    obj <- UpdateAssay5(obj)
  }
}

Message('>>>>> Check default columns...')
obj <- CheckSeuratMetaData(obj, parameter$default_colnames)

Message('>>>>> Subset Seurat object...')
obj <- SubsetObjectWrapper(obj, parameter$subset)

Message('>>>>> Finally check the Seurat object...')
obj <- ShrinkSeuratObject(obj, parameter$assays, parameter$scale.data, parameter$misc.counts)

obj
str(obj@meta.data)
print(obj@misc$colors)

Message('>>>>> Save object to: ', obj_out)
save(obj, file = obj_out)

Message('>>>>> Output some lists: ', dirname(obj_out))
setwd(dirname(obj_out))
for (i in c("Samples", "Groups", "Clusters")) {
  write(levels(obj[[]][, i]), paste0(i, ".list"))
  colors <- obj@misc$colors[[i]]
  rename.df <- data.frame(
    V1 = names(colors),
    V2 = names(colors),
    V3 = colors
  )
  writeTable(rename.df, paste0(i, "_conf.xls"), col.names = FALSE)
}



