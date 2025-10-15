
source("site_lib/pkg_load.R")

pkg_conf <- yaml::yaml.load_file("Eula_pkgs.yaml")

load_pkg("Eula.Seurat", pkg_conf = pkg_conf)
