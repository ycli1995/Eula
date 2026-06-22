
source("pkg_load.R")

pkg_conf <- yaml::yaml.load_file("Eula_pkgs.yaml")

load_pkg("Eula.CNV", pkg_conf = pkg_conf)
