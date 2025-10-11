
.warn <- function(..., call. = FALSE, immediate. = TRUE) {
  warning(..., call. = call., immediate. = immediate.)
}

get_depend_tree <- function(name, conf) {
  name <- name[1]
  queue <- c()
  for (i in conf[[name]]$import) {
    queue <- c(queue, get_depend_tree(i, conf))
  }
  queue <- c(name, queue)
  return(rev(unique(queue)))
}

get_imported_by <- function(name, conf) {
  dep.map <- suppressWarnings(
    stack(sapply(names(conf), function(x) conf[[x]]$import))
  )
  imported.by <- dep.map$ind[dep.map$values %in% name]
  unique(imported.by)
}

get_source_path <- function(conf) {
  dependencies <- c()
  queue <- unlist(conf$dependencies)
  queue <- union(queue, names(conf$dependencies))
  while (length(queue)) {
    i <- queue[1]
    queue <- queue[-1]
    if (!i %in% dependencies) {
      dependencies <- c(dependencies, i)
      queue <- c(conf$dependencies[[i]], queue)
    }
  }
  if (!all(dependencies %in% names(conf$path))) {
    bad <- paste(dependencies[!import %in% names(conf$path)], collapse = ", ")
    stop("The following scripts do not have path: \n  ", bad)
  }
  conf$path[dependencies]
}

attach_lib <- function(name, conf, verbose = TRUE) {
  name <- name[1]
  sub.conf <- conf[[name]]
  import <- sub.conf$import
  if (any(!import %in% names(conf))) {
    stop(
      "The following dependency lib(s) not found in 'conf':\n  ",
      setdiff(import, names(conf)) %>% paste()
    )
  }
  if (name %in% search()) {
    detach_lib(name, conf, verbose)
  }
  for (i in import) {
    attach_lib(i, conf, verbose = verbose)
  }
  for (i in sub.conf$library) {
    if (!paste0("package:", i) %in% search()) {
      .warn('You may also need package [', i, '], try loading it ...')
      suppressMessages(library(i, character.only = TRUE))
    }
  }
  my.env <- new.env()
  for (file in get_source_path(sub.conf)) {
    source(file, chdir = TRUE, local = my.env)
  }
  if (verbose) {
    message("Attach lib '", name, "'")
  }
  # my.env <- asNamespace(my.env)
  attach(my.env, name = name)
  attachNamespace(name)
  invisible(NULL)
}

detach_lib <- function(name, conf, verbose = TRUE) {
  name <- name[1]
  imported.by <- get_imported_by(name, conf)
#  if (any(imported.by %in% search())) {
#    for (i in imported.by) {
#      detach_lib(i, conf, verbose)
#    }
#    return(invisible(NULL))
#  }
  sub.conf <- conf[[name]]
  import <- sub.conf$import
  for (i in import) {
    if (i %in% search()) {
      detach_lib(i, conf)
    }
  }
  if (name %in% search()) {
    if (verbose) {
      message("Detach '", name, "'")
    }
    detach(name, character.only = TRUE)
  }
  invisible(NULL)
}

# Backup Eula ##################################################################
backup_lib <- function(name, conf, outdir, verbose = TRUE) {
  name <- name[1]
  sub.conf <- conf[[name]]
  import <- sub.conf$import

  outdir2 <- file.path(outdir, "Eula/R")
  dir.create(outdir2, FALSE, TRUE)
  outdir2 <- tools::file_path_as_absolute(outdir2)
  dir.create(file.path(outdir2, "lib"), FALSE, TRUE)

  for (i in import) {
    backup_lib(i, conf, outdir, verbose = verbose)
  }
  for (file in get_source_path(sub.conf)) {
    from.file <- tools::file_path_as_absolute(file)
    to.dir <- file.path(outdir2, file)
    if (verbose) {
      message("Backup:\n  from: ", from.file, "\n  to: ", to.dir)
    }
    invisible(file.copy(from.file, to.dir, overwrite = TRUE))
  }
  invisible(NULL)
}

backup_conf <- function(name, conf, outdir, verbose = TRUE) {
  name <- name[1]
  sub.conf <- conf[[name]]
  import <- sub.conf$import

  outdir2 <- file.path(outdir, "Eula/R")
  dir.create(outdir2, FALSE, TRUE)
  outdir2 <- tools::file_path_as_absolute(outdir2)

  imported <- get_depend_tree(name, conf)
  conf.used <- conf[imported]
  yaml::write_yaml(conf.used, file.path(outdir2, "Eula_lib.yaml"))

  for (i in imported) {
    for (file in conf.used[[i]]$src) {
      from.file <- tools::file_path_as_absolute(file)
      to.dir <- dirname(file.path(outdir2, file))
      invisible(file.copy(from.file, to.dir, recursive = TRUE))
    }
  }
  invisible(file.copy(
    tools::file_path_as_absolute(paste0(name, ".R")),
    dirname(file.path(outdir2, paste0(name, ".R"))),
    recursive = TRUE
  ))
  invisible(file.copy(
    tools::file_path_as_absolute("site_lib"),
    outdir2,
    recursive = TRUE
  ))
  invisible(NULL)
}
