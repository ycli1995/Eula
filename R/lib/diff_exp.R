
suppressMessages(library(future))
suppressMessages(library(future.apply))

.Wilcox <- function(data.use, cells.1, cells.2, limma = FALSE, ...) {
  data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
  my.sapply <- ifelse(nbrOfWorkers() == 1, sapply, future_sapply)

  presto.check <- requireNamespace("presto", quietly = TRUE)
  limma.check <- requireNamespace("limma", quietly = TRUE)

  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(group.info[, "group"])

  if (presto.check[1] && (!limma)) {
    data.use <- data.use[, rownames(group.info), drop = FALSE]
    res <- presto::wilcoxauc(X = data.use, y = group.info[, "group"])
    res <- res[1:(nrow(res)/2), ]
    return(data.frame(p_val = res$pval, row.names = rownames(data.use)))
  }

  if (!limma) {
    p_val <- my.sapply(seq_len(nrow(data.use)), function(i) {
      wilcox.test(data.use[i, ] ~ group.info[, "group"], ...)$p.value
    })
    return(data.frame(p_val, row.names = rownames(data.use)))
  }

  if (!limma.check) {
    stop(
      "To use the limma implementation of the Wilcoxon Rank Sum Test,
       please install the limma package:
       --------------------------------------------
       install.packages('BiocManager')
       BiocManager::install('limma')
       --------------------------------------------"
    )
  }
  overflow <- is.na(suppressWarnings(ncol(data.use) * ncol(data.use)))
  if (overflow) {
    stop("Use 'limma' with ", ncol(data.use), " cells may cause overflow.")
  }
  j <- seq_along(cells.1)
  p_val <- my.sapply(seq_len(nrow(data.use)), function(i) {
    min(1, 2 * min(limma::rankSumTestWithCorrelation(j, data.use[i, ])))
  })
  return(data.frame(p_val, row.names = rownames(data.use)))
}

.TTest <- function(data.use, cells.1, cells.2, ...) {
  my.sapply <- ifelse(nbrOfWorkers() == 1, sapply, future_sapply)
  p_val <- unlist(my.sapply(
    seq_len(nrow(data.use)),
    function(i) t.test(data.use[i, cells.1], data.use[i, cells.2], ...)$p.value
  ))
  data.frame(p_val, row.names = rownames(data.use))
}

.MAST <- function(data.use, cells.1, cells.2, latent.vars = NULL, ...) {
  if (!requireNamespace("MAST", quietly = TRUE)) {
    stop("Please install MAST - learn more at https://github.com/RGLab/MAST")
  }

  group.info <- data.frame(row.names = c(cells.1, cells.2))
  latent.vars <- latent.vars %||% group.info
  group.info[cells.1, "condition"] <- "Group1"
  group.info[cells.2, "condition"] <- "Group2"
  group.info[, "condition"] <- factor(group.info[, "condition"])

  latent.vars.names <- c("condition", colnames(latent.vars))
  latent.vars <- cbind(latent.vars, group.info)
  latent.vars$wellKey <- rownames(latent.vars)
  fdat <- data.frame(primerid = rownames(data.use))
  rownames(fdat) <- fdat[, 1]
  
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(data.use),
    check_sanity = FALSE,
    cData = latent.vars,
    fData = fdat
  )
  SummarizedExperiment::colData(sca)$condition <- relevel(
    factor(SummarizedExperiment::colData(sca)$condition),
    ref = "Group1"
  )
  fmla <- as.formula(paste0(" ~ ", paste(latent.vars.names, collapse = "+")))
  
  zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
  summaryDt <- MAST::summary(zlmCond, doLRT = 'conditionGroup2')$datatable
  summaryDt <- summaryDt[summaryDt[, "component"] == "H", , drop = FALSE]
  data.frame(p_val = summaryDt[, 4], row.names = summaryDt[, 1])
}

