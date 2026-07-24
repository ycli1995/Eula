seurat_find_markers <- function(
  obj,
  logfc.threshold = 0,
  min.pct = 0,
  min.diff.pct = -Inf,
  verbose = FALSE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 0,
  min.cells.group = 3,
  ...
) {
  Seurat::FindMarkers(
    object = obj,
    logfc.threshold = logfc.threshold,
    min.pct = min.pct,
    min.diff.pct = min.diff.pct,
    verbose = verbose,
    only.pos = only.pos,
    max.cells.per.ident = max.cells.per.ident,
    random.seed = random.seed,
    latent.vars = latent.vars,
    min.cells.feature = min.cells.feature,
    min.cells.group = min.cells.group,
    ...
  )
}

test_de_helper <- function(obj, test.use, slot = "data", ...) {
  ans <- suppressMessages(
    seurat_find_markers(
      obj,
      ident.1 = "0",
      ident.2 = "1",
      test.use = test.use
    )
  )

  m <- SeuratObject::GetAssayData(obj[[Seurat::DefaultAssay(obj)]], slot)
  cells.1 <- colnames(m)[SeuratObject::Idents(obj) == "0"]
  cells.2 <- colnames(m)[SeuratObject::Idents(obj) == "1"]

  out <- suppressMessages(
    testDE(
      object = m,
      cells.1 = cells.1,
      cells.2 = cells.2,
      test.use = test.use
    )
  )
  out <- out[rownames(ans), , drop = FALSE]
  expect_equal(out$p_val, ans$p_val)
  expect_equal(out$p_val_adj, ans$p_val_adj)
}

test_that("differWilcox works", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  obj <- SeuratObject::pbmc_small
  test_de_helper(obj, "wilcox")
})

test_that("differWilcox - limma works", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  obj <- SeuratObject::pbmc_small
  test_de_helper(obj, "wilcox_limma")
})

test_that("differTTest works", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  obj <- SeuratObject::pbmc_small
  test_de_helper(obj, "t")
})

test_that("differROC works", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("ROCR")

  obj <- SeuratObject::pbmc_small
  ans <- seurat_find_markers(
    obj,
    ident.1 = "0",
    ident.2 = "1",
    test.use = "roc"
  )

  m <- SeuratObject::GetAssayData(obj)
  cells.1 <- colnames(m)[SeuratObject::Idents(obj) == "0"]
  cells.2 <- colnames(m)[SeuratObject::Idents(obj) == "1"]

  out <- differROC(m, cells.1, cells.2)[rownames(ans), , drop = FALSE]
  expect_equal(out$AUC, ans$myAUC)
  expect_equal(out$power, ans$power)
})

test_that("differMAST works", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("MAST")

  obj <- SeuratObject::pbmc_small
  test_de_helper(obj, "MAST")
})

test_that("differDESeq2 works", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("DESeq2")

  obj <- SeuratObject::pbmc_small
  m <- SeuratObject::GetAssayData(
    obj[[SeuratObject::DefaultAssay(obj)]],
    "counts"
  )
  obj[[SeuratObject::DefaultAssay(obj)]] <- SeuratObject::SetAssayData(
    obj[[SeuratObject::DefaultAssay(obj)]],
    "counts",
    as.matrix(m) + 1
  )
  test_de_helper(obj, "DESeq2", slot = "counts")
})
