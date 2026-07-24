
test_that("foldChange works", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  obj <- SeuratObject::pbmc_small
  ans <- Seurat::FoldChange(obj, ident.1 = "0", ident.2 = "1")

  m <- SeuratObject::GetAssayData(obj)
  cells.1 <- colnames(m)[SeuratObject::Idents(obj) == "0"]
  cells.2 <- colnames(m)[SeuratObject::Idents(obj) == "1"]
  out <- foldChange(m, cells.1 = cells.1, cells.2 = cells.2, verbose = FALSE)

  out <- out[rownames(ans), ]
  expect_equal(out[[1]], ans[[1]])
  expect_equal(round(out[["pct.1"]], 3), ans[["pct.1"]])
  expect_equal(round(out[["pct.2"]], 3), ans[["pct.2"]])
})
