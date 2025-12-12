
test_that("normalizeFileName", {
  expect_equal(normalizeFileName("A/B"), "A_B")
  expect_equal(normalizeFileName("A B"), "A_B")
  expect_equal(normalizeFileName("A//B"), "A__B")
  expect_equal(normalizeFileName("A  B"), "A__B")
  expect_equal(normalizeFileName("A(B"), "A_B")
  expect_equal(normalizeFileName("A(B)"), "A_B_")
  expect_equal(normalizeFileName("A(B))"), "A_B__")
})

test_that("splitArgs", {
  args <- list(
    arg1 = "A,B,C",
    arg2 = "AA",
    arg3 = NULL,
    arg4 = logical()
  )
  new.args <- lapply(args, splitArgs)
  nochange <- setdiff(names(args), "arg1")
  for (i in nochange) {
    expect_equal(new.args[[i]], args[[i]])
  }
  expect_equal(new.args[["arg1"]], c("A", "B", "C"))
})

test_that("getArgList", {
  genes <- c("gene1", "gene2", "gene3")
  genes_file <- withr::local_tempfile()
  write(genes, genes_file)

  args <- list(
    empty = NULL,
    empty.char = character(),
    empty.logi = logical(),
    empty.num = numeric(),
    coma.char = "AA,B,CC",
    multi.char = letters[1:5],
    multi.logi = rep(TRUE, 3),
    multi.num = 1:3,
    list.args = list(AA = 1, BB = 2),
    gene.list = genes_file
  )
  new.args <- lapply(args, getArgList)
  nochange <- setdiff(names(args), c("coma.char", "gene.list"))
  for (i in nochange) {
    expect_equal(new.args[[i]], args[[i]])
  }
  expect_equal(new.args[["coma.char"]], c("AA", "B", "CC"), ignore_attr = TRUE)
  expect_equal(new.args[["gene.list"]], genes, ignore_attr = TRUE)
})

test_that("getDefaultArgs", {
  params <- list(features = letters, cells = 1:10)
  defaults <- list(
    features = LETTERS,
    group.by = "group",
    split.by = character(),
    use.name = TRUE,
    assay = NULL,
    reductions = NULL
  )
  params <- getDefaultArgs(params, defaults)
  list2env(params, envir = environment())

  expect_null(assay)
  expect_null(reductions)
  expect_true(use.name)
  expect_length(split.by, 0)
  expect_equal(group.by, "group")
  expect_equal(cells, 1:10)
  expect_equal(features, letters)
})
