
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
  args <- list(arg1 = "A,B,C", arg2 = "AA", arg3 = NULL, arg4 = logical())

  args <- lapply(args, splitArgs)

  expect_equal(args[["arg1"]], c("A", "B", "C"))
  expect_equal(args[["arg2"]], "AA")
  expect_null(args[["arg3"]])
  expect_type(args[["arg4"]], "logical")
  expect_length(args[["arg4"]], 0)
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
  params <- getDefaultArgs(defaults, params)
  list2env(params, envir = environment())

  expect_null(assay)
  expect_null(reductions)
  expect_true(use.name)
  expect_length(split.by, 0)
  expect_equal(group.by, "group")
  expect_equal(cells, 1:10)
  expect_equal(features, letters)
})
