
test_that("recodeFactor", {
  x <- c(rep.int("cluster1", 3), rep.int("cluster2", 2), rep.int("C3", 4))
  x <- as.factor(x)

  # one-to-one
  map1 <- list(c2 = "cluster2", c1 = "cluster1", c3 = "cluster3")
  x1 <- recodeFactor(x, map1)
  expect_equal(levels(x1), c("c2", "c1", "C3"))
  x1 <- recodeFactor(x, map1, keep.orders = TRUE)
  expect_equal(levels(x1), c("C3", "c1", "c2"))

  # one-to-many
  map1 <- list(c2 = c("cluster2", "cluster1"), c3 = "C3")
  x1 <- recodeFactor(x, map1)
  expect_equal(levels(x1), c("c2", "c3"))
  x1 <- recodeFactor(x, map1, keep.orders = TRUE)
  expect_equal(levels(x1), c("c3", "c2"))

  # Duplicated old levels
  map1 <- list(c2 = c("cluster2", "C3"), c3 = "C3")
  x1 <- recodeFactor(x, map1)
  expect_equal(levels(x1), c("c2", "cluster1"))
  x1 <- recodeFactor(x, map1, keep.orders = TRUE)
  expect_equal(levels(x1), c("c2", "cluster1"))
})

test_that("norm_fname", {
  expect_equal(norm_fname("A/B"), "A_B")
  expect_equal(norm_fname("A B"), "A_B")
  expect_equal(norm_fname("A//B"), "A__B")
  expect_equal(norm_fname("A  B"), "A__B")
  expect_equal(norm_fname("A(B"), "A_B")
  expect_equal(norm_fname("A(B)"), "A_B_")
  expect_equal(norm_fname("A(B))"), "A_B__")
})

test_that("fetch_default_params", {
  params <- list(features = letters, cells = 1:10)
  defaults <- list(
    features = LETTERS,
    group.by = "group",
    split.by = character(),
    use.name = TRUE,
    assay = NULL,
    reductions = NULL
  )
  params <- fetch_default_params(defaults = defaults, params = params)
  list2env(params, envir = environment())

  expect_null(assay)
  expect_null(reductions)
  expect_true(use.name)
  expect_length(split.by, 0)
  expect_equal(group.by, "group")
  expect_equal(cells, 1:10)
  expect_equal(features, letters)
})
