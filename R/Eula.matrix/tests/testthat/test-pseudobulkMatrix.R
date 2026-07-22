
test_pseudobulkMatrix <- function(
    mat,
    groups,
    method = c("aggregate", "average"),
    reverse = FALSE,
    bulk = FALSE,
    sparse = FALSE,
    drop = TRUE,
    ...
) {
  method <- match.arg(method)
  groups <- as.factor(groups)
  if (length(unique(groups)) == 1) {
    expect_error(
      out <- pseudobulkMatrix(
        object = mat,
        cell.groups = groups,
        method = method,
        reverse = TRUE,
        bulk = bulk,
        sparse = sparse,
        drop = drop,
        ...
      )
    )
    reverse <- FALSE
  }
  out <- pseudobulkMatrix(
    object = mat,
    cell.groups = groups,
    method = method,
    reverse = reverse,
    bulk = bulk,
    sparse = sparse,
    drop = drop,
    ...
  )
  expect_equal(rownames(out), rownames(mat))
  if (sparse) {
    expect_true(inherits(out, "CsparseMatrix"))
  } else {
    expect_true(is.matrix(out))
  }
  if (bulk) {
    expect_in("bulk", colnames(out))
    a <- out[, "bulk", drop = TRUE]
    b <- switch(
      EXPR = method,
      "average" = Matrix::rowMeans(mat),
      "aggregate" = Matrix::rowSums(mat)
    )
    expect_equal(a, b, ignore_attr = TRUE)
    out <- out[, - which(colnames(out) %in% "bulk"), drop = FALSE]
  }
  if (drop) {
    groups <- droplevels(groups)
  }
  expect_equal(colnames(out), levels(groups))
  for (i in levels(groups)) {
    select <- groups == i
    if (reverse) {
      select <- !select
    }
    a <- out[, i, drop = TRUE]
    b <- switch(
      EXPR = method,
      "average" = Matrix::rowMeans(mat[, select, drop = FALSE]),
      "aggregate" = Matrix::rowSums(mat[, select, drop = FALSE])
    )
    expect_equal(a, b, ignore_attr = TRUE)
  }
  out
}

test_that("`pseudobulkMatrix` for >1 levels.", {
  m <- Matrix::rsparsematrix(10, 20, density = 0.5)
  dimnames(m) <- list(letters[seq_len(nrow(m))], LETTERS[seq_len(ncol(m))])

  labels <- c("AA", "BB", "CC")
  g <- sample(labels, ncol(m), replace = TRUE)
  args <- list(
    method = c("aggregate", "average"),
    reverse = c(TRUE, FALSE),
    bulk = c(TRUE, FALSE),
    sparse = c(TRUE, FALSE),
    drop = c(TRUE, FALSE)
  )
  args <- expand.grid(args, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(args))) {
    out <- test_pseudobulkMatrix(
      mat = m,
      groups = g,
      method = args$method[i],
      reverse = args$reverse[i],
      bulk = args$bulk[i],
      sparse = args$sparse[i],
      drop = args$drop[i]
    )
    out <- test_pseudobulkMatrix(
      mat = as.matrix(m),
      groups = g,
      method = args$method[i],
      reverse = args$reverse[i],
      bulk = args$bulk[i],
      sparse = args$sparse[i],
      drop = args$drop[i]
    )
  }
})
