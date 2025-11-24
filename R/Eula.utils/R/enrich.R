
#' @importFrom qvalue qvalue
#' @export
runEnrich <- function(
    genes,
    term2gene,
    gene2term,
    bg.genes = NULL,
    p.adjust.method = "BH",
    min.size = 10,
    max.size = 500,
    force.background = FALSE,
    return.math = TRUE
) {
  if (!is.character(genes)) {
    stop("'genes' must be a character vector.")
  }
  if (length(genes) == 0) {
    stop("Number of input gene is 0.")
  }
  if (!is.list(term2gene)) {
    stop("'term2gene' must be a list: term id -> gene id")
  }
  if (!is.list(gene2term)) {
    stop("'gene2term' must be a list: gene id -> term id")
  }
  genes <- unique(genes)
  q.gene2term <- gene2term[genes]
  q.gene2term <- q.gene2term[lengths(q.gene2term) > 0]
  if (length(q.gene2term) == 0) {
    message("No gene can be mapped to any term.")
    message("--> return NULL...")
    return(NULL)
  }
  all.genes.in.term <- unique(unlist(term2gene))
  if (is.null(bg.genes)) {
    bg.genes <- all.genes.in.term
  } else {
    if (!is.character(bg.genes)) {
      stop("'bg.genes' must be a character vector or NULL.")
    }
    if (!force.background) {
      bg.genes <- intersect(bg.genes, all.genes.in.term)
    }
  }
  q.gene2term_df <- data.frame(
    gene.id = rep(names(q.gene2term), lengths(q.gene2term)),
    term.id = unlist(q.gene2term)
  )
  q.gene2term_df <- unique(q.gene2term_df)

  q.term2gene <- with(
    q.gene2term_df,
    split(as.character(gene.id), as.character(term.id))
  )
  q.term2gene <- lapply(q.term2gene, intersect, bg.genes)
  q.terms <- names(q.term2gene)

  term2gene <- term2gene[q.terms]
  term2gene <- lapply(term2gene, intersect, bg.genes)
  term2gene <- .select_geneset_by_size(term2gene, min.size, max.size)
  if (length(term2gene) == 0) {
    message("No gene sets have size between ", min.size, " and ", max.size, ".")
    message("--> return NULL...")
    return(NULL)
  }

  q.term2gene <- q.term2gene[names(term2gene)]
  fg.count <- lengths(q.term2gene)
  bg.count <- lengths(term2gene)
  n.bg.genes <- length(bg.genes)
  n.fg.genes <- length(q.gene2term)
  results <- data.frame(
    term_id = names(q.term2gene),
    gene_ratio = if (return.math) {
      fg.count / n.fg.genes
    } else {
      paste(fg.count, "/", n.fg.genes)
    },
    bg_ratio = if (return.math) {
      bg.count / n.bg.genes
    } else {
      paste(bg.count, "/", n.bg.genes)
    },
    fg_count = fg.count,
    bg_count = bg.count,
    genes = unlist(lapply(q.term2gene, paste, collapse = "/")),
    p_val = phyper(
      q = fg.count - 1,
      m = bg.count,
      n = n.bg.genes - bg.count,
      k = n.fg.genes,
      lower.tail = FALSE
    ),
    stringsAsFactors = FALSE
  )
  results$p_val_adj <- p.adjust(results$p_val, method = p.adjust.method)
  results$q_val <- tryCatch(
    qvalue(p = results$p_val, lambda = 0.05, pi0.method = "bootstrap")$qvalues,
    error = function(e) NA
  )
  results[order(results$p_val), ]
}

.select_geneset_by_size <- function(genesets, min.size, max.size) {
  if (is.na(min.size) || is.null(min.size)) {
    min.size <- 1
  }
  if (is.na(max.size) || is.null(max.size)) {
    max.size <- Inf
  }
  size <- lengths(genesets)
  genesets[size >= min.size & size <= max.size]
}
