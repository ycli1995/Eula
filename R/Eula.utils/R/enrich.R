
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
  q_gene2term <- gene2term[genes]
  q_gene2term <- q_gene2term[lengths(q_gene2term) > 0]
  if (length(q_gene2term) == 0) {
    message("No gene can be mapped to any term.")
    message("--> return NULL...")
    return(NULL)
  }
  all_genes_in_term <- unique(unlist(term2gene))
  if (is.null(bg.genes)) {
    bg.genes <- all_genes_in_term
  } else {
    if (!is.character(bg.genes)) {
      stop("'bg.genes' must be a character vector or NULL.")
    }
    if (!force.background) {
      bg.genes <- intersect(bg.genes, all_genes_in_term)
    }
  }
  q_gene2term_df <- data.frame(
    gene_id = rep(names(q_gene2term), lengths(q_gene2term)),
    term_id = unlist(q_gene2term)
  )
  q_gene2term_df <- unique(q_gene2term_df)

  q_term2gene <- with(
    q_gene2term_df,
    split(as.character(gene_id), as.character(term_id))
  )
  q_term2gene <- lapply(q_term2gene, intersect, bg.genes)
  q_terms <- names(q_term2gene)

  term2gene <- term2gene[q_terms]
  term2gene <- lapply(term2gene, intersect, bg.genes)
  term2gene <- .select_geneset_by_size(term2gene, min.size, max.size)
  if (length(term2gene) == 0) {
    message("No gene sets have size between ", min.size, " and ", max.size, ".")
    message("--> return NULL...")
    return(NULL)
  }

  q_term2gene <- q_term2gene[names(term2gene)]
  fg_count <- lengths(q_term2gene)
  bg_count <- lengths(term2gene)
  n_bg.genes <- length(bg.genes)
  n_fg_genes <- length(q_gene2term)
  results <- data.frame(
    term_id = names(q_term2gene),
    gene_ratio = if (return.math) {
      fg_count / n_fg_genes
    } else {
      paste(fg_count, "/", n_fg_genes)
    },
    bg_ratio = if (return.math) {
      bg_count / n_bg.genes
    } else {
      paste(bg_count, "/", n_bg.genes)
    },
    fg_count = fg_count,
    bg_count = bg_count,
    genes = unlist(lapply(q_term2gene, paste, collapse = "/")),
    p = phyper(
      q = fg_count - 1,
      m = bg_count,
      n = n_bg.genes - bg_count,
      k = n_fg_genes,
      lower.tail = FALSE
    ),
    stringsAsFactors = FALSE
  )
  results$p_adj <- p.adjust(results$p, method = p.adjust.method)
  results$q <- tryCatch(
    qvalue(p = results$p, lambda = 0.05, pi0.method = "bootstrap")$qvalues,
    error = function(e) NA
  )
  results[order(results$p), ]
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
