#' Plot expression of a gene/feature
#'
#' This function returns a raw ggplot object that you can decorate
#' further
#'
#' @export
mbl_plot_expression <- function(y, gene, group, ...) {
  stopifnot(is(y, "DGEList") || is(y, "EList"))
  assert_string(gene)
  assert_character(group, min.len = 1, max.len = 2)
  gidx <- mbl_fuzzy_find_gene_row(y, gene)

  dat <- mbl_tidy(y[gidx,])
  dat <- .with_aes_columns(dat, group)

  ggplot(dat, aes(.group_by, cpm)) +
    geom_boxplot(outlier.size = 0) +
    geom_jitter(width = 0.25) +
    ylab("log2(cpm)") +
    xlab(paste(group, collapse = "_")) +
    theme(axis.text.x=element_text(angle=90, hjust=1))
}

mbl_heatmap <- function(x, ...) {

}

# Utilify functions ============================================================
# Some serious voodoo is going on here
.with_aes_columns <- function(x, aesthetic, ...) {
  assert_class(x, "data.frame")
  assert_character(aesthetic, min.len = 1, max.len = min(3, ncol(x)))
  assert_subset(aesthetic, colnames(x))

  if (length(aesthetic) == 1L) {
    x[[".group_by"]] <- x[[aesthetic]]
  } else {
    is.cat <- sapply(x[, aesthetic], is.categorical)
    assert_true(all(is.cat))
    x <- tidyr::unite_(x, ".group_by", aesthetic, remove = FALSE)
  }

  x
}

#' Fast and loose way to identify which row a query search for a string across all columns of
mbl_fuzzy_find_gene_row <- function(y, query) {
  stopifnot(is(y, "DGEList") || is(y, "EList"))
  assert_string(query)
  oquery <- query

  gi <- y$genes
  gi$.rn. <- rownames(y)

  query <- tolower(query)

  idx <- NULL
  for (cname in colnames(gi)) {
    vals <- tolower(gi[[cname]])
    i <- which(vals == query)
    if (length(i) > 1) {
      warning("More than one row for `", query, "` found -- taking first one",
              call. = TRUE)
      idx <- i[1L]
      break
    } else if (length(i) == 1) {
      idx <- i
    }
  }
  if (is.null(idx)) {
    stop("Could not find gene using provided query: `", oquery, "`")
  }
  idx
}
