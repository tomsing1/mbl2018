#' Plot expression of a gene/feature
#'
#' This function returns a raw ggplot object that you can decorate
#' further
#'
#' @export
#' @param y the DGEList or voomed object your data is in
#' @param gene the name of the gene identifier or symbol you want to plot
#' @param group what column in the $samples (or $targets, for vm) pheno table
#'   to use to plot split expression across the x axis
#' @param color_by name of column (like `group` param) to color your points by
#' @return a ggplot object
mbl_plot_expression <- function(y, gene, group, color_by = NULL, ...) {
  stopifnot(is(y, "DGEList") || is(y, "EList"))
  assert_string(gene)
  assert_character(group, min.len = 1, max.len = 2)
  gidx <- mbl_fuzzy_find_gene_row(y, gene)

  dat <- mbl_tidy(y[gidx,])

  # Add column to data.frame to indicate group of observations
  dat <- .with_aes_columns(dat, group, ".group_by")

  if (!is.null(color_by)) {
    assert_character(color_by, min.len = 1, max.len = 2)
    dat <- .with_aes_columns(dat, color_by, ".color_by")
  }

  gg <- ggplot(dat, aes(.group_by, cpm)) +
    geom_boxplot(outlier.size = 0) +
    ylab("log2(cpm)") +
    xlab(paste(group, collapse = "_")) +
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    ggtitle(paste(gene, "Expression"))

  if (is.character(color_by))   {
    gg <- gg+ geom_jitter(aes(color = .color_by), width = 0.25)
  } else {
    gg <- gg + geom_jitter(width = 0.25)
  }

  gg
}

mbl_heatmap <- function(x, ...) {

}

# Utilify functions ============================================================
# Some serious voodoo is going on here
.with_aes_columns <- function(x, aesthetic, out_column, ...) {
  assert_class(x, "data.frame")
  assert_character(aesthetic, min.len = 1, max.len = min(3, ncol(x)))
  assert_subset(aesthetic, colnames(x))

  if (length(aesthetic) == 1L) {
    x[[out_column]] <- x[[aesthetic]]
  } else {
    is.cat <- sapply(x[, aesthetic], is.categorical)
    assert_true(all(is.cat))
    x <- tidyr::unite_(x, out_column, aesthetic, remove = FALSE)
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
