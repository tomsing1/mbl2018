#' Perform a principle components analysis over a DGEList
#'
#' This function performs a PCA over a DGEList object, and returns the
#' loadings (the position of the samples on the new axes/PCs), and the perecent
#' of variance explained by each PC.
#'
#' Refer to XXX section of the tutorial for more information.
#' FIXME: Provide link to PCA section of tutorial
#'
#' @export
#' @rdname mbl_pca
#'
#' @param x a `DGEList` of expression data
#' @param ntop the number of genes to include in the PCA decomposition. Genes
#'   are ordered by ones with the highest variance to lowest.
#' @param center a logical value indicating whether the variables should be
#'   shifted to be zero centered. Alternately, a vector of length equal the
#'   number of columns of x can be supplied. The value is passed to scale.
#' @param scale. a logical value indicating whether the variables should be
#'   shifted to be zero centered. Alternately, a vector of length equal the
#'   number of columns of x can be supplied. The value is passed to scale.
#' @return a list with `$loadings`, `$percentVar`, ...
mbl_pca <- function(x, ntop = 500, center = TRUE, scale. = FALSE, ...,
                    pcs = 1:10) {
  UseMethod("mbl_pca", x)
}

#' @export
#' @method mbl_pca DGEList
mbl_pca.DGEList <- function(x, ntop = 500, center = TRUE, scale. = FALSE, ...,
                            pcs = 1:10, prior.count = 3, log = TRUE) {
  assert_class(x, "DGEList")
  m <- cpm(x, prior.count = 3, log = TRUE, normalized.lib.sizes = TRUE)
  mbl_pca(m, ntop = ntop, center = center, scale. = scale.,
          row_covariates = x$genes, col_covariates = x$samples, ..., pcs = pcs)
}

#' @export
#' @method mbl_pca EList
mbl_pca.EList <- function(x, ntop = 500, center = TRUE, scale. = FALSE, ...,
                          pcs = 1:10) {
  assert_class(x, "EList")
  mbl_pca(x$E, ntop = ntop, center = center, scale. = scale.,
          row_covariates = x$genes, col_covariates = x$targets, ..., pcs = pcs)
}

#' @method mbl_pca matrix
mbl_pca.matrix <- function(x, ntop = 500, center = TRUE, scale. = FALSE, ...,
                           row_covariates = NULL, col_covariates = NULL,
                           pcs = 1:10) {
  assert_class(x, "matrix")
  assert_integerish(ntop, lower = 2, upper = nrow(x))
  assert_integerish(pcs)
  bad.pcs <- pcs < 1 | pcs > ncol(x) | pcs > nrow(x)
  if (any(bad.pcs)) {
    # warning()
    pcs <- pcs[!bad.pcs]
  }

  if (is(row_covariates, "data.frame")) {
    assert_true(nrow(x) == nrow(row_covariates))
    assert_character(rownames(row_covariates))
    assert_true(all(rownames(x) == rownames(row_covariates)))
  }
  if (is(col_covariates, "data.frame")) {
    assert_true(ncol(x) == nrow(col_covariates))
    assert_character(rownames(col_covariates))
    assert_true(all(colnames(x) == rownames(col_covariates)))
  }

  rv <- matrixStats::rowVars(x)
  take <- head(order(rv, decreasing = TRUE), ntop)

  pca <- prcomp(t(x[take,]))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  names(percentVar) <- paste0("PC", seq(percentVar))

  dat <- as.data.frame(pca$x)

  percentVar <- percentVar[colnames(dat)]

  if (is(col_covariates, "data.frame")) {
    dat <- cbind(dat, col_covariates[rownames(dat),,drop = FALSE])
  }

  result <- list(data = dat, percentVar = percentVar, prcomp = pca)
  class(result) <- c("mblpca")
  result
}

#' @noRd
#' @export
#' @method print mblpca
print.mblpca <- function(x, ...) {
  # dx <- deparse(x)
  # browser()
  # mc <- match.call()
  cat(format(x, ...), "\n")
}

format.mblpca <- function(x, ...) {
  # browser()
  # xs <- substitute(x)
  out <- paste(
    "===========================================================\n",
    "A simplified PCA result object from mbl2018 package\n",
    "-----------------------------------------------------------\n",
    "\n",
    "  * Rotated data (samples) with sample covariates are found\n",
    "    in: `result$data`\n",
    "  * Percent variance explained per PC found in:\n",
    "    `result$percentVar`\n",
    "  * The original result from the 'real R' PCA call is found\n",
    "    in: `result$prcomp`\n",
    "\n",
    " Type `?mbl_pca` for more help\n",
    "===========================================================\n",
    sep = "")
  out
}

