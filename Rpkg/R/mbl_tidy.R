#' Converts your expression data into a (huge) "tidy" data.frame
#'
#' Use this function to extract all of your expression data into a data.frame
#' that can be used with dplyr and ggplot to make arbitrariy informative
#' plots.
#'
#' @export
#' @seealso [mbl_plot_expression()]
#'
#' @param x The expression object to tidy
#' @return a (huge) data.frame with your expression data. Each row holds the
#'   expression of one gene in one sample. The columns include all of the
#'   gene- and sample-level metadata for the obseration.
#'
#' @examples
#' # Make a boxplot with points of Fxyd6 in the cheek wildtype/knockout
#' y <- mbl_load_rnaseq("mouse", dataset = "mbl")
#' ydat <- mbl_tidy(y) # all of the rnaseq data
#' gdat <- filter(ydat, source == "cheek", symbol == "Fxyd6")
#' ggplot(gdat, aes(x = genotype, y = cpm)) +
#'   geom_boxplot() +
#'   geom_point()
mbl_tidy <- function(x, ...) {
  UseMethod("mbl_tidy", x)
}

#' @rdname mbl_tidy
#' @export
#' @importFrom edgeR cpm
#' @importFrom reshape2 melt
#' @method mbl_tidy DGEList
mbl_tidy.DGEList <- function(x, normalized.lib.sizes = TRUE, prior.count = 3, ...) {
  mats <- list(
    cpm = cpm(x, normalized.lib.sizes = normalized.lib.sizes,
              log = TRUE, prior.count = prior.count),
    count = x$counts)

  mbl_tidy.core(mats, genes = x$genes, samples = x$samples)
}

#' @rdname mbl_tidy
#' @export
#' @method mbl_tidy EList
mbl_tidy.EList <- function(x, ...)  {
  mats <- list(cpm = x$E)
  if (is.matrix(x$weights)) {
    mats$weight <- x$weights
    rownames(mats$weight) <- rownames(x)
    colnames(mats$weight) <- colnames(x)
  } else {
    names(mats)[1L] <- "value"
  }

  mbl_tidy.core(mats, genes = x$genes, samples = x$targets)
}

mbl_tidy.core <- function(mats, genes, samples, ...) {
  if (is.matrix(mats)) mats <- list(value = mats)
  stopifnot(is.list(mats))
  stopifnot(all(sapply(mats, is.matrix)))
  assert_named(mats, type = "unique")

  rnames <- rownames(mats[[1]])
  snames <- colnames(mats[[1]])
  genes$.gene_id <- rnames
  gid.col <- sapply(genes, function(xx) all(xx == rnames))
  gid.col <- colnames(genes)[which(gid.col)[1L]]
  if (gid.col != ".gene_id") genes$.gene_id <- NULL

  samples$.sample_id <- snames
  sid.col <- sapply(samples, function(xx) all(xx == snames))
  sid.col <- colnames(samples)[which(sid.col)[1L]]
  if (sid.col != ".sample_id") samples$.sample_id <- NULL

  adat.all <- lapply(names(mats), function(mname) {
    m <- mats[[mname]]
    stopifnot(all.equal(rownames(m), rnames))
    m <- melt(m)
    m <- transform(m, Var1 = as.character(Var1), Var2 = as.character(Var2))
    colnames(m) <- c(gid.col, sid.col, mname)
    m
  })
  adat <- do.call(cbind, adat.all)
  # if there were multiple matrices, there will be multiple sample_id columns
  # so we remove those
  adat <- adat[, !duplicated(colnames(adat))]
  out <- inner_join(adat, genes, by = gid.col)
  out <- inner_join(out, samples, by = sid.col)
  out
}
