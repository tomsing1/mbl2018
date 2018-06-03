#' Converts a DGEList into a tidy data.frame
#'
#' This makes downstream exploratory data analysis a bit easier.
#'
#' @export
#' @importFrom edgeR cpm
#' @importFrom reshape2 melt
#' @importFrom broom tidy
#' @method tidy DGEList
tidy.DGEList <- function(x, normalized.lib.sizes = TRUE, prior.count = 3, ...) {
  cpms <- cpm(x, normalized.lib.sizes = normalized.lib.sizes,
              log = TRUE, prior.count = prior.count)
  genes <- x$genes
  genes$.gene_id <- rownames(x)
  gid.col <- sapply(genes, function(xx) all(xx == rownames(x)))
  gid.col <- colnames(genes)[which(gid.col)[1L]]
  if (gid.col != ".gene_id") genes$.gene_id <- NULL

  samples <- x$samples
  samples$.sample_id <- colnames(x)
  sid.col <- sapply(samples, function(xx) all(xx == colnames(x)))
  sid.col <- colnames(samples)[which(sid.col)[1L]]
  if (sid.col != ".sample_id") samples$.sample_id <- NULL

  m <- melt(cpms)
  m <- transform(m, Var1 = as.character(Var1), Var2 = as.character(Var2))
  colnames(m) <- c(gid.col, sid.col, "cpm")
  m$count <- melt(x$counts)[[3L]]

  out <- inner_join(m, genes, by = gid.col)
  out <- inner_join(out, samples, by = sid.col)
  out
}
