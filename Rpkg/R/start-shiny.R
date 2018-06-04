#' Export a DGEList as a csv file
#'
#' This file is suitable for:
#' * degust.erc.monash.edu
#' * https://kcvi.shinyapps.io/START/ (https://github.com/jminnier/STARTapp)
#'
#' no pdata
#'
#' @export
#' @param x DGEList
#' @return data.frame
as_csv <- function(x, ...) {
  stopifnot(is(x, "DGEList"))
  gi <- select(x$genes, gene.id = ens_gene, gene.symbol = ext_gene)
  cnts <- x$counts

  new.names <- gsub("_", " ", as.character(x$samples$group))
  new.names <- tools::toTitleCase(new.names)
  new.names <- gsub(" ", "", new.names)
  new.names <- make.unique(new.names, sep = "_")
  new.names <- ifelse(!grepl("_\\d+$", new.names),
                      paste0(new.names, "_0"),
                      new.names)

  colnames(cnts) <- new.names
  out <- cbind(gi, cnts)
  out
}
