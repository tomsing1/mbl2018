#' Loads RNA-seq datasets
#'
#' This function will load an R object that contains all of the expression data
#' for a given organism from the the data that was generated before the
#' workshop, or during the first few days of it.
#'
#' @export
#' @param organism either "mouse", "fly", or "fish"
#' @param dataset either "provided" or "generated"
#' @return a DGEList of the data.
mbl_load_rnaseq <- function(organism = c("mouse", "fly", "fish"),
                            dataset = c("provided", "generated"), ...) {

}
