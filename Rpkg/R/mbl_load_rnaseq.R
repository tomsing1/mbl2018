#' Loads RNA-seq dataset generated for MBL2018 Neurobiology Course.
#'
#' @description
#' This loads curated DGEList objects for the RNA-seq datasets generated here
#' for the specified organism.
#'
#' By default this will return a DGEList that has both the data generated
#' prior to the course (the "may" dataset), and the data you generated here
#' into one object. There is a `"dataset"` column in the `$samples` data.frame
#' that can be used to split the two datasets.
#'
#' Alternatively you can load on specific dataset or another by specifing
#' the `dataset` parameter to the function.
#'
#' @export
#' @param organism either "mouse", "fly", or "fish"
#' @param dataset either `"all"`, `"mbl"`, or `"may"`.
#' @return a DGEList of the RNAseq data.
mbl_load_rnaseq <- function(organism = c("mouse", "fly", "fish"),
                            dataset = c("all", "may", "mbl"), ...) {
  organism <- match.arg(organism)
  dataset <- match.arg(dataset)

  url <- "https://s3.amazonaws.com/mbl.data/mapping/%s/%s/kallisto-DGEList.rds"
  url <- sprintf(url, dataset, organism)
  readRDS(url(url))
}
