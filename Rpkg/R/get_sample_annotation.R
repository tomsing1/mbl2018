#' Retrieve the sample annotation
#' 
#' Returns the sample annotations for a particular dataset.
#' 
#' @export
#' @param organism either "mouse", "fly", or "fish". If `NULL`, annotations
#'   for all experiments are returned.
#' @param dataset either "provided" or "generated"
#' @return a data.frame of sample annotations
get_sample_annotation <- function(organism = NULL,
                                  dataset = c("provided", "generated")) {
  if (!is.null(organism)) {
    organism <- match.arg(organism, c("mouse", "fly", "fish"))
  }
  dataset <- match.arg(dataset)
  
  out <- read.csv(
    paste(
      "https://raw.githubusercontent.com/tomsing1/mbl2018/master",
      "fastq_processing/may/sample_annotation.csv",    
      sep="/"),
    stringsAsFactors = FALSE)
  rownames(out) <- out$sample_id
  
  out$group <- ifelse(out$species == "fly",
                      paste(out$genotype, out$treatment, sep = "_"),
                      paste(out$source, out$genotype, sep = "_"))

  if (!is.null(organism)) {
    out <- subset(out, species == organism)
  }
  
  out
}