#' Retrieve the sample annotation
#' 
#' Returns the sample annotations for a particular dataset.
#' 
#' @export
#' @param organism either "mouse", "fly", or "fish".
#' @param dataset either "provided" or "generated"
#' @return a data.frame of sample annotations
get_sample_annotation <- function(organism = c("mouse", "fly", "fish"),
                                  dataset = c("provided", "generated"),
                                  for_kallisto = FALSE,
                                  kallisto_parent_dir = NULL) {
  organism <- match.arg(organism)
  dataset <- match.arg(dataset)
  
  if (dataset == "provided") {
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
    out <- subset(out, species == organism)
  } else {
    stop("We don't have the generated data yet.")
  }
  
  if (for_kallisto) {
    out <- rename(out, sample = "sample_id", condition = "group")
    if (is.null(kallisto_parent_dir)) {
      kallisto_parent_dir <- file.path("/data/may-kallisto", organism)
    }
    out$path <- file.path(kallisto_parent_dir, out$sample)
    out <- select(out, sample, condition, path, source, treatment, genotype)
    rmcols <- sapply(out, function(x) all(is.na(x)))
    out <- out[!rmcols]
    
    bad.files <- !file.exists(out$path)
    if (any(bad.files)) {
      stop("Can't find kallisto directories: ",
           paste(out$path[bad.files]), collapse =",")
    }
  }
  out
}
