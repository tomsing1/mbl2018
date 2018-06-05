#' Retrieve the sample annotation
#'
#' Returns the sample annotations for a particular dataset.
#'
#' @export
#' @param organism either "mouse", "fly", or "fish".
#' @param dataset either "may" or "generated"
#' @return a data.frame of sample annotations
mbl_get_sample_annotation <- function(organism = c("mouse", "fly", "fish"),
                                      dataset = c("may", "generated"),
                                      kallisto_parent_dir = NULL,
                                      trim = TRUE) {
  organism <- match.arg(organism)
  dataset <- match.arg(dataset)

  if (dataset == "may") {
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

  if (organism %in% c("mouse", "fly")) {
    out$animal_id <- paste(out$genotype, out$replicate, sep = "_")
  } else {
    out$animal_id <- paste("fish", out$replicate, sep = "_")
  }

  out <- select(out, group, sample_id, source, genotype, treatment, animal_id,
                everything())
  if (trim) {
    out <- select(out, group:animal_id)
  }

  if (is.character(kallisto_parent_dir)) {
    out$kallisto_fn = file.path(kallisto_parent_dir, out$sample_id,
                                "abundance.h5")
  }

  # remove columns that are all NA
  rm.me <- sapply(out, function(vals) all(is.na(vals)))
  out <- out[, !rm.me, drop = FALSE]
  out
}
