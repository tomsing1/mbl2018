# This script copies and renames FASTQ files from a temporary S3 bucket into
# the final destination, e.g. subfolders for each species.

library(readr)
library(aws.s3)
library(here)
library(dplyr)

kSampleAnno <- here::here("fastq_processing", "may", "sample_annotation.csv")
kS3Bucket <- "mbl.data"
kS3InPrefix <- "tmp"
kS3OutPrefix <- "reads/may"

Main <- function() {
  sample_anno <- read_csv(kSampleAnno, col_types = cols())
  stopifnot(!any(duplicated(sample_anno$fastq)))
  purrr::walk2(sample_anno$Key, sample_anno$fastq, .f = function(old, new) {
    message(sprintf("Copying file %s to %s", old, new))
    aws.s3::copy_object(from_object = paste0(kS3InPrefix, "/", old),
                        to_object = paste0(kS3OutPrefix, "/", new),
                        from_bucket = kS3Bucket, to_bucket = kS3Bucket)
  })
}

if (!interactive()) {
  Main()
}
