library(edgeR)
library(aws.s3)
library(dplyr)
library(readr)
aws.signature::use_credentials()

# Lots of assumptions in this script
# 1. You have run:
#    $ aws s3 sync s3://mbl.data/star-alignments/may /data/alignments

align.dir <- "/data/alignments/may"
anno.base <- "references/%s/gene_table.csv"

genDGEList <- function(organism = c("mouse", "fly", "fish")) {
  gene.info <- s3read_using(read_csv, object = sprintf(anno.base, organism),
                            bucket = "mbl.data")
}
