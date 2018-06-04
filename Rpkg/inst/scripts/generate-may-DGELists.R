# install.packages("aws.s3")
library(edgeR)
library(aws.s3)
library(dplyr)
library(readr)
library(tidyr)
aws.signature::use_credentials()

# Lots of assumptions in this script
# 1. You have run:
#    $ aws s3 sync s3://mbl.data/star-alignments/may /data/alignments
stopifnot(dir.exists('/data/alignments'))

align.dir <- "/data/alignments/may"

genDGEList <- function(organism = c("mouse", "fly", "fish")) {
  count.base <- file.path("/data/alignments", organism)
  sample.names <- dir(count.base, full.names = FALSE)
  
  # create count table and gather alignment info per sample ==============
  count.fn <- file.path(count.base, sample.names, "ReadsPerGene.out.tab")
  counts.all <- lapply(count.fn, function(fn) {
    cnts <- read_tsv(fn, col_names = FALSE, skip = 4, col_types = "ciii")
    cnts <- data.frame(count = cnts[[2]], row.names = cnts[[1]])
    colnames(cnts) <- basename(dirname(fn))
    header <- read_tsv(fn, col_names = FALSE, n_max = 4,  col_types = "ciii")
    meta <- as.data.frame(t(header[, 2]))
    colnames(meta) <- header[[1]]
    list(counts = cnts, meta = meta)
  })
  counts <- do.call(cbind, lapply(counts.all, "[[", "counts"))
  rownames(counts) <- rownames(counts.all[[1]]$counts)
  
  ameta <- bind_rows(lapply(counts.all, "[[", "meta"))
  rownames(ameta) <- colnames(counts)
  
  si <- cbind(
    get_sample_annotation(organism, "provided")[rownames(ameta),],
    ameta)
  
  if (any(is.na(si$sample_id))) stop("sample information matching went south")
  
  gene.info <- get_gene_annotation(organism)
  gxref <- match(rownames(counts), gene.info$ensembl_gene_id)
  if (any(is.na(gxref))) stop("Gene matching went south")
  
  gi <- gene.info[gxref, c("ensembl_gene_id", "symbol", "biotype")]
  stopifnot(all.equal(rownames(counts), rownames(gi)))

  out <- DGEList(counts, samples = si, genes = gi)
  # organism specific group specification
  if (organism == "fly") {
    out$samples$group <- paste(si$genotype, si$treatment, sep = "_")
  } else {
    out$samples$group <- paste(si$source, si$genotype, sep = "_")
  }
  
  calcNormFactors(out)
}

if (FALSE) {
  y.mouse <- genDGEList("mouse")
  s3saveRDS(y.mouse, "s3://mbl.data/star-alignments/may/mouse/DGEList.rds")
  
  y.fly <- genDGEList("fly")
  s3saveRDS(y.fly, "s3://mbl.data/star-alignments/may/fly/DGEList.rds")
  
  y.fish <- genDGEList("fish")
  s3saveRDS(y.fish, "s3://mbl.data/star-alignments/may/fish/DGEList.rds")
}

