library(tximport)
library(edgeR)
library(aws.s3)
aws.signature::use_credentials("mbl18")

dat.dir <- "/data/kallisto"
organism <- "mouse"

build_kallisto_dge <- function(organism, kdir) {
  if (FALSE) {
    organism <- "mouse"
    kdir <- file.path(dat.dir, organism)
  }
  txi <- mbl_get_transcript_annotation(organism)
  si <- mbl_get_sample_annotation(organism, "provided",
                                  kallisto_parent_dir = kdir)
  res <- tximport(si$kallisto_fn, type = "kallisto",
                  tx2gene = txi[, c("target_id", "ens_gene")],
                  ignoreTxVersion = TRUE,
                  countsFromAbundance = "lengthScaledTPM")

  colnames(res$abundance) <- si$sample_id

  gi <- mbl_get_gene_annotation(organism)
  out <- DGEList(counts = res$abundance, samples = select(si, -kallisto_fn))
  out$genes <- gi[rownames(out),]
  out$genes$length <- rowMeans(res$length)
  out
}

if (FALSE) {
  ymouse <- build_kallisto_dge("mouse", file.path(dat.dir, "mouse"))
  s3saveRDS(ymouse, "s3://mbl.data/mapping/may/mouse/kallisto-DGEList.rds")

  yfly <- build_kallisto_dge("fly", file.path(dat.dir, "fly"))
  s3saveRDS(yfly, "s3://mbl.data/mapping/may/fly/kallisto-DGEList.rds")

  yfish <- build_kallisto_dge("fish", file.path(dat.dir, "fish"))
  s3saveRDS(yfish, "s3://mbl.data/mapping/may/fish/kallisto-DGEList.rds")
}

