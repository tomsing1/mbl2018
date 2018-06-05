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
  si <- mbl_get_sample_annotation(organism, "may",
                                  kallisto_parent_dir = kdir)
  res <- tximport(si$kallisto_fn, type = "kallisto",
                  tx2gene = txi[, c("target_id", "ens_gene")],
                  ignoreTxVersion = TRUE,
                  countsFromAbundance = "lengthScaledTPM")

  colnames(res$counts) <- si$sample_id

  gi <- mbl_get_gene_annotation(organism)
  out <- DGEList(counts = res$counts, samples = select(si, -kallisto_fn))
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

if (FALSE) {
  # Saving sample information for organisms for students to use later
  si.mouse <- mbl_get_sample_annotation("mouse")
  s3write_using(si.mouse, write.csv,
                object = "s3://mbl.data/mapping/may/mouse/sample-information.csv")

  si.fly <- mbl_get_sample_annotation("fly")
  s3write_using(si.fly, write.csv,
                object = "s3://mbl.data/mapping/may/fly/sample-information.csv")

  si.fish <- mbl_get_sample_annotation("fish")
  s3write_using(si.fish, write.csv,
                object = "s3://mbl.data/mapping/may/fish/sample-information.csv")
}

if (FALSE) {
  sim <- s3read_using(read.csv, row.names = 1,
                      object = "s3://mbl.data/mapping/may/mouse/sample-information.csv")
}
