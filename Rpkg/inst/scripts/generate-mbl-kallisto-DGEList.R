library(tximport)
library(edgeR)
library(aws.s3)
library(rhdf5)
library(dplyr)

aws.signature::use_credentials("mbl18")

dat.dir <- "/data/kallisto"
organism <- "mouse"

sample.info <- read.csv("https://s3.amazonaws.com/mbl.data/reads/mbl/sample_annotations.csv",
                        stringsAsFactors = FALSE)
sample.info$kallisto_fn <- file.path(
  dat.dir, sample.info$species,
  sample.info$sample_id,
  "abundance.h5")
all(file.exists(sample.info$kallisto_fn))


build_kallisto_dge <- function(organism, kdir) {
  txi <- mbl_get_transcript_annotation(organism)

  if (organism == "mouse") {
    # may mouse: sample_id, source, geneotype, animal_id
    si <- sample.info %>%
      subset(species == organism) %>%
      transform(group = paste(source, genotype, treatment, sep = "_"),
                animal_id = paste(genotype, replicate, sep = "_"),
                stringsAsFactors = FALSE) %>%
      transform(group = sub("_NA", "", group), stringsAsFactors = FALSE) %>%
      select(group, sample_id, source, genotype, treatment, animal_id, kallisto_fn)
  } else if (organism == "fish") {
    # may fish: group, sample_id, source, genotype, animal_id
    si <- sample.info %>%
      subset(species == organism) %>%
      transform(group = paste(source, genotype, sep = "_"),
                treatment = NA_character_,
                animal_id = paste(genotype, replicate, sep = "_"),
                stringsAsFactors = FALSE) %>%
      select(group, sample_id, source, genotype, animal_id, kallisto_fn)
  } else {
    # may fly: group, sample_id, source, genotype, treatment, animal_id
    si <- sample.info %>%
      subset(species == organism) %>%
      transform(group = paste(source, genotype, sep = "_"),
                treatment = NA_character_,
                animal_id = NA_character_,
                stringsAsFactors = FALSE) %>%
      select(group, sample_id, source, genotype, animal_id, kallisto_fn)
  }

  rownames(si) <- si$sample_id

  res <- tximport(si$kallisto_fn, type = "kallisto",
                  tx2gene = txi[, c("target_id", "ens_gene")],
                  ignoreTxVersion = TRUE,
                  countsFromAbundance = "lengthScaledTPM")

  colnames(res$counts) <- si$sample_id

  gi <- mbl_get_gene_annotation(organism)
  out <- DGEList(counts = res$counts, samples = select(si, -kallisto_fn))
  out$samples$group <- as.character(out$samples$group)
  gnz <- gi[rownames(out),]
  if (!all.equal(gnz$ens_gene, rownames(out))) {
    browser()
  }

  out$genes <- gnz <- gi[rownames(out),]
  out$genes$length <- rowMeans(res$length[rownames(out),])
  out
}

if (FALSE) {
  ymouse <- build_kallisto_dge("mouse", file.path(dat.dir, "mouse"))
  s3saveRDS(ymouse, "s3://mbl.data/mapping/mbl/mouse/kallisto-DGEList.rds")

  yfly <- build_kallisto_dge("fly", file.path(dat.dir, "fly"))
  s3saveRDS(yfly, "s3://mbl.data/mapping/mbl/fly/kallisto-DGEList.rds")

  yfish <- build_kallisto_dge("fish", file.path(dat.dir, "fish"))
  s3saveRDS(yfish, "s3://mbl.data/mapping/mbl/fish/kallisto-DGEList.rds")
}

assemble_dataset <- function(organism) {
  may <- mbl_load_rnaseq(organism, "may")
  mbl <- mbl_load_rnaseq(organism, "mbl")

  if (!all.equal(rownames(may$samples), colnames(may$counts))) {
    stop("may data not concordant")
  }
  if (!all.equal(rownames(mbl$samples), colnames(mbl$counts))) {
    stop("mbl data not concordant")
  }
  may$samples$dataset <- "may"
  mbl$samples$dataset <- "mbl"

  # stopifnot(all.equal(colnames(may$samples), colnames(mbl$samples)))
  stopifnot(all.equal(rownames(may), rownames(mbl)))

  cnts <- cbind(may$counts, mbl$counts)
  smpls <- bind_rows(may$samples, mbl$samples) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    select(-lib.size, -norm.factors)

  # there are duplicate sample identifiers
  rn <- paste(smpls$dataset, smpls$sample_id, sep="_")
  rownames(smpls) <- rn
  colnames(cnts) <- rn

  out <- DGEList(cnts, samples = smpls, genes = mbl$genes)
}

# Assemble all datasets
if (FALSE) {
  all.mouse <- assemble_dataset("mouse")
  s3saveRDS(all.mouse, "s3://mbl.data/mapping/all/mouse/kallisto-DGEList.rds")

  all.fly <- assemble_dataset("fly")
  s3saveRDS(all.fly, "s3://mbl.data/mapping/all/fly/kallisto-DGEList.rds")

  all.fish <- assemble_dataset("fish")
  s3saveRDS(all.fish, "s3://mbl.data/mapping/all/fish/kallisto-DGEList.rds")

}

