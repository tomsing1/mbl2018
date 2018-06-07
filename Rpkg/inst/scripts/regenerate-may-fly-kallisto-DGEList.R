library(dplyr)
library(tximport)
library(edgeR)
library(aws.s3)
aws.signature::use_credentials("mbl18")

organism <- "fly"

sample.info <- read.csv("https://s3.amazonaws.com/mbl.data/reads/may/sample_annotations.csv",
                        stringsAsFactors = FALSE)

si <- subset(sample.info, species == organism) %>%
  transform(genotype = sub("wt$", "wildtype", genotype),
            stringsAsFactors = FALSE,
            animal_id = paste(genotype, replicate, sep = "_"),
            treatment = gsub(" ", "_", treatment)) %>%
  # transform(source = gsub(" ", "_", source),
  #           genotype = gsub(" ", "_", genotype),
  #           treatment = gsub(" ", "_", treatment)) %>%
  transform(genotype = gsub(" ", "_", genotype),
            animal_id = gsub(" ", "_", animal_id),
            group = paste(source, genotype, treatment, sep = "_"),
            stringsAsFactors = FALSE) %>%
  transform(group = gsub(" ", "_", group)) %>%
  select(group, sample_id, source, genotype, treatment, animal_id)

si$kallisto_fn <- with(si, paste0("/data/mayfly/", sample_id, "/abundance.h5"))
all(file.exists(si$kallisto_fn))

txi <- mbl_get_transcript_annotation(organism)
res <- tximport(si$kallisto_fn, type = "kallisto",
                tx2gene = txi[, c("target_id", "ens_gene")],
                ignoreTxVersion = TRUE,
                countsFromAbundance = "lengthScaledTPM")
colnames(res$counts) <- si$sample_id
gi <- mbl_get_gene_annotation(organism)
out <- DGEList(counts = res$counts, samples = select(si, -kallisto_fn))
gnz <- gi[rownames(out),]
stopifnot(all.equal(gnz$ens_gene, rownames(out)))

out$genes <- gnz
out$genes$length <- rowMeans(res$length[rownames(out),])
s3saveRDS(out, "s3://mbl.data/mapping/may/fly/kallisto-DGEList.rds")


xxx = s3readRDS("s3://mbl.data/mapping/may/fly/kallisto-DGEList.rds")
