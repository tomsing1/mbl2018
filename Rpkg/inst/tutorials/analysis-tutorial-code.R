# Assemble own kallisto dgelist ================================================
library(dplyr)
library(edgeR)

pheno <- read.csv("~/may-mouse-kallisto/sample-information.csv",
                  stringsAsFactors = FALSE, row.names = 1)
txinfo <- mbl_get_transcript_annotation("mouse")

library(tximport)

# Find files to load ===========================================================
# dirs <- dir("~/may-mouse-kallisto")
# dirs <- dir("~/may-mouse-kallisto", pattern = "mm")
# files <- paste("~/may-mouse-kallisto", dirs, "abundance.tsv", sep = "/")


# I will cheat because I want the files to be processed in the same order
# that they appear in my `pheno` table
files2 <- paste("~/may-mouse-kallisto",
                pheno$sample_id,
                "abundance.tsv",
                sep ="/")
file.exists(files2)

XXX <- select(txinfo,target_id,ens_gene)
txi <- tximport(files2, type = "kallisto",
                tx2gene = XXX, ignoreTxVersion = TRUE,
                countsFromAbundance = "lengthScaledTPM")
is(txi)
names(txi)

head(txi$counts)
dim(txi$counts)
nrow(txi$counts)

head(txi$length)

# want to label columns with their sample id
colnames(txi$counts) <- pheno$sample_id
my <- DGEList(txi$counts, samples = pheno)

# Let's put the gene data.frame in there!
# where is the gene data.frame?
#
#   filter() to select cases based on their values.
#   arrange() to reorder the cases.
#   select() and rename() to select variables based on their names.
#   mutate() and transmute() to add new variables that are functions of existing variables.
#   summarise() to condense multiple values to a single value.
#   sample_n() and sample_frac() to take random samples.
#   group_by()

ginfo <- distinct(txinfo, ens_gene, .keep_all = TRUE)
head(txi$counts)

# how to index R objects
ginfo[1:3,]
ginfo[1:3, 1:2]

ginfo[1:3, c("ens_gene", "ext_gene")]
head(rownames(ginfo))

rownames(ginfo) <- ginfo$ens_gene
head(rownames(ginfo))

gsorted <- ginfo[rownames(txi$counts), ]

all.equal(rownames(gsorted), rownames(txi$counts))


my <- DGEList(txi$counts, samples = pheno, genes = gsorted)
saveRDS(my, "~/mouse-DGEList.rds")

some <- my[100:110, 1:3]

is(my)


gi <- txinfo %>%
  group_by(ens_gene) %>%
  summarize(ext_gene = ext_gene[1], ntx = n(), biotype = transcript_biotype[1]) %>%
  as.data.frame

# Are the genes in gi in the same order as the rows of `my`?
head(rownames(my))
head(gi$ens_gene)

# There are many ways to fix this.
# You can index the rows of a 2d object (data.frame) by index:
gi[c(1, 3, 5),]

# Or you can index by rowname.
# Let's set the rownames() of our thing to the ensembl gene ids
rownames(gi) <- gi$ens_gene
gi[c(1, 3, 5),]

gi[head(rownames(my)),]

##all.equal??
