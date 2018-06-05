# Assemble own kallisto dgelist
library(dplyr)
library(edgeR)

pheno <- read.csv("~/may-mouse-kallisto/sample-information.csv",
                  stringsAsFactors = FALSE, row.names = 1)
txinfo <- mbl_get_transcript_annotation("mouse")

# Find files to load ===========================================================
dirs <- dir("~/may-mouse-kallisto")
dirs <- dir("~/may-mouse-kallisto", pattern = "mm")
files <- paste("~/may-mouse-kallisto", dirs, "abundance.tsv", sep = "/")


# I will cheat because I want the files to be processed in the same order
# that they appear in my `pheno` table
files2 <- paste("~/may-mouse-kallisto",
                pheno$sample_id,
                "abundance.tsv",
                sep ="/")
file.exists(files2)

txi <- tximport(files2, type = "kallisto",
                tx2gene = tx2gene, ignoreTxVersion = TRUE,
                countsFromAbundance = "lengthScaledTPM")
is(txi)
names(txi)
head(txi$counts)
head(txi$length)

# want to label columns with their sample id
colnames(txi$counts) <- pheno$sample_id
my <- DGEList(txi$counts, samples = pheno)

# Let's put the gene data.frame in there!
# where is the gene data.frame?
gi <- txinfo %>%
  group_by(ens_gene) %>%
  summarize(ext_gene = ext_gene[1], ntx = n(), biotype = transcript_biotype[1]) %>%
  as.data.frame





# Are they in the same order?
# use rownames() <-
rownames(gi) <- gi$ens_gene

my$genes <- gi[rownames(my), ]
all.equal(my$genes$ens_gene, rownames(my))


my <-
files <- paste("~/may-mouse-kallisto", )
# Need
y <- mbl_load_rnaseq(organism = "mouse", dataset = "may")
y <- calcNormFactors(y)

des0 <- model.matrix(~ 0 + group, y$samples)
num.expr.per.row <- rowSums(y$counts == 0)

vma <- voom(y, des0, plot = TRUE)

vmf <-
