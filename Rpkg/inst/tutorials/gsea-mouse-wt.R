# This script shows you how to run a differential expression and GSEA analysis
# using multiGSEA.
#
# There are instructions at the end of this script to show you how to install
# multiGSEA on your machine: it's easy!
library(mbl2018)
library(edgeR)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

# Analysis of old vs young mice
ym.all <- mbl_load_rnaseq("mouse", "mbl")
ym <- ym.all[, ym.all$samples$genotype == "wildtype"]
ym <- calcNormFactors(ym)

design <- model.matrix(~ 0 + group, ym$samples)

ymf <- ym[filterByExpr(ym, design),]

# Step 0a: Define the questions you want to ask. You are restricted to asking
#          questions using only the terms found in the columns of `design`:
#          colnames(design)
colnames(design)
cm <- makeContrasts(
  cheek_palate = groupcheek_wildtype - grouppalate_wildtype,
  levels = design)

# Step 0b: Use voomWithQualityWeights to prep data
vm <- voomWithQualityWeights(ymf, design, plot = TRUE)


# Step 1: Fit the mean expression of each gene per group in your experiment.
fit <- lmFit(vm, design)

# Step 2: Calculate differences in means for your questions in "cm"
fit2 <- contrasts.fit(fit, cm)

# Step 3: do something about the variance
fit2 <- eBayes(fit2)

reschk <- topTable(fit2, "cheek_palate", n = Inf)

# multiGSEA ====================================================================
library(multiGSEA)
library(multiGSEA.shiny)

# Get Gene Sets
gdb <- getMSigGeneSetDb(c("h", "c2", "c3", "c5", "c7"), "mouse", "ensembl")

mg <- multiGSEA(gdb, vm, vm$design, cm[, "cheek_palate"],
                methods = c("camera", "goseq"),
                min.gs.size = 5,
                feature.bias = setNames(vm$genes$length, rownames(vm)))

explore(mg)

# Installation =================================================================
# If multiGSEA and friends aren't installed on your machine, you should be
# able to install them like so:
install.packages("remotes")
source("http://bioconductor.org/biocLite.R")
biocLite(
  c("lianos/multiGSEA.shiny", "lianos/GeneSetDb.MSigDB.Mmusculus.v61"))

