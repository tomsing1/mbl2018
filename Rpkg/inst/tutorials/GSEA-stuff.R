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
ym.all <- mbl_load_rnaseq("mouse", "all")
# ym <- ym[, ym$samples$group %in% c("cheek_wildtype_Old", "cheek_wildtype_Young")]
age.samples <- grepl("Old|Young", ym$samples$group)
ym <- ym[, age.samples]
ym <- calcNormFactors(ym)

design <- model.matrix(~ 0 + group, ym$samples)
ymf <- ym[filterByExpr(ym, design),]

# Step 0a: Define the questions you want to ask. You are restricted to asking
#          questions using only the terms found in the columns of `design`:
#          colnames(design)
colnames(design)
cm <- makeContrasts(
  age_cheek = groupcheek_wildtype_Old - groupcheek_wildtype_Young,
  age_trig = grouptrigeminal_wildtype_Old - grouptrigeminal_wildtype_Young,
  age = ((grouptrigeminal_wildtype_Old + groupcheek_wildtype_Old) / 2) -
    ((grouptrigeminal_wildtype_Young + groupcheek_wildtype_Young) / 2),
  levels = design)

# Step 0b: Use voomWithQualityWeights to prep data
vm <- voomWithQualityWeights(ymf, design, plot = TRUE)


# Step 1: Fit the mean expression of each gene per group in your experiment.
fit <- lmFit(vm, design)

# Step 2: Calculate differences in means for your questions in "cm"
fit2 <- contrasts.fit(fit, cm)

# Step 3: do something about the variance
fit2 <- eBayes(fit2)

reschk <- topTable(fit2, "age_cheek", n = Inf)
restrig <- topTable(fit2, "age_trig", n = Inf)
resall <- topTable(fit2, "age", n = Inf)
# multiGSEA ====================================================================
library(multiGSEA)
library(multiGSEA.shiny)

# Get Gene Sets
gdb <- getMSigGeneSetDb(c("h", "c2", "c5", "c7"), "mouse", "ensembl")

mg <- multiGSEA(gdb, vm, vm$design, cm[, "age"], methods = "camera",
                min.gs.size = 5)


