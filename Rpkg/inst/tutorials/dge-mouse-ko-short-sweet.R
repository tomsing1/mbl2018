# This script is a companion to the heavily commented and verbose
# `dge-mouse-ko.R` script in this directory.
#
# This file has minimal comments, the point of which are to show that:
#
# 1. these moves aren't really complex
# 2. we are essentially performing the same three steps we explained we were
#    doing with a t-test.
#
# Recall that we genarlized the steps required to perform a t-test as:
#
# 1. Estimate the mean expression of the genes in your experimental groups.
# 2. Calculate the differences in the mean expression in your
#    comparisons (contrasts) of interest.
# 3. Incorporate the variance of your measurements to estimate statistical
#    significance of the
#
# Keep an eye out to see where these are happening in this script.
#
# Finally, at the very end, there is an extremely condensed version of the
# code.
library(mbl2018)
library(edgeR)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

# Fetch the data (DGEList) for the experiment
ym <- mbl_load_rnaseq("mouse", "may")

# Account for differences in read depth among the samples
ym <- calcNormFactors(ym)

# Define the design of the experiment in a way that allows us to estimate the
# mean expression of each gene per replicate group/condition:
design <- model.matrix(~ 0 + group, data = ym$samples)

# Filter out lowly expressed genes
expressed <- filterByExpr(ym, design)
ymf <- ym[expressed, ]

# Call the voom function which will provide us an object that is similar
# to a DGEList, but can now be analyzed using linear models. Note that in the
# following `vm` object:
#
# 1. the expression values in `vm` object have been log2 transformed.
# 2. the `$samples`` data.frame is now called `$targets`.
vm <- voom(ymf, design, plot = TRUE)

# STEP 1 OF T-TEST "WORKFLOW"
# ---------------------------
# estimate the mean expression of the genes in our replicate groups, as they
# are defined in the `vm$design` matrix.
fit <- lmFit(vm, vm$design)

# STEP2a OF T-TEST WORKFLOW
# -------------------------
# Define the differences in expression you want to evaluate by writing a
# mathematical function using the terms found in `colnames(fit$coefficients)`.
#
# Note that you can define any number of comparisons you are interested in
# analyzing. This command defines two comparisons to estimate:
#
# 1. The log fold change of gene expression in the cheek_knockout samples vs
#    the cheek_wildtype samples. You will reference this comparsion by the name
#    "cheek_KO"
# 2. The log fold change of gene expression between wildtype cheek samples vs
#    the wildetype palate samples. You will reference this comparison by the
#    name "cheek_palate".
cm <- makeContrasts(
  cheek_KO = groupcheek_knockout - groupcheek_wildtype,
  cheek_palate = groupcheek_wildtype - grouppalate_wildtype,
  levels = vm$design)

# STEP 2b OF T-TEST WORKFLOW
# --------------------------
# now calculate the differences in the mean expression of the groups in the
# comparisons (contrasts) enumerated above
fit2 <- contrasts.fit(fit, cm)

# STEP 3: OF T-TEST WOFKLOW
# -------------------------
# This call below does a bit of magic to estimate the variance in our data.
# It uses all the data to get a better estimate than the naive approach of just
# using the measurements from each gene individually.
fit2 <- eBayes(fit2)

# PARY TIME
# ---------
# Now we can get the results we worked hard for.

# This command gets our cheek_KO results. res.ko will have a data.frame that
# has for all your genes their logFC, pvalues, adj.P.Val (FDR), etc.
res.ko <- topTable(fit2, "cheek_KO", number = Inf)

# This command gets our tissue-vs-tissue comparison:
res.wt <- topTable(fit2, "cheek_palate", number = Inf)

# FINISHED =====================================================================

# Now super short and sweet

# Load and prep data
ym <- mbl_load_rnaseq("mouse", "may")
ym <- calcNormFactors(ym)

# Define experimental design and remove lowly expressed genes
design <- model.matrix(~ 0 + group, data = ym$samples)

# Filter out lowly expressed genes
expressed <- filterByExpr(ym, design)
ymf <- ym[expressed, ]

# Use voom to log transform data and prep it for use with linear models
# vm <- voom(ymf, design, plot = TRUE)

# Instead of normal voom, we are using voomWithQualityWeights to do our best
# to be robust to low quality samples
cols <- mbl_create_color_map(ym$samples$group) # This was not covered during
cols <- cols[ym$samples$group]                 # the tutorial
vm <- voomWithQualityWeights(ymf, design, col = cols, plot = TRUE)

# Define the questions you want to ask using the columns of your design matrix.
# ie. what comparisons to you want to measure the biological effect and
# p-values for
cm <- makeContrasts(
  cheek_KO = groupcheek_knockout - groupcheek_wildtype,
  cheek_palate = groupcheek_wildtype - grouppalate_wildtype,
  levels = vm$design)

fit <- lmFit(vm, vm$design)    # Step 1: estimate mean expression of genes
fit2 <- contrasts.fit(fit, cm) # Step 2: calculate mean differences across groups
fit2 <- eBayes(fit2)           # Step 3: deal with the variance in our data

# This command gets our cheek_KO results. res.ko will have a data.frame that
# has for all your genes their logFC, pvalues, adj.P.Val (FDR), etc.
res.ko <- topTable(fit2, "cheek_KO", number = Inf)

# This command gets our tissue-vs-tissue comparison:
res.wt <- topTable(fit2, "cheek_palate", number = Inf)
