library(mbl2018)
library(edgeR)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

# Data Preparation =============================================================
ym <- mbl_load_rnaseq("mouse", "may")
ym <- calcNormFactors(ym)

# We need to define the design of the experiment, ie. what defines the replciate
# groups? The samples are different based on:
# 1. "source": where was the sample taken from?
# 2. "genotype": is this a WT or knockout?
#
# The combincation of these two factors is stored in the "group" column:
View(ym$samples)
table(ym$samples$group)

# Prepping for t-tests (using linear models) ===================================

# use the values of group to define the design for the linear model.
design <- model.matrix(~ 0 + group, data = ym$samples)
colSums(design)

# Recall that t-tests assume "normality" somewhere.
# RNA-seq COUNT data violates this assumption, the "voom" function learns
# how we can account for these violations and still do t-tests.
#
# Before we use it, it is a good idea to remove lowly expressed genes from
# downstream analysis. The `filterByExpr()` function is provided by the edgeR
# package to help here.
#
# filterByExpr finds the genes that "have have sufficiently large counts to be
# retained in a statistal analysis". You can refer to its help page
# (?filterByExpr) for more information.
expressed <- filterByExpr(ym, design)
head(expressed)

ymf <-

vm <- voom(ym, design, plot = TRUE)

fit <- lmFit()

ymf <- ym[filterByExpr()]
