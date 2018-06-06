library(mbl2018)
library(edgeR)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

# Data Preparation =============================================================
ym <- mbl_load_rnaseq("mouse", "may")

# This `is` function tells you what something is. From left to right, we go
# from "most specific thing" to least, just like we are:
# Humans -> Mammals -> Bipeds -> Vertebrates -> Mammal -> Animals
is(ym)

# Recall what a DGEList object has in it

# Brief discussion of

ym <- calcNormFactors(ym)
dim(ym)

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

# What do the rows and columns mean of the design?
# How many rows should there be? How many columns?
View(de)
colSums(design)

# Before we continue, it is a good idea to remove lowly expressed genes from
# downstream analysis. The `filterByExpr()` function is provided by the edgeR
# package to help here.
#
# filterByExpr finds the genes that "have have sufficiently large counts to be
# retained in a statistal analysis". You can refer to its help page
# (?filterByExpr) for more information.
expressed <- filterByExpr(ym, design)
head(expressed)  # what does expressed look like (a logical vector)
summary(expressed) # what is the distrubtion of TRUE/FALSE?

# We can index into the DGEList as a 2d object
#   * First dimension is rows (genes)
#   * Second dimenstion is samples (columns)
# We want to keep only the rows (genes) that are expressed, and all of the
# samples (Note the index for the second dimension is blank)
ymf <- ym[expressed, ]

# Recall that t-tests assume "normality" somewhere.
# RNA-seq COUNT data violates this assumption, the "voom" function learns
# how we can account for these violations and still do t-tests.
vm <- voom(ymf, design, plot = TRUE)
is(vm)
names(vm)

# Let's now fit the linear model.
# To do this, we need the data (vm$E) and the design (vm$design). the $weights
# are used here.
#
# When the design is parameterized in the way we have descriged (~ 0 + group)
# the linear model learns the mean expression of each gene in each group.
fit <- lmFit(vm, vm$design)

# what's in `fit`?
names(fit)
# Look at ?lmFit to learn what these things are (under the Values section)

# The mean expression values per group are stored in the vm$coefficients
# matrix
head(fit$coefficients)

# Let's see if this is true. I made a helper function to create a plot that
# shows the expression values per group
mbl_plot_expression(vm, "ENSMUSG00000000037", "group")

# STEP 1 FOR T-TEST IS COMPLETE
# We know the mean expression of each gene per group

# Now that we have learned the mean expression, we just have to define the
# effect we want to measure.
# Recall that our data is in log-space now, so if you want to divide you
# need to actually subtract.

# So, to estimate the fold change of the cheek knocout we subtract the
# wildtype cheek expression value from the knockout value.
#
# You can manually do it like this:
head(fit$coefficients[, 1] - fit$coefficients[, 2])

# Adding and subtracting across the columns of the coefficient is called
# "defining a contrast", so let's do that:
cm <- makeContrasts(
  cheek_KO = groupcheek_knockout - groupcheek_wildtype,
  levels = vm$design)

# Now run the contrast
fit2 <- contrasts.fit(fit, cm)

# STEP 2 IS COMPLETE. We now have measured the difference in mean expression
# of the contrasts defined in `cm`
head(fit2$coefficients)

# This part is magic: we are trying to compensate for the fact that we have
# few replicates, 2 in one condition, 3 in the error, and being able to estimate
# the variance of our expression estimates with these small observations is
# not robust
fit2 <- eBayes(fit2)

# Lets get the results!
res <- topTable(fit2, "cheek_KO", number = Inf)
head(res)

mbl_plot_expression(vm, "ENSMUSG00000068075", "group")

# What does the pvalue distribution look like?
ggplot(res, aes(P.Value)) +
  geom_histogram(bins = 50)

# Let's take a birds eye view of our data with Principal Components Analysis.
#
# What is PCA? (whiteboard)
#
# We have provided a simple PCA function to make this easy
pca <- mbl_pca(vm)


# What type of contrasts whould give results?
# Try one for yourself.




# GSEA =========================================================================
library(multiGSEA)
library(multiGSEA.shiny)

gdb <- getMSigGeneSetDb(c("h", "c2", "c5"), species = "mouse",
                        id.type = "ensembl")
mgko <- multiGSEA(gdb, vm, vm$design, cm[, "cheek_KO"], methods = "camera")




