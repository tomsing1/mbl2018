library(mbl2018)
library(edgeR)
library(dplyr)
library(ggplot2)
theme_set(theme_classic())

ym <- mbl_load_rnaseq("fly", "all")
ym <- calcNormFactors(ym)

# Find genes differentially regulated in neurons during development. ===========
samples.ndev <- c("neuron_wildtype_72h", "neuron_wildtype_120h")
ym.ndev <- ym[, ym$samples$group %in% samples.ndev]
ym.ndev$design <- model.matrix(~ 0 + group, ym.ndev$samples)
colnames(ym.ndev$design) <- sub("group" ,"", colnames(ym.ndev$design))

ym.ndev <- ym.ndev[filterByExpr(ym.ndev, ym.ndev$design),]
vm.ndev <- voomWithQualityWeights(ym.ndev, ym.ndev$design)

cm <- makeContrasts(
  neuro_age = neuron_wildtype_120h - neuron_wildtype_72h,
  levels = vm.ndev$design)

res.ndev <- lmFit(vm.ndev, vm.ndev$design) %>%
  contrasts.fit(cm) %>%
  eBayes() %>%
  topTable("neuro_age", n = Inf)

# Find genes

res.ndev <-
ym.neuro.deb <-
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
