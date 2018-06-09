# mbl2018

This repository contains the materials used for the genomics portion of the
Neurobiology course at MBL, June 2018.

In here you will find:

1. A number of scripts used to setup the amazon machines used for processing
   the FASTQ files from QC to quantitation with Kallisto, as well as being
   configured with RStudio for data analysis.

2. An `Rpkg` subdirectory, which contains the `mbl2018` R package used during
   the workshop.

## Getting Help

Please feel free to contact either Thomas Sandmann or Steve Lianoglou if you
need help with any of the topics covered during our workshop.

## The mbl2018 R package

The `mbl2018` R package contains utility functions to help with:

1. Loading the expression data (`mbl_load_rnaseq()`).
2. Running a PCA analysis on your data (`mbl_pca()`).
3. Convenience function to plot one gene at a time (`mbl_plot_expression()`)
4. Converting your expression data into a (huge) "tidy table" so that you can
   use ggplot2 to create arbitrarily awesome visualizations (`mbl_tidy()`)

You will also find a number of tutorials found in the `inst/tutorials`
directory (see below).

After you load the `mbl2018` library using the `library()` function
(see below), you can find more help on how to use these functions by
prefixing their names with `?` in your R terminal and reading their help
page, ie.

```r
library(mbl2018)
?mbl_load_rnaseq
```

Will bring up a help page for thie `mbl_load_rnaseq` function.

### Package Installation

You will first need to ensure that you have a few helper packages installed
in R  so that you can easily install this package from github. Run these
commands in your R workspace to ensure that the installation of `mbl2018`
goes smoothly:

```r
install.packages(c("remotes", "devtools", "gplots"))
```

You can now install the `mbl2018` package by running these commands in your
R workspace:

```r
source("https://bioconductor.org/biocLite.R")
biocLite("tomsing1/mbl2018", subdir = "Rpkg")
```

After you install this package, you can make its functions available to you
by loading it with the `library()` function.

```r
library(mbl2018)
```

You are now ready to follow the tutorials and "recipes" found in the
`inst/tutorials` directory.

### Tutorials

There are a number of tutorial files in the `inst/tutorials` directory.

Perhaps the mose useful files for your study after the course will be:

1. The `inst/tutorials/exploratory-data-analysis-recipes.Rmd` file, which
   takes you through the steps required to load data, perform a PCA and
   differential expression analysis, and plot some results.
2. The `inst/tutorials/analysis-rnaseq-tutorial.Rmd` file was not completed.
   It does, however, include:
   i. a diagram of the `DGEList` data structure, to remind you of how the
      genomics data is stored and made accessible to you in the objects that are
      returned from the `mbl_load_rnaseq()` function.
   ii. links to other gene expression analysis tutorials online that are
       *extremely* helpful to read, and were used as inspiration for some
       of the topics covered in our tutorials and workshop.
3. The `inst/tutorials/dge-mouse-ko-short-sweet.R` file is a "slim" R script
   that performs a simple differential expression analysis. It includes comments
   that point out which steps of the R code correspond to the three basic
   steps of the "t-test workflow" we were hammering on during the workshop.
4. The `inst/tutorials/gsea-mouse-wt.R` script shows you how to do a 
   differential gene expression analysis and gene set enrichment analysis using
   the [multiGSEA](https://github.com/lianos/multiGSEA) package.

