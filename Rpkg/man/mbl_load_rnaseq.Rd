% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mbl_load_rnaseq.R
\name{mbl_load_rnaseq}
\alias{mbl_load_rnaseq}
\title{Loads RNA-seq dataset generated for MBL2018 Neurobiology Course.}
\usage{
mbl_load_rnaseq(organism = c("mouse", "fly", "fish"), dataset = c("all",
  "may", "mbl"), ...)
}
\arguments{
\item{organism}{either "mouse", "fly", or "fish"}

\item{dataset}{either `"all"`, `"mbl"`, or `"may"`.}
}
\value{
a DGEList of the RNAseq data.
}
\description{
Curated datasets are stored on the `s3://mbl.data` Amazon S3 bucket.
This function allows you to load the different datasets we generated for the
model organisms studied here

By default this will return a DGEList that has both the data generated
prior to the course (the "may" dataset), and the data you generated here
into one object. There is a `"dataset"` column in the `$samples` data.frame
that can be used to split the two datasets.

Alternatively you can load on specific dataset or another by specifing
the `dataset` parameter to the function.
}
\examples{
# To load all of the mouse data:
ym.all <- mbl_load_rnaseq("mouse", dataset = "all")

# Load only the "mbl" mouse dataset
ym.mbl <- mbl_load_rnaseq("mouse", dataset = "mbl")
}
