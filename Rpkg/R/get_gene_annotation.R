#' Retrieves the gene-level anntation file for this organism
#' 
#' @export
#' @param organism either "mouse", "fly", or "fish"
#' @return the gene-level annotation data.frame for the organism
get_gene_annotation <- function(organism = c("mouse", "fly", "fish")) {
  organism <- match.arg(organism)
  anno.base <- "s3://mbl.data/references/%s/gene_table.csv"
  gene.info <- s3read_using(read_csv, object = sprintf(anno.base, organism))
  
  gi <- as.data.frame(gene.info)
  gi <- rename_(gi, ensembl_gene_id = "ens_gene", symbol = "ext_gene",
                biotype = "transcript_biotype")
  gi <- subset(gi, !duplicated(ensembl_gene_id))
  rownames(gi) <- gi$ensembl_gene_id
  gi
}
