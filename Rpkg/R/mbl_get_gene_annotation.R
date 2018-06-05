#' Retrieves the gene-level anntation file for this organism
#'
#' @export
#' @param organism either "mouse", "fly", or "fish"
#' @return the gene-level annotation data.frame for the organism
mbl_get_gene_annotation <- function(organism = c("mouse", "fly", "fish"),
                                    rm.description = TRUE) {
  # Ensure user asked for a valid organism
  organism <- match.arg(organism)

  # Retrieve transcript annotation table and manipulate it into one suitable
  # for gene-level annotations.
  tx.info <- mbl_get_transcript_annotation(organism, rm.description = FALSE)

  # This code block induces an ordering on the "transcript_biotype" columns
  # such that when the table is ordered by this column (via the `arrange()`
  # call further down, "protein_coding" comes first, "lincRNA" second, and
  # so on.
  btype <- c("protein_coding", "lincRNA", "snoRNA", "miRNA")
  btype <- unique(c(btype, tx.info$transcript_biotype))
  tx.info$transcript_biotype <- factor(tx.info$transcript_biotype, btype)

  # Now we aggregate the transcript table by gene identifier, and only take
  # one row per gene. This code block uses the "dplyr" package to manipulate
  # the tx.info data.frame into the final gene level data.frame we want.
  gi <- tx.info %>%
    arrange(ens_gene, transcript_biotype) %>%
    group_by(ens_gene) %>%
    mutate(n_txs = n()) %>%
    ungroup %>%
    distinct(ens_gene, .keep_all = TRUE) %>%
    rename(biotype = transcript_biotype) %>%
    select(ens_gene, ext_gene, n_txs, biotype, description) %>%
    as.data.frame

  rownames(gi) <- gi$ens_gene

  if ("description" %in% colnames(gi) && rm.description) {
    gi$description <- NULL
  }

  gi
}

#' @export
#' @rdname mbl_get_gene_annotation
#' @importFrom readr read_csv read_tsv
mbl_get_transcript_annotation <- function(organism = c("mouse", "fly", "fish"),
                                          rm.description = TRUE) {
  organism <- match.arg(organism)
  anno.base <- "s3://mbl.data/references/%s/gene_table.csv"
  t.info <- s3read_using(read_csv, object = sprintf(anno.base, organism))

  ti <- as.data.frame(t.info)
  rownames(ti) <- ti$target_id
  if ("description" %in% colnames(ti) && rm.description) {
    ti$description <- NULL
  }

  ti
}
