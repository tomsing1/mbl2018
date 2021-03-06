library(readr)
library(stringr)
library(dplyr)
library(purrr)
library(aws.s3)
library(here)
library(tidyr)

kSampleKey <- here::here("fastq_processing", "mbl", "sample_key.txt")
kOutPath <- here::here("fastq_processing", "mbl", "sample_annotation.csv")
kS3Bucket <- "mbl.data"
kS3InPrefix <- "tmp/mbl"
kS3OutPrefix <- "reads/mbl"

kSpecies <- c(mm = "mouse",
              dm = "fly",
              dr = "fish")
kSource <- c(ch = "cheek",
             pa = "palate",
             tg = "trigeminal",
             wp = "whisker pad",
             ne = "neuron",
             dg = "DRG",
             kr = "keratinocyte",
             ep = "epidermis",
             ed = "epidermis_dorsal",
             ev = "epidermis_ventral")
kGenotype <- c(ko = "knockout",
               wt = "wildtype",
               co = "cora",
               pi = "PI4K",
               ed = "eda_mutant",
               ko = "CKO",
               yo = NA,
               ol = NA
               )
kTreatment <- c(
  yo = "Young",
  ol = "Old"
)

Decode <- function(sample_id, species) {
  switch(species, 
         "mouse" = {
           data.frame(
             sample_id = sample_id,
             species = "mouse",
             source = kSource[str_sub(sample_id, 3, 4)],
             genotype = kGenotype[str_sub(sample_id, 5, 6)],
             replicate = str_sub(sample_id, 7, 8),
             treatment = kTreatment[str_sub(sample_id, 5, 6)],
             stringsAsFactors = FALSE
           )},
           "fish" = {
             data.frame(
               sample_id = sample_id,
               species = "fish",
               source = kSource[str_sub(sample_id, 3, 4)],
               genotype = kGenotype[str_sub(sample_id, 5, 6)],
               treatment = NA,
               replicate = str_sub(sample_id, 7, 8),
               stringsAsFactors = FALSE
             )},
             "fly" = {
               data.frame(
                 sample_id = sample_id,
                 species = "fly",
                 source = kSource[str_sub(sample_id, 3, 4)],
                 genotype = kGenotype[str_sub(sample_id, 5, 6)],
                 replicate = str_sub(sample_id, 9, 10),
                 treatment = NA,
                 stringsAsFactors = FALSE
               )
         })
}

ListFastqFiles <- function() {
  fastq <- get_bucket_df(kS3Bucket, prefix = kS3InPrefix) %>%
    dplyr::filter(grepl("fastq.gz$", Key),
                  grepl("unmapped.1.fastq.gz$", Key)) %>%
    dplyr::select(Key) %>%
    dplyr::mutate(barcode = stringr::str_split_fixed(Key, fixed("."), 
                                                     n = 6)[, 3],
                  lane = as.integer(stringr::str_split_fixed(Key, fixed("."), 
                                                   n = 6)[, 2]),
                  Key = basename(Key))
}

Main <- function() {
  
  key <- read_tsv(kSampleKey, col_types = cols(
    sample_id = col_character(),
    library = col_character(),
    lane = col_integer(),
    barcode = col_character()
  )) %>%
    dplyr::mutate(species = kSpecies[stringr::str_sub(sample_id, 1, 2)])

  fastq_files <- ListFastqFiles()
  
  sample_anno <- purrr::map2_df(key$sample_id, key$species, .f = Decode) %>%
    dplyr::mutate(fastq = paste0(species, "/", sample_id, "_R1.fastq.gz")) %>%
    dplyr::select(-species) %>%
    dplyr::left_join(key, by = "sample_id") %>%
    dplyr::left_join(fastq_files, by = c("barcode", "lane")) %>%
    tidyr::replace_na(list(genotype = "wildtype"))
  write_csv(sample_anno, path = kOutPath)
  aws.s3::s3write_using(
    sample_anno, FUN = write_csv, 
    object = paste0(kS3OutPrefix, "/", "sample_annotations.csv"),
    bucket = kS3Bucket)
}

if (!interactive()) {
  Main()
}
