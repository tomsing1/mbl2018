library(readr)
library(stringr)
library(dplyr)
library(purrr)
library(aws.s3)
library(here)

kSampleKey <- here::here("fastq_processing", "may", "sample_key.txt")
kOutPath <- here::here("fastq_processing", "may", "sample_annotation.csv")
kS3Bucket <- "mbl.data"
kS3InPrefix <- "tmp/may"
kS3OutPrefix <- "reads/may"

kIgnore <- c("^dmunknown")  # prefixes for irrelevant samples
kSpecies <- c(mm = "mouse",
              dm = "fly",
              dr = "fish")
kSource <- c(ch = "cheek",
             pa = "palate",
             tg = "trigeminal",
             wp = "whisker pad",
             ne = "neuron",
             dg = "DRG",
             kr = "keratinocyte")
kGenotype <- c(ko = "knockout",
               wt = "wildtype",
               ao = "silencing channel",
               ap = "activating channel")
kTreatment <- c( hi = "histamine inhibition",
                 no = "no treatment",
                 at = "ATP activation")

Decode <- function(sample_id, species) {
  switch(species, 
         "mouse" = {
           data.frame(
             sample_id = sample_id,
             species = "mouse",
             source = kSource[str_sub(sample_id, 3, 4)],
             genotype = kGenotype[str_sub(sample_id, 5, 6)],
             treatment = NA,
             replicate = str_sub(sample_id, 7, 8),
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
                 treatment = kTreatment[str_sub(sample_id, 7, 8)],
                 replicate = str_sub(sample_id, 9, 10),
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
                  barcode = str_replace(barcode, "_", "-"),
                  Key = basename(Key))
}

DecodeAging <- function(sample_id) {
  data.frame(
    sample_id = sample_id,
    species = "fly",
    source = ifelse(grepl("^cd4a", sample_id), "neuron", "epidermis"),
    genotype = 'wt',
    treatment = ifelse(grepl("120h", sample_id), "120h", "72h"),
    replicate = NA,
    stringsAsFactors = FALSE
  )
} 

Main <- function() {
  
  # identify available FASTQ files
  fastq_files <- ListFastqFiles()
  
  # read keys for all samples
  all_samples <- read_tsv(kSampleKey,
                          col_types = cols(
                            sample_id = col_character(),
                            library = col_character(),
                            lane = col_integer(),
                            barcode = col_character()
                          ))
  
  # identify the two experiments
  aging <- dplyr::filter(all_samples, grepl("^c4da|^epi", sample_id,
                                            perl = TRUE))
  other <- dplyr::filter(all_samples, !grepl("^c4da|^epi", sample_id,
                                             perl = TRUE))
  
  #---- process experiment 'other'
  key <- other %>%
    dplyr::filter(!grepl(paste(kIgnore, collapse = "|"), sample_id)) %>%
    dplyr::mutate(species = kSpecies[stringr::str_sub(sample_id, 1, 2)])
  
  
  other_anno <- purrr::map2_df(key$sample_id, key$species, .f = Decode) %>%
    dplyr::mutate(fastq = paste0(species, "/", sample_id, "_R1.fastq.gz")) %>%
    dplyr::select(-species) %>%
    dplyr::left_join(key, by = "sample_id") %>%
    dplyr::left_join(fastq_files, by = c("barcode", "lane"))
  
  #---- process experiment 'aging'
  aging_anno <- DecodeAging(aging$sample_id) %>%
    dplyr::mutate(fastq = paste0(species, "/", sample_id, "_R1.fastq.gz")) %>%
    dplyr::left_join(all_samples, by = "sample_id") %>%
    dplyr::left_join(fastq_files, by = c("barcode", "lane"))
  
  sample_anno <- dplyr::bind_rows(other_anno, aging_anno) %>%
    write_csv(path = kOutPath)
  aws.s3::s3write_using(
    sample_anno, FUN = write_csv, 
    object = paste0(kS3OutPrefix, "/", "sample_annotations.csv"),
    bucket = kS3Bucket)
}

if (!interactive()) {
  Main()
}
