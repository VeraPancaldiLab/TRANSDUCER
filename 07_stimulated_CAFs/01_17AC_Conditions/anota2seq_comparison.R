library(tidyverse)
library(anota2seq)
################################PARAMETERS######################################
filter_samples = "17AC" # NULL | 17AC | 02136
exclude_samples = c("Batch_A_17AC_FAKi_Input", "Batch_A_17AC_TGF_F8")
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/07_stimulated_CAFs/01_17AC_Conditions/")

# Data loading
counts <- read_tsv("../00_Data/rawcounts.tsv")
manip_info <- read_tsv("../00_Data/Data_RNA_sample_Jacobo.tsv")

## filter of samples
if (!is.null(filter_samples)){
  counts <- dplyr::select(counts, Geneid, names(counts)[str_detect(names(counts), pattern = filter_samples)])
}

if (!is.null(exclude_samples)){
  counts <- dplyr::select(counts, Geneid, names(counts)[!names(counts) %in% exclude_samples])
}