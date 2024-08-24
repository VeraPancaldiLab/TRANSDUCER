library(tidyverse)
library(edgeR)
library(biomaRt)
################################################################################
setwd("~/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/180122_Various/")

# Data loading

read_tsv("../../00_Data/Processed_data/rawHost_Cyt.tsv") -> host_cyt_raw
read_tsv("../../00_Data/Processed_data/rawTumor_Cyt.tsv") -> tumor_cyt_raw

cleanTumour <- c(
  "PDAC014T",
  "PDAC018T",
  "PDAC029T",
  "PDAC013T",
  "PDAC020T",
  "PDAC021T",
  "PDAC003T",
  "PDAC017T",
  "PDAC031T",
  "PDAC009T",
  "PDAC012T",
  "PDAC032T",
  "PDAC028T",
  "PDAC030T",
  "PDAC006T",
  "PDAC027T",
  "PDAC025T",
  "PDAC011T",
  "PDAC016T",
  "PDAC007T",
  "PDAC024T",
  "PDAC008T",
  "PDAC026T",
  "PDAC019T",
  "PDAC015T",
  "PDAC001T",
  "PDAC022T"
) # 27

cleanHost <- c(
  "PDAC018T",
  "PDAC029T",
  "PDAC013T",
  "PDAC020T",
  "PDAC021T",
  "PDAC017T",
  "PDAC031T",
  "PDAC012T",
  "PDAC032T",
  "PDAC028T",
  "PDAC030T",
  "PDAC027T",
  "PDAC025T",
  "PDAC011T",
  "PDAC024T",
  "PDAC008T",
  "PDAC026T",
  "PDAC019T",
  "PDAC015T",
  "PDAC022T"
) # 20


# Preprocessing
tumor_cyt_raw %>% column_to_rownames("EnsemblID") %>%
  dplyr::select(cleanTumour) %>%
  dplyr::filter_all(any_vars(. != 0)) %>%
  edgeR::DGEList() %>% edgeR::calcNormFactors(method = "TMM") %>%
  edgeR::cpm() ->tumor_cyt_norm # log transformation affects deconvolution performance, so this is non logged

host_cyt_raw %>% column_to_rownames("EnsemblID") %>%
  dplyr::select(cleanHost) %>%
  dplyr::filter_all(any_vars(. != 0)) %>%
  edgeR::DGEList() %>% edgeR::calcNormFactors(method = "TMM") %>%
  edgeR::cpm() -> host_cyt_norm # log transformation affects deconvolution performance, so this is non logged


# Ensembl ID fetch
ensembl_hs <- useEnsembl(biomart = "genes",
                         dataset = "hsapiens_gene_ensembl")
annot_ensembl_hs <- getBM(attributes = c('ensembl_gene_id',
                                         'external_gene_name'), mart = ensembl_hs)

ensembl_mm <- useEnsembl(biomart = "genes",
                          dataset = "mmusculus_gene_ensembl")
annot_ensembl_mm <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_name'), mart = ensembl_mm)

# Writing
tumor_cyt_norm %>% as_tibble(rownames = "EnsemblIDs") %>% write_tsv("01_Input/TumourCyt_foranalysis.tsv")
annot_ensembl_hs %>% write_tsv("01_Input/annotation_ensembl_hs.tsv")

host_cyt_norm %>% as_tibble(rownames = "EnsemblIDs") %>% write_tsv("01_Input/HostCyt_foranalysis.tsv")
annot_ensembl_mm %>% write_tsv("01_Input/annotation_ensembl_mm.tsv")
