library(tidyverse)
library(edgeR)
library(biomaRt)
################################################################################
setwd("~/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/180122_Various/")

# Data loading

read_tsv("../../00_Data/Processed_data/rawHost_Cyt.tsv") -> cyt_raw

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
cyt_raw %>% column_to_rownames("EnsemblID") %>%
  dplyr::select(cleanHost) %>%
  dplyr::filter_all(any_vars(. != 0)) %>%
  edgeR::DGEList() %>% edgeR::calcNormFactors(method = "TMM") %>%
  edgeR::cpm() -> cyt_norm # log transformation affects deconvolution performance, so this is non logged

ensembl75 <- useEnsembl(biomart = "genes",
                          dataset = "mmusculus_gene_ensembl")


annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_name'), mart = ensembl75)

annot_ensembl75 %>% write_tsv("01_Input/annotation_ensembl.tsv")

cyt_norm %>% as_tibble(rownames = "EnsemblIDs") %>% write_tsv("01_Input/HostCyt_norm.tsv")
