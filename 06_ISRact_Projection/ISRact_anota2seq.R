library(tidyverse)
library(anota2seq)
library(anota2seqUtils)
library(biomaRt)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/081221_TranslationEfficacy/")

# Data loading
normTumor_cyt <- read_tsv("../../00_Data/Processed_data/normTumor_Cyt.tsv")
normTumor_pol <- read_tsv("../../00_Data/Processed_data/normTumor_Pol.tsv")
sample_info <- read_tsv("../../00_Data/Processed_data/sample_info.tsv") %>% column_to_rownames("sample")

# Joint dataframe creation
sn_cyt <- paste(colnames(normTumor_cyt[-1]), "cyt", sep = "_")
sn_pol <- paste(colnames(normTumor_pol[-1]), "pol", sep = "_")

normTumor_cyt.tmp <- normTumor_cyt
normTumor_pol.tmp <- normTumor_pol

colnames(normTumor_cyt.tmp) <- c("EnsemblID",sn_cyt)
colnames(normTumor_pol.tmp) <- c("EnsemblID",sn_pol)

norm_tumor <- inner_join(normTumor_cyt.tmp, normTumor_pol.tmp, by= "EnsemblID")

# Translate to geneID
## get gene ID with Biomart
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)

annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id'), mart = ensembl75)

translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

norm_tumor <- dplyr::mutate(norm_tumor, geneID = translate["EnsemblID"])
