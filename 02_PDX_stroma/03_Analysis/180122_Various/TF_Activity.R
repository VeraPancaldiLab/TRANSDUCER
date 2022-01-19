library(tidyverse)
library(dorothea)
library(viper)
################################################################################
setwd("~/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/180122_Various/")

# Data loading
read_tsv("01_Input/HostCyt_norm.tsv") -> cyt_norm
cyt_genenames <- cyt_norm
read_tsv("01_Input/annotation_ensembl.tsv") -> annot_ensembl75

translate <- deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_name")])

cyt_norm$EnsemblIDs %>% translate[.] -> cyt_genenames$EnsemblIDs
cyt_genenames <- cyt_genenames[!is.na(cyt_genenames$EnsemblIDs),]
print(paste(round(100*nrow(cyt_genenames)/nrow(cyt_norm),2), "% of genes kept after ensemblID translation"))

# TF activity
## load Dorothea Regulons
### GTex
data(dorothea_mm, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))
vipertype = "Gtex"

### PANcancer
# data(dorothea_hs_pancancer, package = "dorothea")
# regulons <- dorothea_hs_pancancer %>%
#   dplyr::filter(confidence %in% c("A", "B","C"))
# vipertype = "PANcancer"

## VIPER
minsize = 5
ges.filter = FALSE

tf_activities_stat <- dorothea::run_viper(cyt_genenames, regulons,
                                          options =  list(minsize = minsize, eset.filter = ges.filter, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))

