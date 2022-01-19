library(tidyverse)
library(dorothea)
library(viper)
################################################################################
setwd("~/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/180122_Various/")

# Data loading
read_tsv("01_Input/HostCyt_norm.tsv") -> cyt_norm
cyt_genenames_ <- cyt_norm
read_tsv("01_Input/annotation_ensembl.tsv") -> annot_ensembl

translate <- deframe(annot_ensembl[c("ensembl_gene_id", "external_gene_name")])

cyt_norm$EnsemblIDs %>% translate[.] -> cyt_genenames_$EnsemblIDs
cyt_genenames_ <- cyt_genenames_[!is.na(cyt_genenames_$EnsemblIDs),]
distinct(cyt_genenames_, EnsemblIDs, .keep_all=TRUE) -> cyt_genenames
print(paste(round(100*nrow(cyt_genenames)/nrow(cyt_norm),2), "% of genes kept after ensemblID translation and replicate removal"))


# TF activity
## load Dorothea Regulons
### GTex
data(dorothea_mm, package = "dorothea")
regulons <- dorothea_mm %>%
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

cyt_genenames %>% column_to_rownames("EnsemblIDs") %>% dorothea::run_viper(regulons,
                                          options =  list(minsize = minsize, eset.filter = ges.filter, 
                                                          cores = 1, verbose = FALSE, nes = TRUE)) -> tf_activities_stat

tf_activities_stat %>% as_tibble(rownames = "TF") %>% write_tsv("01_Input/TFact_stroma_Gtex.tsv")
