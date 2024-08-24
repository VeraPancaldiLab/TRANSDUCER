library(tidyverse)
library(dorothea)
library(viper)
################################################################################
setwd("~/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/180122_Various/")

# Data loading
read_tsv("01_Input/TumourCyt_foranalysis.tsv") -> cyt_tumor
cyt_tumor_genenames_ <- cyt_tumor
read_tsv("01_Input/annotation_ensembl_hs.tsv") -> annot_ensembl_hs

read_tsv("01_Input/HostCyt_foranalysis.tsv") -> cyt_host
cyt_host_genenames_ <- cyt_host
read_tsv("01_Input/annotation_ensembl_mm.tsv") -> annot_ensembl_mm

# Gene name conversion
translate_hs <- deframe(annot_ensembl_hs[c("ensembl_gene_id", "external_gene_name")])

cyt_tumor$EnsemblIDs %>% translate_hs[.] -> cyt_tumor_genenames_$EnsemblIDs
cyt_tumor_genenames_ <- cyt_tumor_genenames_[!is.na(cyt_tumor_genenames_$EnsemblIDs), ]
distinct(cyt_tumor_genenames_, EnsemblIDs, .keep_all = TRUE) -> cyt_tumor_genenames
print(paste(round(100 * nrow(cyt_tumor_genenames) / nrow(cyt_tumor), 2), "% of Human genes kept after ensemblID translation and replicate removal"))


translate_mm <- deframe(annot_ensembl_mm[c("ensembl_gene_id", "external_gene_name")])

cyt_host$EnsemblIDs %>% translate_mm[.] -> cyt_host_genenames_$EnsemblIDs
cyt_host_genenames_ <- cyt_host_genenames_[!is.na(cyt_host_genenames_$EnsemblIDs), ]
distinct(cyt_host_genenames_, EnsemblIDs, .keep_all = TRUE) -> cyt_host_genenames
print(paste(round(100 * nrow(cyt_host_genenames) / nrow(cyt_host), 2), "% of Mice genes kept after ensemblID translation and replicate removal"))


# TF activity
## PANcancer (Human)
data(dorothea_hs_pancancer, package = "dorothea")
regulons_hs <- dorothea_hs_pancancer %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))
minsize <- 5
ges.filter <- FALSE

cyt_tumor_genenames %>%
  column_to_rownames("EnsemblIDs") %>%
  dorothea::run_viper(regulons_hs,
    options = list(
      minsize = minsize, eset.filter = ges.filter,
      cores = 1, verbose = FALSE, nes = TRUE
    )
  ) -> tf_activities_stat_hs


## GTex (stroma)
data(dorothea_mm, package = "dorothea")
regulons_mm <- dorothea_mm %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))
minsize <- 5
ges.filter <- FALSE

cyt_host_genenames %>%
  column_to_rownames("EnsemblIDs") %>%
  dorothea::run_viper(regulons_mm,
    options = list(
      minsize = minsize, eset.filter = ges.filter,
      cores = 1, verbose = FALSE, nes = TRUE
    )
  ) -> tf_activities_stat_mm

# Writing
tf_activities_stat_hs %>%
  as_tibble(rownames = "TF") %>%
  write_tsv("02_Output/TFact_tumor_PanCan.tsv")

tf_activities_stat_mm %>%
  as_tibble(rownames = "TF") %>%
  write_tsv("02_Output/TFact_stroma_Gtex.tsv")
