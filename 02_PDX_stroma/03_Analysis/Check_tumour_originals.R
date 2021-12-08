library(tidyverse)
library(Hmisc)
library(corrplot)
###########################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis")
mine <- read_tsv("../00_Data/Processed_data/normHost_Cyt.tsv")
load("../00_Data/Remy_processed_data/all_RNA/Stroma/Mcpmallrna.RData")
colnames(Mcpmallrna) <- paste(colnames(Mcpmallrna), "Remy", sep = "_")

Mcpmallrna %>%
  as_tibble(rownames = "EnsemblID") %>%
  dplyr::filter(EnsemblID %in% mine$EnsemblID) -> Mcpmallrna


all(Mcpmallrna$EnsemblID == mine$EnsemblID) %>% stopifnot()


corr <- rcorr(mine %>% column_to_rownames("EnsemblID") %>% data.matrix(), Mcpmallrna %>% column_to_rownames("EnsemblID") %>% data.matrix(), "pearson")
corr$r <- corr$r[colnames(mine)[-1], colnames(Mcpmallrna)[-1]]
corr$P <- corr$P[colnames(mine)[-1], colnames(Mcpmallrna)[-1]]
corrplot(corr$r, order="original", , is.corr = FALSE, type = "full", p.mat = corr$P,
         sig.level = 0.05, method="color", insig = "label_sig")
