library(tidyverse)
library(Hmisc)
library(corrplot)
library(lmtest)
library(ggridges)
library(Biobase)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/081221_TranslationEfficacy/")

# Data loading
normHost_cyt <- read_tsv("../../00_Data/Processed_data/normHost_Cyt.tsv")
normHost_pol <- read_tsv("../../00_Data/Processed_data/normHost_Pol.tsv")
sample_info <- read_tsv("../../00_Data/Processed_data/sample_info.tsv") %>% column_to_rownames("sample")

# Joint dataframe creation
sn_cyt <- paste(colnames(normHost_cyt[-1]), "cyt", sep = "_")
sn_pol <- paste(colnames(normHost_pol[-1]), "pol", sep = "_")

normHost_cyt.tmp <- normHost_cyt
normHost_pol.tmp <- normHost_pol

colnames(normHost_cyt.tmp) <- c("EnsemblID",sn_cyt)
colnames(normHost_pol.tmp) <- c("EnsemblID",sn_pol)

norm_host <- inner_join(normHost_cyt.tmp, normHost_pol.tmp, by= "EnsemblID")

# eset creation
## PhenoData
sn <- c(colnames(normHost_cyt[-1]), colnames(normHost_pol[-1]))

frac <- c(rep("cyt", length(sn_cyt)),
          rep("pol", length(sn_pol)))

rna_conc <- c(sample_info[colnames(normHost_cyt[-1]), "RNAconc_cyt"],
              sample_info[colnames(normHost_pol[-1]), "RNAconc_pol"])

tumor_ica <- bind_rows(sample_info[colnames(normHost_cyt[-1]),str_subset(colnames(sample_info), "ICA")],
                       sample_info[colnames(normHost_pol[-1]),str_subset(colnames(sample_info), "ICA")],)

rownames(tumor_ica) <- c(sn_cyt, sn_pol)
  
phenotype.tmp <- data.frame(sample = factor(sn),
                        fraction = factor(frac),
                        rna_conc = as.double(rna_conc))

rownames(phenotype.tmp) <- colnames(norm_host)[-1]

phenotype <- merge(phenotype.tmp, tumor_ica, by = "row.names",) %>% column_to_rownames("Row.names")

phenoData <- new("AnnotatedDataFrame", data = phenotype)
norm_host %>% column_to_rownames("EnsemblID") %>% data.matrix() -> assayData

dset <- ExpressionSet(assayData = assayData[,rownames(phenotype)], phenoData = phenoData)
