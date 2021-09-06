setwd("/home/jacobo/Documents/phd_project_perse/preprocessing/signatures/cell_state/fibrotic_ecm/")
library(tidyverse)
library(GEOquery)
library(lumi)

# data loading
# ## no norm
# no_norm <- data.frame(read_tsv("data/GSE45686_non_normalized.txt"))
# rownames(no_norm) <- no_norm$ID_REF
# no_norm <- subset(no_norm, select = -c(ID_REF))
# colnames(no_norm) <- paste(rep(c("SAMPLE", "PVAL"),40), sort(rep(1:40, 2)))
# 
# signals <- no_norm[grepl("SAMPLE", colnames(no_norm))]
# significances <- no_norm[grepl("PVAL", colnames(no_norm))]
# 
# ## metadata
# metadata <- gse[["GSE45686_series_matrix.txt.gz"]]@phenoData@data
# 
# metadata <- metadata[c("description", "title", "characteristics_ch1.2", "characteristics_ch1.1", "characteristics_ch1.4", "batch:ch1")]
# colnames(metadata) <- str_replace(colnames(metadata), 'characteristics_ch1.2', 'cell_type')
# colnames(metadata) <- str_replace(colnames(metadata), 'characteristics_ch1.1', 'ecm_type')
# colnames(metadata) <- str_replace(colnames(metadata), 'characteristics_ch1.4', 'fraction')
# colnames(metadata) <- str_replace(colnames(metadata), 'batch:ch1', 'batch')
# 
# 
# metadata$cell_type <- str_replace(metadata$cell_type, 'cell line: ', '')
# metadata$ecm_type <- str_replace(metadata$ecm_type, "ecm type: ", '')
# metadata$fraction <- str_replace(metadata$fraction, "molecule subtype: ", '')

## normalized
gse <- getGEO("GSE45686", GSEMatrix=FALSE)
samples <- GSMList(gse)
detectionP <- data.frame("ID_REF" = Table(samples[[1]])["ID_REF"])
exprV <- data.frame("ID_REF" = Table(samples[[1]])["ID_REF"])
metadata <- data.frame()
for (s in samples) {
  if (all(exprV$ID_REF != Table(s)["ID_REF"])) {
    stop("not matching refference IDs!")
  }
  exprV[, Meta(s)$description] = Table(s)["VALUE"]
  detectionP[, Meta(s)$description] = Table(s)["Detection Pval"]
  metadata[Meta(s)$description, c("cell_type","ecm_type","cell_line","batch","fraction")] <- Meta(s)$characteristics_ch1
}
rownames(exprV) <- exprV$ID_REF
rownames(detectionP) <- detectionP$ID_REF
exprV$ID_REF <- NULL
detectionP$ID_REF <- NULL

## metadata
metadata$cell_type <- str_replace(metadata$cell_type, 'cell type: ', '')
metadata$ecm_type <- str_replace(metadata$ecm_type, "ecm type: ", '')
metadata$cell_line <- str_replace(metadata$cell_line, 'cell line: ', '')
metadata$batch <- str_replace(metadata$batch, 'batch: ', '')
metadata$fraction <- str_replace(metadata$fraction, "molecule subtype: ", '')

## separation of fractions
poly.v <- metadata[metadata$fraction == "polysome-associated RNA",] %>% row.names()
total.v <- metadata[metadata$fraction == "total RNA",] %>% row.names()

# Initial distribution
boxplot(exprV[sample(rownames(exprV), 10000),], #set poly.v or total.v to only see separately
        main = "sample vs sample. 10000 random genes", xaxt = "n")

boxplot(t(exprV[sample(rownames(exprV), 100),]),
        main = "gene vs gene 100 random genes", xaxt = "n")



