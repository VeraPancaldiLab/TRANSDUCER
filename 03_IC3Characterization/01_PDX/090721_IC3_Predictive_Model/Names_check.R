setwd("Documents/03_Sauyeun_paper/01_PDX/090721_IC3_Predictive_Model/")
library(tidyverse)
library(edgeR)
library(Hmisc)
library(corrplot)

# data load
Hcpmallrna_TMM <-  readRDS("../Hcpmallrna.rds") # Remy sent it to me explicitly to check (TMM)
raw_counts <- read_tsv("01_Input/geneCount_raw_28s_totalRNA.tsv") %>% column_to_rownames("EnsemblID") %>% .[rownames(Hcpmallrna_TMM),]

# Raw count norm
dge <- DGEList(raw_counts)
dge <- calcNormFactors(dge,
                        method = "TMM")
dge_TMM <- cpm(dge)

# sample name conversion (comment to see unsup)
equivalences <- read_tsv("../Remy_processed_data/sample_names_equivalences.tsv") %>% .[,c("fastq_name", "CITID")]
equivalences$CITID <- paste(equivalences$CITID, "check", sep = "_")
translate <- deframe(equivalences)

dge_TMM %>% colnames(.) %>% translate[.] -> colnames(dge_TMM)


# correlation
corr <- rcorr(dge_TMM, Hcpmallrna_TMM, "spearman")
corr$r <- corr$r[colnames(dge_TMM), colnames(Hcpmallrna_TMM)]
corr$P <- corr$P[colnames(dge_TMM), colnames(Hcpmallrna_TMM)]
corrplot(corr$r, order="original", , is.corr = FALSE, type = "full", p.mat = corr$P,
         sig.level = 0.05, method="color", insig = "label_sig")
