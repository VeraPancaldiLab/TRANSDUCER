library(tidyverse)
library(biomaRt)

setwd("/home/jacobo/Documents/03_Sauyeun_paper/02_CCLE/")

# Raw counts -> TMM
# normalization should be performed 
# In the context of the specific to analyze, so should be done once samples are substarcted 
r.cl <- read.delim(file="01_Input/CCLE_RNAseq_genes_counts_20180929.gct", skip=2)

##Dup gene merging (sum as these are reads)
r.cl <- aggregate(r.cl[c(-1,-2)],
                    list(r.cl$Description),
                    FUN = sum)

rownames(r.cl) <- r.cl[,"Group.1"]
r.cl[,"Group.1"] <- NULL

saveRDS(r.cl, "02_Output/CCLE_symbols_counts.RDS")


## Data loading
tpm.cl <- read_tsv("01_Input/CCLE_RNAseq_rsem_genes_tpm_20180929.txt")

### Ensembl
cl_ens <- as.data.frame(tpm.cl[, -c(1,2)])
ensembl_ids <- tpm.cl$gene_id
ensembl_ids <- sub('\\.[0-9]*$', '', ensembl_ids) #remove version of enseml.ids
rownames(cl_ens) <- ensembl_ids

### Gene Symbol
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   host="grch37.ensembl.org", path="/biomart/martservice", #Yeah, its the grch37 assembly
                   dataset = "hsapiens_gene_ensembl")

ensembl_table <- getBM(attributes = c('ensembl_gene_id', "hgnc_symbol", 'external_gene_name'),
                       mart = ensembl)

cl_sym <- merge(cl_ens, ensembl_table, by.x="row.names", by.y="ensembl_gene_id",
                all.x = T)

print(paste(nrow(cl_sym)/nrow(cl_ens)*100, "% rows kept after merging for gene symbol"))

cl_sym <- aggregate(cl_sym,
                    list(cl_sym$external_gene_name),
                    FUN = mean)
print(paste(nrow(cl_sym)/nrow(cl_ens)*100, "% rows kept after aggregating for duplicated external gene names"))

row.names(cl_sym) <- cl_sym[,"Group.1"]
cl_sym[,c("Group.1","Row.names", "external_gene_name", "hgnc_symbol")] <- NULL

## Data export
saveRDS(cl_sym, "02_Output/CCLE_symbols_TPM.RDS")
saveRDS(cl_ens, "02_Output/CCLE_ensembl_TPM.RDS")