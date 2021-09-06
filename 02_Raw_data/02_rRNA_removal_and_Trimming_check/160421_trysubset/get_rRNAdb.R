# Title: get rRNA sequences from mice and hsa from ensembl 75 version
# Author details: Author: Jacobo Solorzano, Contact details: jacobo.solorzano@inserm.fr
# Script and data info: This script basically gets the rRNAs FASTAs to construct a rRNA db 
###################### to be used in rRNA read depletion
# INPUT: just to choose the mart, filters, and attributes from your licking
# OUTPUT: a fasta file "PDXrRNAs_ensembl75.fasta" with mmu and hsa ensembl 75 rRNAs
# Creation 09/04/21
# Modified -

library("biomaRt")
library("Biostrings")
setwd("~/Documents/02_TRANSDUCER/02_Raw_data/02_rRNA_removal_and_Trimming_check/")
# Choose the mart where you want your data to come
listMarts()

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   version = "75")

# Choose Dataset (Species)
listDatasets(mart)
hsa <- useDataset("hsapiens_gene_ensembl",mart=mart)
mmu <- useDataset("mmusculus_gene_ensembl",mart=mart)

# Choose Filters and attributes to be fetched
lfil <- listFilters(hsa)
filters <- c("biotype")
values <- c("rRNA")

latt <- listAttributes(hsa)
attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "cdna") 

# Fetch Results
hsa_rRNAs <- getBM(attributes=attributes, 
      filters = filters, 
      values = values, 
      mart = hsa)

mmu_rRNAs <- getBM(attributes=attributes, 
                   filters = filters, 
                   values = values, 
                   mart = mmu)

rRNAs <- rbind(hsa_rRNAs, mmu_rRNAs)
  
# Export as FASTA
seqs <- DNAStringSet(rRNAs[,"cdna"])

names(seqs) <- rRNAs$ensembl_transcript_id

writeXStringSet(seqs, "01_Input/PDXrRNAs_ensembl75.fasta")
