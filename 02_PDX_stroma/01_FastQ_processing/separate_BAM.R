#!/usr/bin/env Rscript
wd <- "/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/01_FastQ_processing"
setwd(wd)
# Try with the version of 2017 rsubread
#install.packages("https://bioconductor.org/packages/3.5/bioc/src/contrib/Rsubread_1.26.1.tar.gz", repos=NULL, type="source") 

#library(devtools)
#install_github("cit-bioinfo/SMAP")
source("01_Input/SMAP/R/SMAPcount.R")
#library("SMAP")
library(tidyverse)

# Split Bams into mice and human counts
counts = SMAPcount(SMAPBAM = "01_Input/output_bams",
                   GTF = "02_Output/combined_reference/combined.gtf")
countsTumor = counts$FCcountsTumor
countsHost = counts$FCcountsHost
save.image("02_Output/Backup.SMAP.RData")
load("02_Output/Backup.SMAP.RData")

# Split counts into Total and Polysome
## Edit sample names 
colnames(countsTumor) %>% str_remove("X01_Input.output_bams.") %>% 
  str_remove("_Aligned.sortedByCoord.out.bam") -> colnames(countsTumor)

colnames(countsHost) %>% str_remove("X01_Input.output_bams.") %>% 
  str_remove("_Aligned.sortedByCoord.out.bam") -> colnames(countsHost)

## create object with Total/Polysome sample names
all(colnames(countsTumor) %in% colnames(countsHost)) %>% stopifnot()

colnames(countsTumor) %>% str_detect("_F8.9", negate = T) %>% 
  colnames(countsTumor)[.] -> names_total

colnames(countsTumor) %>% str_detect("_F8.9") %>%
  colnames(countsTumor)[.] -> names_polysome

all(str_remove(names_polysome, "_F8.9") == names_total) %>% stopifnot()
sample_names = tibble(total = names_total,
                      polysome = names_polysome)

## Split Tumour
countsTumor %>% as_tibble(rownames = "EnsemblID") %>%
  dplyr::select(c("EnsemblID", sample_names$total)) %>%
  write_tsv("02_Output/countsTumor_total.tsv")

countsTumor %>% as_tibble(rownames = "EnsemblID") %>%
  dplyr::select(c("EnsemblID", sample_names$polysome)) %>%
  rename_at(sample_names$polysome, ~ sample_names$total) %>% # this renames with official sample names
  write_tsv("02_Output/countsTumor_polysome.tsv") 

## Split Stroma
countsHost %>% as_tibble(rownames = "EnsemblID") %>%
  dplyr::select(c("EnsemblID", sample_names$total)) %>%
  write_tsv("02_Output/countsHost_total.tsv")  

countsHost %>% as_tibble(rownames = "EnsemblID") %>%
  dplyr::select(c("EnsemblID", sample_names$polysome)) %>%
  rename_at(sample_names$polysome, ~ sample_names$total) %>% # this renames with official sample names
  write_tsv("02_Output/countsHost_polysome.tsv")

## Example data
# d = system.file("extdata","Example_STAR_OUTPUT", package = "SMAP")
# dcounts = SMAPcount(SMAPBAM=d, GTF = "/home/jacobo/Documents/02_TRANSDUCER/03_Deconvolution/01_FastQ_processing/02_Output/combined_reference/combined.gtf")
# countsTumor = counts$FCcountsTumor
# countsHost = counts$FCcountsHost

# # Normalizations
# ## TMM
# library(edgeR)
# countsTumor.edgeR <- DGEList(countsTumor)
# countsTumor.edgeR <- calcNormFactors(countsTumor.edgeR,
#                                      method = "TMM")
# tail(cpm(countsTumor.edgeR))
# tail(cpm(countsTumor))
# tail(countsTumor/countsTumor.edgeR$samples$norm.factors)
# 
# countsTumor.tmm <- cpm(countsTumor.edgeR)
# countsTumor.tmm.l <- cpm(countsTumor.edgeR, log = TRUE) # calculate mean dif between these two
# 
# ## Upper Quartile
# library(edgeR)
# countsTumor.edgeR <- DGEList(countsTumor)
# countsTumor.edgeR <- calcNormFactors(countsTumor.edgeR,
#                                      method = "upperquartile")
# 
# countsTumor.uq <- cpm(countsTumor.edgeR)
# countsTumor.uq.l <- cpm(countsTumor.edgeR, log = TRUE)
# 
# ## TPM
# #https://support.bioconductor.org/p/91218/
# tpm <- function(counts,len) {
#   x <- counts/len
#   return(t(t(x)*1e6/colSums(x)))
# }
