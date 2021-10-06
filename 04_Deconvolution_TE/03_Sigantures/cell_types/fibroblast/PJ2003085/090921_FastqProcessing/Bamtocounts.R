library(tidyverse)
library(Rsubread)
library(edgeR)

################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/04_Deconvolution_TE/03_Sigantures/cell_types/fibroblast/PJ2003085/090921_FastqProcessing")

# coding genes/all genes
featureType <- "CDS"
#featureType <- "exon"

bams <- list.files("02_Output/output_bams/", "*.bam", full.names = TRUE)
counts_etc <- featureCounts(files = bams,
                            GTF.featureType = featureType,
                            annot.ext = "01_Input/GRCh38.104/Homo_sapiens.GRCh38.104.gtf",
                            isGTFAnnotationFile = TRUE,
                            useMetaFeatures = TRUE, # to count gene IDs
                            isPairedEnd = TRUE,
                            nthreads = 14
              )

counts <- counts_etc$counts


colnames(counts) %>%
  str_remove("_Aligned.sortedByCoord.out.bam") -> 
  colnames(counts)


counts %>% as_tibble(rownames = "EnsemblID") %>%
  write_tsv(paste("02_Output/rawcounts_", featureType, ".tsv", sep = ""))

# Prefilter
counts_fil <- counts[which(rowSums(counts!=0) > 2 &
                           rowSums(select(data.frame(counts),ends_with("Pol")) > 0) &
                           rowSums(select(data.frame(counts),ends_with("tot")) > 0)),]

counts_fil %>% as_tibble(rownames = "EnsemblID") %>%
  write_tsv(paste("02_Output/filteredcounts_", featureType, ".tsv", sep = ""))

# Normalization
counts_es <- DGEList(counts_fil)
counts_tmm <- calcNormFactors(counts_es,
                               method = "TMM")
counts_tmm <- cpm(counts_tmm)
counts_tmm %>% as_tibble(rownames = "EnsemblID") %>%
  write_tsv(paste("02_Output/tmmcounts_", featureType, ".tsv", sep = ""))

