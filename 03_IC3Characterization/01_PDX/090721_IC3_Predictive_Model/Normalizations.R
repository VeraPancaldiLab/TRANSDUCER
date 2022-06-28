setwd("~/Documents/02_TRANSDUCER/03_IC3Characterization/01_PDX/090721_IC3_Predictive_Model/")
source("~/Documents/01_General_tools/00_Other_peoples/RNAseq/counts_to_tpm/counts_to_tpm.R")
library(tidyverse)
library(edgeR)

# load raw counts
sauyeun_rc <- read_tsv("01_Input/geneCount_raw_28s_totalRNA.tsv") %>% column_to_rownames("EnsemblID")
pacaomics_rc <- read_tsv("01_Input/Human-Tumor_rawcount_Transcriptome.tsv") %>% column_to_rownames("EnsemblID")

# Generate gene lengths from gtf (http://genomespot.blogspot.com/2019/01/using-gtf-tools-to-get-gene-lengths.html)
get_GRCh37_genelegths <- "python3 ~/Documents/01_General_tools/00_Other_peoples/RNAseq/gtf_to_genelengths/gtftools.py -l 01_Input/GRCh37.genelengths 01_Input/Homo_sapiens.GRCh37.87.gtf"
get_GRCh38_genelegths <- "python3 ~/Documents/01_General_tools/00_Other_peoples/RNAseq/gtf_to_genelengths/gtftools.py -l 01_Input/GRCh38.genelengths 01_Input/Homo_sapiens.GRCh38.104.gtf"
system(get_GRCh37_genelegths)
system(get_GRCh38_genelegths)
GRCh37_genelegths <- read_tsv("01_Input/GRCh37.genelengths") %>% column_to_rownames("gene")
GRCh38_genelegths <- read_tsv("01_Input/GRCh38.genelengths") %>% column_to_rownames("gene")


# Normalize
## TPM (fix)
sauyeun_tpm <- simple_tpm(sauyeun_rc, lengths = GRCh38_genelegths[rownames(sauyeun_rc), "mean", drop=F])
pacaomics_tpm <- simple_tpm(pacaomics_rc, lengths = GRCh38_genelegths[rownames(pacaomics_rc), "mean"])

## TMM
### 
sauyeun_es <- DGEList(sauyeun_rc)
sauyeun_tmm <- calcNormFactors(sauyeun_es,
                                method = "TMM")
sauyeun_tmm <- cpm(sauyeun_tmm)
sauyeun_tmm %>% as_tibble(rownames = "EnsemblID") %>% write_tsv("02_Output/sauyeun_tmm.tsv")
###
pacaomics_es <- DGEList(pacaomics_rc)
pacaomics_tmm <- calcNormFactors(pacaomics_es,
                               method = "TMM")
pacaomics_tmm <- cpm(pacaomics_tmm)
pacaomics_tmm %>% as_tibble(rownames = "EnsemblID") %>% write_tsv("02_Output/pacaomics_tmm.tsv")

