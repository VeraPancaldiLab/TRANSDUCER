setwd("~/Documents/03_Sauyeun_paper/01_PDX/090721_IC3_Predictive_Model/")
source("~/Documents/01_General_tools/00_Other_peoples/RNAseq/counts_to_tpm/counts_to_tpm.R")

# load raw counts
rc_sauyeun <- read_tsv("01_Input/geneCount_raw_28s_totalRNA.tsv") %>% column_to_rownames("EnsemblID")
rc_pacaomics <- read_tsv("01_Input/Human-Tumor_rawcount_Transcriptome.tsv") %>% column_to_rownames("EnsemblID")

# Generate gene lengths from gtf (http://genomespot.blogspot.com/2019/01/using-gtf-tools-to-get-gene-lengths.html)
get_genelegths <- "python3 ~/Documents/01_General_tools/00_Other_peoples/RNAseq/gtf_to_genelengths/gtftools.py -l 01_Input/GRCh37.genelengths 01_Input/Homo_sapiens.GRCh37.87.gtf"
system(get_genelegths)
gene_lengths <- read_tsv("01_Input/GRCh37.genelengths") %>% column_to_rownames("gene")

# Normalize
tpm_sauyeun <- counts_to_tpm(rc_sauyeun, featureLength = gene_lengths[rownames(rc_sauyeun), "mean"], meanFragmentLength = ??????)
