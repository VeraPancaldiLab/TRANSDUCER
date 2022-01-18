library(tidyverse)
library(edgeR)
library(mMCPcounter)
################################################################################
setwd("~/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/180122_Various/")

# Data loading

read_tsv("../../00_Data/Processed_data/rawHost_Cyt.tsv") -> cyt_raw

cleanHost <- c(
  "PDAC018T",
  "PDAC029T",
  "PDAC013T",
  "PDAC020T",
  "PDAC021T",
  "PDAC017T",
  "PDAC031T",
  "PDAC012T",
  "PDAC032T",
  "PDAC028T",
  "PDAC030T",
  "PDAC027T",
  "PDAC025T",
  "PDAC011T",
  "PDAC024T",
  "PDAC008T",
  "PDAC026T",
  "PDAC019T",
  "PDAC015T",
  "PDAC022T"
) # 20

# Preprocessing
cyt_raw %>% column_to_rownames("EnsemblID") %>%
  dplyr::select(cleanHost) %>%
  dplyr::filter_all(any_vars(. != 0)) %>%
  edgeR::DGEList() %>% edgeR::calcNormFactors(method = "TMM") %>%
  edgeR::cpm() -> cyt_norm # log transformation affects deconvolution performance, so this is non logged

ensembl75 <- useEnsembl(biomart = "genes",
                          dataset = "mmusculus_gene_ensembl",
                          version = 75)


annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id',
                                        'entrezgene'), mart = ensembl75)

# Deconvolution
## mMCPcounter
?mMCPcounter::mMCPcounter.estimate # runnable with EnsemblIDs
mMCPcounter_res <- mMCPcounter.estimate(exp = cyt_norm,
                                               features = "ENSEMBL.ID")
mMCPcounter_res %>% rownames_to_column("cell_types") %>% write_tsv("02_Output/mMCPcounter_results.tsv")

## ImmuCC (CYBERSORT 4 mice)
download.file(url = "https://raw.githubusercontent.com/wuaipinglab/ImmuCC/master/webserver/SignatureMatrix.rnaseq.csv",
              destfile = "01_Input/SignatureMatrix.rnaseq.csv")

system("sed 's/,/\t/g' 01_Input/SignatureMatrix.rnaseq.csv > 01_Input/SignatureMatrix.ImmuCC.tsv")

boxplot(cyt_norm)

cyt_norm %>% as_tibble(rownames = "EnsemblIDs") %>% write_tsv("01_Input/HostCyt_norm.tsv")

### run online and fetch results now
read_tsv("02_Output/CIBERSORT.Output_Job13.txt") %>% as_tibble(.) %>%
  column_to_rownames("Input Sample") %>%
  dplyr::select(!c("P-value", "Pearson Correlation", "RMSE")) -> ImmuCC_res

ImmuCC_res %>% rownames_to_column("samples") %>% write_tsv("02_Output/ImmuCC_results.tsv")
