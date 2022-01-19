library(tidyverse)
library(mMCPcounter)
################################################################################
setwd("~/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/180122_Various/")

# Data loading
read_tsv("01_Input/HostCyt_foranalysis.tsv") %>% column_to_rownames("EnsemblIDs") -> cyt_host

# mMCPcounter
?mMCPcounter::mMCPcounter.estimate # runnable with EnsemblIDs
mMCPcounter_res <- mMCPcounter.estimate(exp = cyt_host,
                                               features = "ENSEMBL.ID")
mMCPcounter_res %>% rownames_to_column("cell_types") %>% write_tsv("02_Output/mMCPcounter_results.tsv")

# ImmuCC (CYBERSORT 4 mice)
download.file(url = "https://raw.githubusercontent.com/wuaipinglab/ImmuCC/master/webserver/SignatureMatrix.rnaseq.csv",
              destfile = "01_Input/SignatureMatrix.rnaseq.csv")

system("sed 's/,/\t/g' 01_Input/SignatureMatrix.rnaseq.csv > 01_Input/SignatureMatrix.ImmuCC.tsv")


cyt_host %>% as_tibble(rownames = "EnsemblIDs") %>% write_tsv("01_Input/HostCyt_norm.tsv")

## run online with these files and frtch results  
read_tsv("02_Output/CIBERSORT.Output_Job13.txt") %>% as_tibble(.) %>%
  column_to_rownames("Input Sample") %>%
  dplyr::select(!c("P-value", "Pearson Correlation", "RMSE")) -> ImmuCC_res

ImmuCC_res %>% rownames_to_column("samples") %>% write_tsv("02_Output/ImmuCC_results.tsv")
