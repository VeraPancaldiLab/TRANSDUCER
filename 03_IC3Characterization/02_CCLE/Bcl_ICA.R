library(tidyverse)
library(pdacmolgrad)
library(JADE)

setwd("/home/jacobo/Documents/03_Sauyeun_paper/02_CCLE/")

# Data loading
interesting_cl <- c("ASPC1_PANCREAS", "HUPT3_PANCREAS", "PANC1005_PANCREAS",
                    "CAPAN2_PANCREAS", "PATU8902_PANCREAS","SW1990_PANCREAS", "PATU8988S_PANCREAS", #this is an outlier
                    "DANG_PANCREAS", "PANC0327_PANCREAS", "HPAC_PANCREAS",
                    "HUPT4_PANCREAS", "PANC0203_PANCREAS","MIAPACA2_PANCREAS",
                    "PATU8988T_PANCREAS", "PANC1_PANCREAS") # No PL45_PANCREAS (absent)

plot_interesting <- gsub('_PANCREAS', '', interesting_cl)

## Metadata
metadata <- read_tsv("01_Input/metadata_Bcl_ST1.tsv")
colnames(metadata) <- c("cell_line", "Ser", "cancer_type",
                        "tissue", "gender", "ethnicity", "age")

metadata$cell_line <- str_replace_all(metadata$cell_line, "[^[:alnum:]]", "") %>%
  toupper()

metadata <- which(metadata$cell_line %in% plot_interesting) %>%
  metadata[., colnames(metadata)]

all(metadata$cell_line == plot_interesting[which(plot_interesting %in% 
                                                   metadata$cell_line)]) # THIS MUST BE TRUE SO YOU CAN REPLACE NAMES
metadata$cell_line <- interesting_cl[which(plot_interesting %in%
                                             metadata$cell_line)] 
## CCLE exp
### TPM
ccle_TPM_sym <- readRDS("02_Output/CCLE_symbols_TPM.RDS")
ccle_TPM_sym <- ccle_TPM_sym[,interesting_cl]
ccle_TPM_sym.non0 <- ccle_TPM_sym[rowSums(ccle_TPM_sym > 0) !=0, ]

### TMM
ccle_counts <- readRDS("02_Output/CCLE_symbols_counts.RDS") #gene symbols, read counts
ccle_counts <- ccle_counts[,interesting_cl]

#### normalization
ccle_counts_dge <- DGEList(ccle_counts)
ccle_TMM_dge <- calcNormFactors(ccle_counts_dge,
                                method = "TMM")

ccle_TMM <- cpm(ccle_TMM_dge)


### CHOOSE NORM
ccle_sym <- ccle_TPM_sym
ccle_sym.non0 <- ccle_TPM_sym.non0


# The ICA
## Functions
### Performn JADE ICA to get a list
###  of S  and A matrixes of a range of components
jade_mat_bestfit <- function(df, range.comp) {
  A_mats <- list()
  S_mats <- list()
  
  for (n.comp in range.comp) {
    jade_result <- JADE(df,
                        n.comp = n.comp
    )
    l_name <- paste(c("nc", n.comp), collapse = "")
      A_mats[[l_name]] <- jade_result[["A"]]
      suffix <- 1:ncol(jade_result[["A"]])
      colnames(A_mats[[l_name]]) <- paste("ic", suffix, sep = "")
    
      S_mats[[l_name]] <- jade_result[["S"]]
  }
  
  return(list("A" = A_mats, "S" = S_mats))
}

## Execution 
icas_list <- jade_mat_bestfit(ccle_sym.non0, 2:15)
icas_list$samples <- colnames(ccle_sym.non0)

## Exportation
saveRDS(object = icas_list, file = "02_Output/Bcl_bestfitICs.RDS")
