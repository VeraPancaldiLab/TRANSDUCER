library(tidyverse)
library(biomaRt)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/100122_ICABoot/")
source("functions.R")
range.comp <- 2:20

# load TEs
TEs <- read_tsv("../081221_TranslationEfficacy/02_Output/TEs.tsv")

# load annotation with Biomart
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "mmusculus_gene_ensembl",
                        version = 75)

#listAttributes(ensembl75, page="feature_page")
annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                          'external_gene_id',
                          'entrezgene',
                          'chromosome_name'), mart = ensembl75)

# filtering
## XY
annot_ensembl75 %>% dplyr::filter(!chromosome_name %in% c("X", "Y")) %>%
  pull(ensembl_gene_id) -> non_sex

TEs %>% dplyr::filter(EnsemblID %in% non_sex) -> TEs_

## mean centring gene wise
TEs_ %>% column_to_rownames("EnsemblID") %>%
  apply(1, function(x) x - mean(x)) %>% t() %>%
  data.frame() -> TEs__

## Inter quartile range (measure variability)
TEs__ %>% apply(1, IQR) -> iqrs
mostvar <- iqrs[iqrs > median(iqrs)]
TEs__ %>% rownames_to_column("EnsemblID") %>%
  dplyr::filter(EnsemblID %in% names(mostvar)) %>%
  column_to_rownames("EnsemblID") -> TEs_icaready

# ICA
## Baseline before bootstrap
TEs_icaready %>% jade_range(range.comp, MARGIN = 1) -> base_res_gene
TEs_icaready %>% jade_range(range.comp, MARGIN = 2) -> base_res_sample

## Bootstrap
gene_boot <- jade_choosencom(TEs_icaready, base_res_gene,
                              MARGIN = 1,
                              iterations = 10,
                              seed = 0
)

sample_boot <- jade_choosencom(TEs_icaready, base_res_sample,
                             MARGIN = 2,
                             iterations = 10,
                             seed = 0
)

boot_plots(s_boot = sample_boot, g_boot = gene_boot, name = "TE")
