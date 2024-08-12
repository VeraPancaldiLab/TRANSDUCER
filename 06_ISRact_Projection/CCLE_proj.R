# ISRactPCA projection into Broad Institute CCLE cohort
################################################################################
#' Import different files
#' Filter data and project ISRactPCA
#' Produce projection scatterplot
#' Compare ISRactPCA with PAMG score
#' PHGDH/CBS doublecheck
################################################################################

# Import libraries
library(tidyverse)
library(biomaRt)
library(pdacmolgrad) # devtools::install_github("RemyNicolle/pdacmolgrad")
library(edgeR)
library(Hmisc)
library(ggrepel)
################################################################################
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")

# Import datasets
## CCLE
### load raw data
ccle_raw <- read_delim(file = "data/CCLE/CCLE_RNAseq_genes_counts_20180929.gct", skip = 2) %>%
  dplyr::rename(EnsemblID = Name, GeneName = Description) %>%
  dplyr::mutate(EnsemblID = str_remove(EnsemblID, "\\.[0-9]*$"))

ccle_info <- read_csv("data/CCLE/primary-screen-cell-line-info.csv") %>%
  dplyr::mutate(ISRact_bin = if_else(ccle_name %in% c("ASPC1_PANCREAS", "PATU8902_PANCREAS", "PATU8988T_PANCREAS", "MIAPACA2_PANCREAS"),
    if_else(ccle_name %in% c("ASPC1_PANCREAS", "PATU8902_PANCREAS"),
      "high_ISRact", "low_ISRact"
    ),
    NA
  ))

## ISRactPCA projection
pca_pdx <- read_rds("data/Classifiers/pca_pdx_ENZO.RDS")
### inherited sample info and top samples from Shin PDX training data
sample_info <- read_delim("data/Sauyeun_PDX/sample_info.tsv",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)

top_samples <- arrange(sample_info, ICA3) %>%
  dplyr::slice(unique(c(1:5, n() - 0:4))) %>%
  mutate(ISRact_bin = ifelse(ICA3 < 0, "low_ISRact", "high_ISRact")) %>%
  dplyr::rename(ISRact = "ICA3") %>%
  arrange(sample)

################################################################################
# PARAMETERS
keep_CCLE <- "PANCREAS" # ALL | PANCREAS
norm_method <- "upperquartile" # TMM | upperquartile
################################################################################

# Processing and Normalization
## CCLE
if (keep_CCLE == "PANCREAS") {
  keep <- dplyr::filter(ccle_info, primary_tissue == "pancreas") %>%
    dplyr::select(ccle_name) %>%
    deframe()
  ccle_raw <- dplyr::select(ccle_raw, EnsemblID, GeneName, all_of(keep))
}
ccle_norm_ <- dplyr::select(ccle_raw, -GeneName) %>%
  column_to_rownames("EnsemblID") %>%
  DGEList() %>%
  calcNormFactors(method = norm_method) %>%
  cpm(log = TRUE)

ccle_norm <- ccle_norm_ %>%
  t() %>%
  as_tibble(rownames = "ccle_name")

### create version with Gene Names
ccle_translate <- deframe(ccle_raw[c("EnsemblID", "GeneName")])

ccle_norm_gn <- as_tibble(ccle_norm_, rownames = "EnsemblID") %>%
  mutate(GeneName = ccle_translate[EnsemblID]) %>%
  dplyr::select(-EnsemblID) %>%
  group_by(GeneName) %>%
  summarise_all(sum)

### Add PAMG status
type_pamg <- projectMolGrad(newexp = column_to_rownames(ccle_norm_gn, "GeneName"), geneSymbols = ccle_norm_gn$GeneName) %>%
  as_tibble(rownames = "ccle_name")

ccle_info <- dplyr::select(type_pamg, ccle_name, PDX) %>%
  right_join(ccle_info, "ccle_name") %>%
  dplyr::rename(PAMG = PDX)


# Projection
## CCLE
### Add missing genes to expression before the PCA
filter_pca <- function(.data, objective) {
  as_tibble(.data, rownames = "tmp") %>%
    dplyr::filter(tmp %in% objective) %>%
    column_to_rownames("tmp") %>%
    data.matrix()
}

pca_pdx$rotation <- filter_pca(pca_pdx$rotation, names(ccle_norm)[-1])
pca_pdx$center <- filter_pca(pca_pdx$center, names(ccle_norm)[-1])
pca_pdx$scale <- filter_pca(pca_pdx$scale, names(ccle_norm)[-1])

ccle_PCAspace <- predict(pca_pdx, ccle_norm) %>%
  as_tibble() %>%
  mutate(ccle_name = ccle_norm$ccle_name, .before = 1) %>%
  left_join(ccle_info[, c("ccle_name", "primary_tissue", "ISRact_bin", "PAMG")], by = "ccle_name") %>%
  dplyr::rename(ISRactPCA = "PC1")

ccle_ISRact_projection <- arrange(ccle_PCAspace, ISRactPCA) %>%
  mutate(ISRactPCA_bin = cut(.$ISRactPCA, breaks = c(quantile(.$ISRactPCA, c(0:3 / 3))), labels = c("low_ISRactPCA", "intermediate_ISRactPCA", "high_ISRactPCA"), include.lowest = TRUE)) %>%
  dplyr::select(ccle_name, ISRactPCA, ISRactPCA_bin)


# Plot pca and projections
## Add ISR status to PCA df
pca_full_df <- pca_pdx[["x"]] %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  dplyr::rename(ISRactPCA = "PC1") %>%
  inner_join(top_samples[, c("sample", "ISRact_bin")])

#-------------------------------------------------------------------------------
# Plot projecion

show_projection <- bind_rows(
  as_tibble(pca_full_df) %>% mutate(dataset = "Shin et al. PDX"),
  as_tibble(ccle_PCAspace %>% mutate(
    ISRact_bin = if_else(is.na(ISRact_bin), "Unknown", ISRact_bin),
    dataset = "CCLE"
  ))
)

projection_scatter_ccle <- dplyr::mutate(show_projection, ISRact_bin = str_replace(ISRact_bin, "ICA3", "ISRact")) %>%
  dplyr::filter(
    dataset %in% c("Shin et al. PDX", "PACAOMICS PDX", "CPTAC", "CCLE"), # "Shin et al. PDX","PACAOMICS PDX", "CPTAC", "CCLE"
    ISRact_bin %in% c("low_ISRact", "high_ISRact", "intermediate_ISRact", "Unknown")
  ) %>%
  ggplot(aes(x = ISRactPCA, y = PC2, color = ISRact_bin, shape = dataset)) +
  geom_point() +
  scale_shape_manual(values = c(`Shin et al. PDX` = 16, `PaCaOmics PDX` = 17, CPTAC = 15, CCLE = 18)) +
  scale_color_manual(values = c(low_ISRact = "seagreen", high_ISRact = "tomato3", intermediate_ISRact = "grey", Unknown = "#619CFF")) +
  theme_bw()

ggsave(projection_scatter_ccle,
  filename = "results/Figures/projection_scatter_ccle.svg",
  width = 7,
  height = 3
)
#-------------------------------------------------------------------------------

# Plot comparisons with Basal/Classical and ISRact
## PAMG vs PC1
correlation_plotter(data = ccle_PCAspace, col1 = "PAMG", col2 = "ISRactPCA", data_name = "CCLE")

## Coloured version for thesis
corr_spearman <- rcorr(ccle_toplot[["ISRactPCA"]], ccle_toplot[["PAMG"]], type = "spearman")
stats <- paste0("Spearman: R = ", round(corr_spearman$r["x", "y"], 2), ", pval = ", round(corr_spearman$P["x", "y"], 4))

scatter_isractpcapamg_ccle <- ggplot(ccle_toplot, aes(x = ISRactPCA, y = PAMG, label = str_remove(ccle_name, "_PANCREAS"))) +
  geom_point(aes(color = ISRact_bin)) +
  scale_color_manual(values = c(low_ISRact = "seagreen", high_ISRact = "tomato3", intermediate_ISRact = "grey", Unknown = "#619CFF")) +
  geom_smooth(method = lm) +
  labs(
    title = paste0("Comparison between ISRactPCA and PAMG in CCLE"),
    subtitle = stats
  ) +
  theme_bw()

ggsave(scatter_isractpcapamg_ccle,
  filename = "results/Figures/scatter_isractpcapamg_ccle.svg",
  width = 5,
  height = 4
)
#-------------------------------------------------------------------------------
# CCLE plot with names and colored according to CBS and PHGDH gene expression
ccle_toplot <- pivot_longer(ccle_norm_gn, -GeneName, names_to = "ccle_name", values_to = "expression") %>%
  pivot_wider(id_cols = "ccle_name", names_from = "GeneName", values_from = "expression") %>%
  inner_join(ccle_PCAspace, by = "ccle_name")

corr_spearman <- rcorr(ccle_toplot[["PHGDH"]], ccle_toplot[["CBS"]], type = "spearman")
stats <- paste0("Spearman: R = ", round(corr_spearman$r["x", "y"], 2), ", pval = ", round(corr_spearman$P["x", "y"], 4))

scatter_phgdhcbs_ccle <- ggplot(ccle_toplot, aes(x = PHGDH, y = CBS, label = str_remove(ccle_name, "_PANCREAS"))) +
  geom_point(aes(color = ISRact_bin)) +
  scale_color_manual(values = c(low_ISRact = "seagreen", high_ISRact = "tomato3", intermediate_ISRact = "grey", Unknown = "#619CFF")) +
  geom_smooth(method = lm) +
  labs(
    title = paste0("Comparison between PHGDH and CBS in CCLE"),
    subtitle = stats
  ) +
  theme_bw()

ggsave(scatter_phgdhcbs_ccle,
  filename = "results/Figures/scatter_phgdhcbs_ccle.svg",
  width = 5,
  height = 4
)
#-------------------------------------------------------------------------------
