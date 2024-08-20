# ISRactPCA generation and ISRactPCA/ICA gene weight analysis

################################################################################
#' Import different files
#' Filter data to produce ISRactPCA and export it
#' Predict ISRActPCA on the full cohort
#' Produce several plots fort thesis figures and beyond
################################################################################

# Import libraries
library(tidyverse)
library(biomaRt)
library(edgeR)
library(factoextra)
library(ggpubr)
library(Hmisc)
library(ggVennDiagram)
library(pdacmolgrad) # devtools::install_github("RemyNicolle/pdacmolgrad")
library(readxl)

## setwd and import local functions
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")
source("src/mostCorrGenes.R")

# Import datasets
## Sauyeun PDX
### metadata
sample_info <- read_delim("data/Sauyeun_PDX/sample_info.tsv",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)

sample_names_equivalences <- read_delim("data/Sauyeun_PDX/sample_names_equivalences.tsv",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)

### raw count processed by Jacobo
rawTumor_Cyt <- read_delim("data/Sauyeun_PDX/rawTumor_Cyt.tsv",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
) %>%
  dplyr::select(-`!`)

### raw count processed by Remy Nicolle
geneCount_raw_28s_totalRNA <- read_delim("data/Sauyeun_PDX/geneCount_raw_28s_totalRNA.tsv",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
) %>%
  rename_with(~ coalesce(sample_names_equivalences$CITID[match(., sample_names_equivalences$fastq_name)], .)) %>%
  dplyr::select(-`!`)

### best subset selected data direct input on PCA
ISRact_data_bestsubset <- read_delim("data/Sauyeun_PDX/5vs5_training_TMM_nonstandardized_bestsubset.tsv",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
) %>%
  dplyr::rename(sample = rowname)

### best overall selected data direct input on PCA
ISRact_data_bestoverall <- read_delim("data/Sauyeun_PDX/5vs5_training_TMM_nonstandardized_bestoverall.tsv",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
) %>%
  dplyr::rename(sample = rowname)


################################################################################
# PARAMETERS
Sauyeun_raw <- rawTumor_Cyt # geneCount_raw_28s_totalRNA | rawTumor_Cyt
cormethod <- "pearson"
keepgenes <- 1000
select_method <- "bestsubset" # bestsubset | bestoverall
norm_method <- "upperquartile" # TMM | upperquartile

# Variables for filter function
filter <- FALSE # exclude just the low purity if FALSE
estimation_method <- "histology" # histology | deconvolution
neoplastic_min <- 0.2
acinar_max <- 0.05
islet_max <- 1
tumor_tissue_min <- 0.8
################################################################################

# Translate EnsemblID to gene names
## Version 75 for PDX data
ensembl75 <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  version = 75
) # listAttributes(ensembl75, page="feature_page")

annot_ensembl75 <- getBM(attributes = c(
  "ensembl_gene_id",
  "external_gene_id"
), mart = ensembl75)

translate <- deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

# Add a Gene names column to Sauyeun_raw
Sauyeun_raw$EnsemblID %>%
  translate[.] %>%
  make.names(unique = TRUE) -> Sauyeun_raw$Genenames

# Processing and Normalization
## Normalize raw counts
Sauyeun_norm_ <- dplyr::select(Sauyeun_raw, -Genenames) %>%
  column_to_rownames("EnsemblID") %>%
  DGEList() %>%
  calcNormFactors(method = norm_method) %>%
  cpm(log = TRUE)

## Export version for comparison
# as_tibble(Sauyeun_norm_, rownames = "EnsemblID") %>% write_tsv("data/Sauyeun_PDX/PDX_PCA_input_data.tsv") #! remove?

## prepare for PCA
Sauyeun_norm <- Sauyeun_norm_ %>%
  t() %>%
  as_tibble(rownames = "sample")

### Visualize it
pivot_longer(Sauyeun_norm,
  cols = 2:length(Sauyeun_norm),
  names_to = "EnsemblID",
  values_to = "Expression",
  names_repair = "unique"
) %>%
  ggplot(data, mapping = aes(x = sample, y = Expression, fill = sample)) +
  geom_violin() +
  theme(legend.position = "none") +
  coord_flip()

## sample and gene filtering
### get 5 and 5 most extreme samples for ICA3
top_samples <- arrange(sample_info, ICA3) %>%
  dplyr::slice(unique(c(1:5, n() - 0:4))) %>%
  mutate(ISRact_bin = ifelse(ICA3 < 0, "low_ISRact", "high_ISRact")) %>%
  arrange(sample)

ggdensity(
  sample_info,
  x = "ICA3",
  rug = TRUE
) + geom_rug(
  data = top_samples,
  inherit.aes = F,
  aes(x = ICA3, color = ISRact_bin)
)

### get the most correlated genes with ICA3
top_genes <- mostCorrGenes(
  Sauyeun_norm,
  select_method,
  keepgenes,
  cormethod,
  top_samples,
  sample_info
)

### Final subset
Sauyeun_norm_subset <- dplyr::select(Sauyeun_norm, sample, top_genes$EnsemblID) %>%
  dplyr::filter(sample %in% top_samples$sample) %>%
  arrange(sample)


# PCA
pca_pdx <- Sauyeun_norm_subset %>%
  column_to_rownames("sample") %>%
  scale() %>%
  prcomp()


## Visualize
fviz_eig(pca_pdx, addlabels = TRUE, ylim = c(0, 70)) # visualize explained variance

fviz_pca_ind(pca_pdx,
  col.ind = top_samples$ISRact_bin, # color by groups
  palette = c("#00FF00", "#FF0000"),
  # addEllipses = TRUE, # Ellipses
  legend.title = "Groups",
  invisible = "quali"
)

## export PCA for projection
# write_rds(pca_pdx, "data/Classifiers/pca_pdx.RDS")
pca_pdx <- read_rds("data/Classifiers/pca_pdx_ENZO.RDS") #read originally produced PCA by Enzo

## Predict Full Shin et al cohort for further comparisons
Shin_PCAspace <- predict(pca_pdx, Sauyeun_norm) %>%
  as_tibble() %>%
  mutate(sample = Sauyeun_norm$sample, .before = 1) %>%
  left_join(top_samples[, c("sample", "ISRact_bin")]) %>%
  mutate(ISRact_bin = replace_na(ISRact_bin, "intermediate_ISRact"))

# write_rds(Shin_PCAspace, "results/pca_pdx_projection_middle.RDS") #! remove?

Shin_ISRactPCA <- arrange(Shin_PCAspace, PC1) %>%
  dplyr::filter(!sample == "!") %>%
  mutate(ISRactPCA_bin = cut(.$PC1, breaks = c(quantile(.$PC1, c(0:3 / 3))), labels = c("low_ISRactPCA", "medium_ISRactPCA", "high_ISRactPCA"), include.lowest = TRUE)) %>%
  dplyr::select(sample, PC1, ISRactPCA_bin, ISRact_bin) %>%
  inner_join(sample_info[, c("sample", "PAMG", "ICA3", "Diabetes")], by = "sample") %>%
  dplyr::rename(ISRact = "ICA3", ISRactPCA = "PC1")

## export object for characterization
write_tsv(Shin_ISRactPCA, "results/Sauyeun_PDX/Shin_ISRactPCA.tsv")

# Plot comparisons with Basal/Classical and ISRact
## ISR vs ISRactPCA
scatter_shin_isractvsisractpca <- Shin_ISRactPCA %>%
  # dplyr::rename(ISRactPCA = "PC1", ISRact = "ICA3") %>%
  correlation_plotter(data = ., col1 = "ISRact", col2 = "ISRactPCA", data_name = "Shin et al. PDX")

ggsave(scatter_shin_isractvsisractpca,
  filename = "results/Figures/scatter_shin_isractvsisractpca.svg",
  width = 2,
  height = 2
)

## PAMG vs ISRactPCA
scatter_shin_pamgvsisractpca <- correlation_plotter(data = Shin_ISRactPCA, col1 = "PAMG", col2 = "ISRactPCA", data_name = "Shin et al. PDX")

ggsave(scatter_shin_pamgvsisractpca,
  filename = "results/Figures/scatter_shin_pamgvsisractpca.svg",
  width = 2,
  height = 2
)

## ISRact vs PAMG
scatter_shin_isractvspamg <- correlation_plotter(data = Shin_ISRactPCA, col1 = "ISRact", col2 = "PAMG", data_name = "Shin et al. PDX")

ggsave(scatter_shin_isractvspamg,
  filename = "results/Figures/scatter_shin_isractvspamg.svg",
  width = 2,
  height = 2
)

#-------------------------------------------------------------------------------
# Plot ISRact distribution

dist_isract <- dplyr::arrange(Shin_ISRactPCA, desc(ISRact)) %>%
  ggplot(aes(x = ISRact)) +
  geom_density() +
  geom_rug(aes(color = ISRact_bin)) +
  scale_colour_manual(values = c(low_ISRact = "seagreen", high_ISRact = "tomato3", intermediate_ISRact = "grey")) +
  theme_classic()

ggsave(dist_isract,
  filename = "results/Figures/dist_isract.svg",
  width = 5,
  height = 2
)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Plot projecion of full cohort
pca_training <- pca_pdx[["x"]] %>% # ! rename as pca_training
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  inner_join(top_samples[, c("sample", "ISRact_bin")])

show_projection <- bind_rows(
  as_tibble(pca_training) %>% mutate(dataset = "Shin et al. PDX"),
  mutate(Shin_PCAspace,
    sample = paste0(sample, "_test"),
    dataset = "Shin et al. PDX"
  )
) %>%
  dplyr::rename(ISRactPCA = "PC1")

projection_scatter_shin_full <- dplyr::mutate(show_projection) %>%
  dplyr::filter(
    dataset %in% c("Shin et al. PDX", "PaCaOmics PDX", "CPTAC", "CCLE"), # "Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"
    ISRact_bin %in% c("low_ISRact", "high_ISRact", "intermediate_ISRact")
  ) %>% # 'low_ISRact', 'high_ISRact', 'intermediate_ISRact'
  ggplot(aes(x = ISRactPCA, y = PC2, color = ISRact_bin, shape = dataset)) +
  geom_point() +
  scale_shape_manual(values = c(`Shin et al. PDX` = 16, `PaCaOmics PDX` = 17, CPTAC = 15, CCLE = 18)) +
  scale_color_manual(values = c(low_ISRact = "seagreen", high_ISRact = "tomato3", intermediate_ISRact = "grey", Unknown = "#619CFF")) + # c('low_ISRact', 'high_ISRact', 'intermediate_ISRact', 'Unknown'))unique(projection_ccle$primary_tissue))
  # ylim(-250,15) +
  theme_bw()

ggsave(projection_scatter_shin_full,
  filename = "results/Figures/projection_scatter_shin_full.svg",
  width = 7,
  height = 3
)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Check ISR marker genes
marker_genes <- c("PHGDH", "CBS", "IL18")

Gene_ISRactPCA_comparison <- as_tibble(Sauyeun_norm_, rownames = "EnsemblID") %>%
  dplyr::mutate(Gene = translate[EnsemblID]) %>%
  dplyr::select(-EnsemblID) %>%
  pivot_longer(-Gene, names_to = "sample") %>%
  dplyr::filter(Gene %in% marker_genes) %>%
  inner_join(Shin_ISRactPCA, by = "sample") %>%
  dplyr::mutate(ISRact_bin = fct(ISRact_bin, levels = c("low_ISRact", "intermediate_ISRact", "high_ISRact"))) %>%
  dplyr::filter(ISRact_bin != "intermediate_ISRact")




for (mgene in marker_genes) {
  mgene_plot <- dplyr::filter(Gene_ISRactPCA_comparison, Gene == mgene) %>%
    ggplot(aes(y = value, x = ISRact_bin, fill = ISRact_bin)) +
    geom_violin() +
    scale_fill_manual(values = c(low_ISRact = "seagreen", high_ISRact = "tomato3")) +
    rotate_x_text(90) +
    geom_boxplot(width = 0.1, fill = "white") +
    labs(title = mgene)

  ggsave(mgene_plot,
    filename = paste0("results/Figures/shin_", mgene, ".svg"),
    width = 3,
    height = 2
  )
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Gene weight comparison with PAMG
PAMG_gw <- as_tibble(pdacmolgrad:::.molGradSys$PDX$gw[pdacmolgrad:::.molGradSys$PDX$k], rownames = "GeneName") %>%
  dplyr::rename(PAMG = "ICA1")

ISRactPCA_gw <- as_tibble(pca_pdx$rotation, rownames = "EnsemblID") %>%
  mutate(GeneName = translate[EnsemblID]) %>%
  dplyr::select(GeneName, PC1) %>%
  dplyr::rename(ISRactPCA = PC1)

## discrete
compare <- list(PAMG_gw$GeneName, ISRactPCA_gw$GeneName)

compare_names <- c("PAMG", "ISRactPCA")
my_breaks <- c(500, 1000, 5000, 10000, 20000)
isractpcavspamg_venn <- ggVennDiagram(
  compare,
  category.names = compare_names,
  label = "count",
  label_geom = "text"
) + scale_fill_gradient(low = "grey90", high = "red", trans = "log", breaks = my_breaks)

ggsave(isractpcavspamg_venn,
  filename = paste0("results/Figures/isractpcavspamg_venn.svg"),
  width = 4,
  height = 2
)
## Quantitative comparison of common genes
PAMG_ISRact_quantc <- full_join(PAMG_gw, ISRactPCA_gw, by = "GeneName") %>%
  mutate(na_PAMG = is.na(PAMG), na_ISRactPCA = is.na(ISRactPCA))

PAMG_ISRact_quantc %>%
  drop_na() %>%
  correlation_plotter(col1 = "ISRactPCA", col2 = "PAMG", data_name = "gene weights")

## Whats the gene weight of ISRactPCA and PAMG unique genes?
## ISRact vs PAMG presence
### positive
PAMG_ISRact_quantc_ISRactpos <- dplyr::filter(PAMG_ISRact_quantc, ISRactPCA > 0 | is.na(ISRactPCA))

PAMG_absence_test <- wilcox.test(
  x = dplyr::filter(PAMG_ISRact_quantc_ISRactpos, na_PAMG == T) %>% dplyr::select(ISRactPCA) %>% deframe() %>% abs(),
  y = dplyr::filter(PAMG_ISRact_quantc_ISRactpos, na_PAMG == F) %>% dplyr::select(ISRactPCA) %>% deframe() %>% abs(),
  alternative = "two.sided"
)

isractpos_commonvsunique <- ggplot(PAMG_ISRact_quantc_ISRactpos, aes(x = na_PAMG, y = ISRactPCA)) +
  geom_violin(aes(fill = na_PAMG)) +
  scale_fill_grey(start = 0.4, end = 0.9) +
  geom_boxplot(width = 0.1) +
  scale_x_discrete(labels = c("common", "unique ISRactPCA")) +
  labs(
    title = "Are ISRactPCA positive contributions different between common and unique genes",
    subtitle = paste0(
      "Wilcoxon test two sided p-value = ",
      formatC(PAMG_absence_test$p.value, format = "e", digits = 2)
    )
  ) +
  theme_bw()

ggsave(isractpos_commonvsunique,
  filename = paste0("results/Figures/isractpos_commonvsunique.svg"),
  width = 4,
  height = 2
)

### negative
PAMG_ISRact_quantc_ISRactneg <- dplyr::filter(PAMG_ISRact_quantc, ISRactPCA < 0 | is.na(ISRactPCA))

PAMG_absence_test <- wilcox.test(
  x = dplyr::filter(PAMG_ISRact_quantc_ISRactneg, na_PAMG == T) %>% dplyr::select(ISRactPCA) %>% deframe() %>% abs(),
  y = dplyr::filter(PAMG_ISRact_quantc_ISRactneg, na_PAMG == F) %>% dplyr::select(ISRactPCA) %>% deframe() %>% abs(),
  alternative = "two.sided"
)

isractneg_commonvsunique <- ggplot(PAMG_ISRact_quantc_ISRactneg, aes(x = na_PAMG, y = ISRactPCA)) +
  geom_violin(aes(fill = na_PAMG)) +
  scale_fill_grey(start = 0.4, end = 0.9) +
  geom_boxplot(width = 0.1) +
  scale_x_discrete(labels = c("common", "unique for ISRactPCA")) +
  labs(
    title = "Are ISRactPCA negative contributions different between common and unique genes",
    subtitle = paste0(
      "Wilcoxon test two sided p-value = ",
      formatC(PAMG_absence_test$p.value, format = "e", digits = 2)
    )
  ) +
  theme_bw()

ggsave(isractneg_commonvsunique,
  filename = paste0("results/Figures/isractneg_commonvsunique.svg"),
  width = 4,
  height = 2
)


## PAMG vs ISRact presence
### positive
PAMG_ISRact_quantc_PAMGpos <- dplyr::filter(PAMG_ISRact_quantc, PAMG > 0 | is.na(PAMG))
ISRactPCA_absence_test <- wilcox.test(
  x = dplyr::filter(PAMG_ISRact_quantc_PAMGpos, na_ISRactPCA == T) %>% dplyr::select(PAMG) %>% deframe() %>% abs(),
  y = dplyr::filter(PAMG_ISRact_quantc_PAMGpos, na_ISRactPCA == F) %>% dplyr::select(PAMG) %>% deframe() %>% abs(),
  alternative = "two.sided"
)
pamgpos_commonvsunique <- ggplot(PAMG_ISRact_quantc_PAMGpos, aes(x = na_ISRactPCA, y = PAMG)) +
  geom_violin(aes(fill = na_ISRactPCA)) +
  scale_fill_grey(start = 0.4, end = 0.9) +
  geom_boxplot(width = 0.1) +
  scale_x_discrete(labels = c("common", "unique for PAMG")) +
  labs(
    title = "Are PAMG positive contributions different between common and unique genes",
    subtitle = paste0(
      "Wilcoxon test two sided p-value = ",
      formatC(ISRactPCA_absence_test$p.value, format = "e", digits = 2)
    )
  ) +
  theme_bw()

ggsave(pamgpos_commonvsunique,
  filename = paste0("results/Figures/pamgpos_commonvsunique.svg"),
  width = 4,
  height = 2
)

### negative
PAMG_ISRact_quantc_PAMGneg <- dplyr::filter(PAMG_ISRact_quantc, PAMG < 0 | is.na(PAMG))
ISRactPCA_absence_test <- wilcox.test(
  x = dplyr::filter(PAMG_ISRact_quantc_PAMGneg, na_ISRactPCA == T) %>% dplyr::select(PAMG) %>% deframe() %>% abs(),
  y = dplyr::filter(PAMG_ISRact_quantc_PAMGneg, na_ISRactPCA == F) %>% dplyr::select(PAMG) %>% deframe() %>% abs(),
  alternative = "two.sided"
)
pamgneg_commonvsunique <- ggplot(PAMG_ISRact_quantc_PAMGneg, aes(x = na_ISRactPCA, y = PAMG)) +
  geom_violin(aes(fill = na_ISRactPCA)) +
  scale_fill_grey(start = 0.4, end = 0.9) +
  geom_boxplot(width = 0.1) +
  scale_x_discrete(labels = c("common", "unique for PAMG")) +
  labs(
    title = "Are PAMG negative contributions different between common and unique genes",
    subtitle = paste0(
      "Wilcoxon test two sided p-value = ",
      formatC(ISRactPCA_absence_test$p.value, format = "e", digits = 2)
    )
  ) +
  theme_bw()

ggsave(pamgneg_commonvsunique,
  filename = paste0("results/Figures/pamgneg_commonvsunique.svg"),
  width = 4,
  height = 2
)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Gene weight comparison with GemPred
GemPred_gw <- read_xlsx("data/Other/NicolleAnnalsOncology2021_TableS3.xlsx", skip = 2) %>%
  dplyr::rename(GeneName = "Gene", GemPred = "Weight")

signature_var <- "ISRactICA" # "ISRactPCA"
if (signature_var == "ISRactPCA") {
  signature_gw <- as_tibble(pca_pdx$rotation, rownames = "EnsemblID") %>%
    mutate(GeneName = translate[EnsemblID]) %>%
    dplyr::select(GeneName, PC1) %>%
    dplyr::rename(ISRactPCA = PC1)
} else if (signature_var == "ISRactICA") {
  signature_gw <- read_rds("../06_Human_Cohort_Projection/01_PDXTranslation_to_PDXTrascription/ISRactICA_IC3.RDS")$S %>%
    as_tibble(., rownames = "EnsemblID") %>%
    mutate(GeneName = translate[EnsemblID]) %>%
    dplyr::select(GeneName, IC.3) %>%
    dplyr::rename(ISRactICA = "IC.3")
}

## discrete
compare <- list(GemPred_gw$GeneName, signature_gw$GeneName)

compare_names <- c("GemPred", signature_var)
my_breaks <- c(500, 1000, 5000, 10000, 20000)
isractpcavsGemPred_venn <- ggVennDiagram(
  compare,
  category.names = compare_names,
  label = "count",
  label_geom = "text"
) + scale_fill_gradient(low = "grey90", high = "red", trans = "log", breaks = my_breaks)

# ggsave(isractpcavsGemPred_venn,
#        filename = paste0("results/Figures/isractpcavsGemPred_venn.svg"),
#        width = 4,
#        height = 2
# )
## Quantitative comparison of common genes
GemPred_ISRact_quantc <- full_join(GemPred_gw, signature_gw, by = "GeneName") %>%
  mutate(na_GemPred = is.na(GemPred), na_ISRact_signature = is.na(get(signature_var)))

GemPred_ISRact_quantc %>%
  drop_na() %>%
  correlation_plotter(col1 = signature_var, col2 = "GemPred", data_name = "gene weights")

## Whats the gene weight of ISRactPCA and GemPred unique genes?
## ISRact vs GemPred presence
### positive
GemPred_ISRact_quantc_ISRactpos <- dplyr::filter(GemPred_ISRact_quantc, get(signature_var) > 0 | is.na(get(signature_var)))

GemPred_absence_test <- wilcox.test(
  x = dplyr::filter(GemPred_ISRact_quantc_ISRactpos, na_GemPred == T) %>% dplyr::select(all_of(signature_var)) %>% deframe() %>% abs(),
  y = dplyr::filter(GemPred_ISRact_quantc_ISRactpos, na_GemPred == F) %>% dplyr::select(all_of(signature_var)) %>% deframe() %>% abs(),
  alternative = "two.sided"
)

isractpos_commonvsunique <- ggplot(GemPred_ISRact_quantc_ISRactpos, aes(x = na_GemPred, y = get(signature_var))) +
  geom_violin(aes(fill = na_GemPred)) +
  scale_fill_grey(start = 0.4, end = 0.9) +
  geom_boxplot(width = 0.1) +
  scale_x_discrete(labels = c("common", paste0("unique ", signature_var))) +
  labs(
    title = paste0("Are ", signature_var, " positive contributions different between common and unique genes"),
    subtitle = paste0(
      "Wilcoxon test two sided p-value = ",
      formatC(GemPred_absence_test$p.value, format = "e", digits = 2)
    )
  ) +
  theme_bw() +
  ylab(signature_var)

# ggsave(isractpos_commonvsunique,
#        filename = paste0("results/Figures/isractpos_commonvsunique.svg"),
#        width = 4,
#        height = 2
# )

### negative
GemPred_ISRact_quantc_ISRactneg <- dplyr::filter(GemPred_ISRact_quantc, get(signature_var) < 0 | is.na(get(signature_var)))

GemPred_absence_test <- wilcox.test(
  x = dplyr::filter(GemPred_ISRact_quantc_ISRactneg, na_GemPred == T) %>% dplyr::select(all_of(signature_var)) %>% deframe() %>% abs(),
  y = dplyr::filter(GemPred_ISRact_quantc_ISRactneg, na_GemPred == F) %>% dplyr::select(all_of(signature_var)) %>% deframe() %>% abs(),
  alternative = "two.sided"
)

isractneg_commonvsunique <- ggplot(GemPred_ISRact_quantc_ISRactneg, aes(x = na_GemPred, y = get(signature_var))) +
  geom_violin(aes(fill = na_GemPred)) +
  scale_fill_grey(start = 0.4, end = 0.9) +
  geom_boxplot(width = 0.1) +
  scale_x_discrete(labels = c("common", paste0("unique ", signature_var))) +
  labs(
    title = paste0("Are ", signature_var, " negative contributions different between common and unique genes"),
    subtitle = paste0(
      "Wilcoxon test two sided p-value = ",
      formatC(GemPred_absence_test$p.value, format = "e", digits = 2)
    )
  ) +
  theme_bw() +
  ylab(signature_var)

# ggsave(isractneg_commonvsunique,
#        filename = paste0("results/Figures/isractneg_commonvsunique.svg"),
#        width = 4,
#        height = 2
# )


## GemPred vs ISRact presence
### positive
GemPred_ISRact_quantc_GemPredpos <- dplyr::filter(GemPred_ISRact_quantc, GemPred > 0 | is.na(GemPred))
ISRact_absence_test <- wilcox.test(
  x = dplyr::filter(GemPred_ISRact_quantc_GemPredpos, na_ISRact_signature == T) %>% dplyr::select(GemPred) %>% deframe() %>% abs(),
  y = dplyr::filter(GemPred_ISRact_quantc_GemPredpos, na_ISRact_signature == F) %>% dplyr::select(GemPred) %>% deframe() %>% abs(),
  alternative = "two.sided"
)
GemPredpos_commonvsunique <- ggplot(GemPred_ISRact_quantc_GemPredpos, aes(x = na_ISRact_signature, y = GemPred)) +
  geom_violin(aes(fill = na_ISRact_signature)) +
  scale_fill_grey(start = 0.4, end = 0.9) +
  geom_boxplot(width = 0.1) +
  scale_x_discrete(labels = c("common", "unique for GemPred")) +
  labs(
    title = "Are GemPred positive contributions different between common and unique genes",
    subtitle = paste0(
      "Wilcoxon test two sided p-value = ",
      formatC(ISRact_absence_test$p.value, format = "e", digits = 2)
    )
  ) +
  theme_bw()

# ggsave(GemPredpos_commonvsunique,
#        filename = paste0("results/Figures/GemPredpos_commonvsunique.svg"),
#        width = 4,
#        height = 2
# )

### negative
GemPred_ISRact_quantc_GemPredneg <- dplyr::filter(GemPred_ISRact_quantc, GemPred < 0 | is.na(GemPred))
ISRact_absence_test <- wilcox.test(
  x = dplyr::filter(GemPred_ISRact_quantc_GemPredneg, na_ISRact_signature == T) %>% dplyr::select(GemPred) %>% deframe() %>% abs(),
  y = dplyr::filter(GemPred_ISRact_quantc_GemPredneg, na_ISRact_signature == F) %>% dplyr::select(GemPred) %>% deframe() %>% abs(),
  alternative = "two.sided"
)
GemPredneg_commonvsunique <- ggplot(GemPred_ISRact_quantc_GemPredneg, aes(x = na_ISRact_signature, y = GemPred)) +
  geom_violin(aes(fill = na_ISRact_signature)) +
  scale_fill_grey(start = 0.4, end = 0.9) +
  geom_boxplot(width = 0.1) +
  scale_x_discrete(labels = c("common", "unique for GemPred")) +
  labs(
    title = "Are GemPred negative contributions different between common and unique genes",
    subtitle = paste0(
      "Wilcoxon test two sided p-value = ",
      formatC(ISRact_absence_test$p.value, format = "e", digits = 2)
    )
  ) +
  theme_bw()

# ggsave(GemPredneg_commonvsunique,
#        filename = paste0("results/Figures/GemPredneg_commonvsunique.svg"),
#        width = 4,
#        height = 2
# )
#-------------------------------------------------------------------------------
