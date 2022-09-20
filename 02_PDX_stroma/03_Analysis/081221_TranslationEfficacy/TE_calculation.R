library(tidyverse)
library(Hmisc)
library(corrplot)
library(lmtest)
library(ggpubr)
library(Biobase)
library(factoextra)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/081221_TranslationEfficacy/")
source("functions.R")

# Data loading
normHost_cyt <- read_tsv("../../00_Data/Processed_data/normHost_Cyt.tsv")
normHost_pol <- read_tsv("../../00_Data/Processed_data/normHost_Pol.tsv")
sample_info <- read_tsv("../../00_Data/Processed_data/sample_info.tsv") %>% column_to_rownames("sample")

# Joint dataframe creation
sn_cyt <- paste(colnames(normHost_cyt[-1]), "cyt", sep = "_")
sn_pol <- paste(colnames(normHost_pol[-1]), "pol", sep = "_")

normHost_cyt.tmp <- normHost_cyt
normHost_pol.tmp <- normHost_pol

colnames(normHost_cyt.tmp) <- c("EnsemblID", sn_cyt)
colnames(normHost_pol.tmp) <- c("EnsemblID", sn_pol)

norm_host <- inner_join(normHost_cyt.tmp, normHost_pol.tmp, by = "EnsemblID")

# Phenotype table creation
sn <- c(colnames(normHost_cyt[-1]), colnames(normHost_pol[-1]))

frac <- c(
  rep("cyt", length(sn_cyt)),
  rep("pol", length(sn_pol))
)

rna_conc <- c(
  sample_info[colnames(normHost_cyt[-1]), "RNAconc_cyt"],
  sample_info[colnames(normHost_pol[-1]), "RNAconc_pol"]
)

diab_stat <- c(
  sample_info[colnames(normHost_cyt[-1]), "Diabetes"],
  sample_info[colnames(normHost_pol[-1]), "Diabetes"]
)

tumor_ica <- bind_rows(
  sample_info[colnames(normHost_cyt[-1]), str_subset(colnames(sample_info), "ICA")],
  sample_info[colnames(normHost_pol[-1]), str_subset(colnames(sample_info), "ICA")],
)

rownames(tumor_ica) <- c(sn_cyt, sn_pol)

phenotype.tmp <- data.frame(
  sample = factor(sn),
  fraction = factor(frac),
  rna_conc = as.double(rna_conc),
  diabetes = factor(diab_stat)
)

rownames(phenotype.tmp) <- colnames(norm_host)[-1]

phenotype <- merge(phenotype.tmp, tumor_ica, by = "row.names", ) %>% column_to_rownames("Row.names")
phenotype <- phenotype[colnames(norm_host)[-1], ]


# PCA analysis
norm_host %>%
  column_to_rownames("EnsemblID") %>%
  data.matrix() %>%
  t() %>%
  prcomp(scale = T) -> norm_pca

var_explained <- norm_pca$sdev^2 / sum(norm_pca$sdev^2)
fviz_eig(norm_pca,
  barfill = "lightgrey",
  barcolor = "black", title = "PDX cyt/pol PCA",
  subtitle = paste(
    "90% variance reached at the",
    which(cumsum(var_explained) > 0.9)[1],
    "th component"
  )
)


all(rownames(norm_pca$x) == rownames(phenotype)) %>% stopifnot()
pca_toplot <- cbind(norm_pca$x, phenotype, by = "row.names")

plot_PCs(pca_toplot, "fraction", 10, 0.5)
plot_PCs(pca_toplot, "diabetes", 10, 0.5)
plot_PCs(pca_toplot, "sample", 10, 0.5)
plot_PCs(pca_toplot, "rna_conc", 10, 0.5)
plot_PCs(pca_toplot, "ICA1", 10, 0.5)
plot_PCs(pca_toplot, "ICA2", 10, 0.5)
plot_PCs(pca_toplot, "ICA3", 10, 0.5)
plot_PCs(pca_toplot, "ICA4", 10, 0.5)
plot_PCs(pca_toplot, "ICA5", 10, 0.5)
plot_PCs(pca_toplot, "ICA6", 10, 0.5)
while (dev.cur() > 1) dev.off()

## is PDAC001T an outlier?
# pca_toplot["PDAC01Tout"] <- "other"
# pca_toplot[pca_toplot["sample"] == "PDAC001T", "PDAC01Tout"] <- "PDAC001T"
# plot_PCs(pca_toplot, "PDAC01Tout", 10, 0.5)

# Similarity through whole genome correlation
statistic <- "pearson"
norm_host %>%
  column_to_rownames("EnsemblID") %>%
  data.matrix() %>%
  rcorr(., type = statistic) -> corr
corrplot(
  corr = corr$r,
  p.mat = corr$P,
  is.corr = F, order = "hclust",
  type = "lower",
  main = paste(statistic, sep = ": ")
)

# TE calculation
all(normHost_cyt$EnsemblID == normHost_pol$EnsemblID) %>% stopifnot()
all(colnames(normHost_cyt) == colnames(normHost_pol)) %>% stopifnot()
lapply(normHost_cyt$EnsemblID, calculaTE) %>%
  as.data.frame(row.names = c(
    colnames(normHost_pol),
    "slope", "Phomo", "Pnorm", "Outlier_P", "Outlier_s"
  )) %>%
  t() %>%
  as_tibble() %>%
  column_to_rownames("EnsemblID") %>%
  mutate_at(vars(-Outlier_s), function(x) as.numeric(as.character(x))) -> TEs.unf

## Pre-filtering tests
### Slope distribution
TEs.unf$slope %>%
  between(-1, 1) %>%
  sum() / nrow(TEs.unf)
TEs.unf %>%
  dplyr::select(slope) %>%
  ggplot(aes(x = slope)) +
  geom_density() +
  ggtitle("lm slope distribution before filtering") +
  theme_classic()

### Outliers
TEs.unf %>%
  dplyr::select(!c(slope, Phomo, Pnorm, Outlier_P, Outlier_s)) %>%
  data.matrix() %>%
  OutlierTest(residualMatrix = .)

TEs.unf %>%
  dplyr::select(c(Outlier_P, Outlier_s)) %>%
  ggplot(aes(x = Outlier_P)) +
  geom_density() +
  scale_x_continuous(trans = "log2") +
  theme_classic()

TEs.unf %>%
  dplyr::select(c(Outlier_P, Outlier_s)) %>%
  dplyr::filter(Outlier_P > 0.05) %>%
  ggplot(aes(x = Outlier_s)) +
  geom_bar() +
  ggtitle("most extreme sample, non filtered genes") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

TEs.unf %>%
  dplyr::select(c(Outlier_P, Outlier_s)) %>%
  dplyr::filter(Outlier_P < 0.05) %>%
  ggplot(aes(x = Outlier_s)) +
  geom_bar() +
  ggtitle("most extreme sample, outlier genes") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## TE filtering
TEs.unf %>%
  dplyr::filter(Phomo > 0.05 & Pnorm > 0.05 & Outlier_P > 0.05) %>%
  dplyr::select(!c(slope, Phomo, Pnorm, Outlier_P, Outlier_s)) -> TEs

## Post-filtering tests
### Slope distribution
TEs.unf[rownames(TEs), ]$slope %>%
  between(-1, 1) %>%
  sum() / nrow(TEs)
TEs.unf[rownames(TEs), ] %>%
  dplyr::select(slope) %>%
  ggplot(aes(x = slope)) +
  geom_density() +
  ggtitle("lm slope distribution after filtering") +
  theme_classic()


### Outliers
TEs %>%
  data.matrix() %>%
  OutlierTest(residualMatrix = .)


### filtering relation with gene expression bias
#### one sample
sample <- "PDAC021T"

stopifnot(normHost_cyt$EnsemblID == normHost_pol$EnsemblID)
bind_cols(normHost_cyt[c("EnsemblID", sample)], normHost_pol[sample]) %>%
  rename(cyt = 2, pol = 3) %>%
  mutate(filtered = if_else(EnsemblID %in% rownames(TEs), "No", "Yes")) %>%
  ggplot(aes(x = cyt, y = pol, color = filtered)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_rug(alpha = 0.1, size = 1.5) +
  theme_bw() +
  theme(legend.position = "top") +
  labs(title = sample)

#### all samples (density)
metric_gene_distributions <- function(metric, title) {
  ### Cytosolic
  normHost_cyt %>%
    column_to_rownames("EnsemblID") %>%
    apply(1, metric) %>%
    as.list() %>%
    as_tibble() %>%
    pivot_longer(cols = starts_with("ENS")) %>%
    mutate(filtered = if_else(name %in% rownames(TEs), "No", "Yes")) %>%
    ggplot(aes(x = value, fill = filtered)) +
    geom_density(alpha = 0.4) +
    ggtitle("Cytosolic") +
    theme_classic() -> bias_cyt

  ### Polysome
  normHost_pol %>%
    column_to_rownames("EnsemblID") %>%
    apply(1, metric) %>%
    as.list() %>%
    as_tibble() %>%
    pivot_longer(cols = starts_with("ENS")) %>%
    mutate(filtered = if_else(name %in% rownames(TEs), "No", "Yes")) %>%
    ggplot(aes(x = value, fill = filtered)) +
    geom_density(alpha = 0.4) +
    ggtitle("Polysome") +
    theme_classic() -> bias_pol

  ### TEs
  TEs.unf %>%
    dplyr::select(!c(slope, Phomo, Pnorm, Outlier_P, Outlier_s)) %>%
    apply(1, metric) %>%
    as.list() %>%
    as_tibble() %>%
    pivot_longer(cols = starts_with("ENS")) %>%
    mutate(filtered = if_else(name %in% rownames(TEs), "No", "Yes")) %>%
    ggplot(aes(x = value, fill = filtered)) +
    geom_density(alpha = 0.4) +
    ggtitle("TEs") +
    theme_classic() -> bias_TE

  metric_densplots <- ggarrange(bias_cyt, bias_pol, bias_TE + rremove("x.text"),
    labels = c("A", "B", "C"),
    ncol = 2, nrow = 2
  )

  annotate_figure(metric_densplots,
    top = text_grob(title, color = "black", face = "bold", size = 14)
  )
}
metric_gene_distributions(min, "Gene min counts")
metric_gene_distributions(max, "Gene max counts")
metric_gene_distributions(median, "Gene median counts")
metric_gene_distributions(IQR, "Variability of expression (genes IQR distribution)")
metric_gene_distributions(sum, "Gene sum of counts")


# Export Filtered Translation Efficacies
TEs %>%
  rownames_to_column("EnsemblID") %>%
  write_tsv("02_Output/TEs.tsv")
