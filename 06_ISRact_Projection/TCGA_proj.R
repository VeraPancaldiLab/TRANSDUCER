# ISRactPCA projection into TCGA
################################################################################
#' Import different files
#' Filter data and project ISRactPCA
#' Produce projection scatterplot
#' Explore relation with PAMG sample score
#' Survival analysis for ISRactPCA
#' Survival analysis for PHGDH/CBS
#' Explore relation with BRCAness and mutation status
################################################################################

# Import library
library(tidyverse)
library(readxl)
library(biomaRt)
library(pdacmolgrad) # devtools::install_github("RemyNicolle/pdacmolgrad")
library(edgeR)
library(Hmisc)
library(ggrepel)
library(DEqMS)
library(matrixStats)
library(msigdbr)
library(fgsea)
library(ggpubr)
library(data.table)
library(ggplot2)
library(survival)
library(survminer)
library(readxl)
library(pheatmap)

################################################################################
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")

# Import datasets
## TCGA
### Data
TCGA <- read_rds("data/TCGA/TCGA_PAAD.rds") %>% as_tibble(rownames = "Gene")

### Metadata
TCGA_clinical <- read_rds("data/TCGA/TCGA_PAAD_survival.rds")

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
# (no variable parameters)
################################################################################

# Translate EnsemblID to gene names
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = "GRCh37") # listAttributes(ensembl75, page="feature_page")
annot_ensembl <- getBM(attributes = c(
  "ensembl_gene_id",
  "external_gene_name"
), mart = ensembl)

gene_to_ensembl <- deframe(annot_ensembl[c("external_gene_name", "ensembl_gene_id")])
ensembl_to_gene <- setNames(names(gene_to_ensembl), gene_to_ensembl)

# Preprocessing
## translate gene names
TCGA_ensembl <- mutate(TCGA, EnsemblID = gene_to_ensembl[Gene]) %>%
  dplyr::select(-c("Gene")) %>%
  group_by(EnsemblID) %>%
  summarise_all(mean)

TCGA_gene <- TCGA %>%
  group_by(Gene) %>%
  summarise_all(mean)

## Create sample_info_TCGA
type_pamg <- projectMolGrad(newexp = column_to_rownames(TCGA_gene, "Gene"), geneSymbols = TCGA_gene$Gene) %>%
  as_tibble(rownames = "sample")

sample_info_TCGA <- dplyr::select(type_pamg, sample, ICGCrnaseq) %>%
  dplyr::rename(PAMG = ICGCrnaseq) %>%
  left_join(rownames_to_column(TCGA_clinical, "sample"), "sample")

# Projecting TCGA on ISRactPCA
#### transpose normalised data
TCGA_toproj <- TCGA_ensembl %>%
  pivot_longer(cols = -EnsemblID, names_to = "case_id", values_to = "Expression") %>%
  dplyr::filter(!str_detect(EnsemblID, "NA")) %>% # remove eventual NA for gene Ensembl
  pivot_wider(names_from = EnsemblID, values_from = Expression, names_repair = "minimal") %>%
  column_to_rownames("case_id")

### remove missing genes from the PCA before projection
filter_pca <- function(.data, objective) {
  as_tibble(.data, rownames = "tmp") %>%
    dplyr::filter(tmp %in% objective) %>%
    column_to_rownames("tmp") %>%
    data.matrix()
}
pca_pdx$rotation <- filter_pca(pca_pdx$rotation, names(TCGA_toproj)[-1])
pca_pdx$center <- filter_pca(pca_pdx$center, names(TCGA_toproj)[-1])
pca_pdx$scale <- filter_pca(pca_pdx$scale, names(TCGA_toproj)[-1])

TCGA_PCAspace <- predict(pca_pdx, TCGA_toproj) %>%
  as_tibble(rownames = "sample") %>%
  left_join(sample_info_TCGA, by = "sample") %>%
  dplyr::rename(ISRactPCA = "PC1") %>%
  mutate(
    ISRactPCA_bin = if_else(ISRactPCA < quantile(ISRactPCA, probs = 0.3333), "low_ISRactPCA",
      if_else(ISRactPCA < quantile(ISRactPCA, probs = 0.6666), "intermediate_ISRactPCA", "high_ISRactPCA")
    ),
    ISRactPCA_bin = fct(ISRactPCA_bin, levels = c("low_ISRactPCA", "intermediate_ISRactPCA", "high_ISRactPCA"))
  )


TCGA_ISRact_projection <- arrange(TCGA_PCAspace, ISRactPCA) %>%
  dplyr::select(sample, ISRactPCA, ISRactPCA_bin, !matches("PC"))

# Plot pca and projections
## Add ISR status to PCA df
pca_training <- pca_pdx[["x"]] %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  dplyr::rename(ISRactPCA = "PC1") %>%
  inner_join(top_samples[, c("sample", "ISRact_bin")])

#-------------------------------------------------------------------------------
# Plot projecion PacaOmics and then projection Proteogenomics

show_projection <- bind_rows(
  as_tibble(pca_training) %>% mutate(dataset = "Sauyeun PDX"),
  as_tibble(TCGA_PCAspace) %>% mutate(
    ISRact_bin = "Unknown",
    dataset = "TCGA"
  )
)

dplyr::mutate(show_projection, ISRact_bin = str_replace(ISRact_bin, "ICA3", "ISRact")) %>%
  dplyr::filter(
    dataset %in% c("Sauyeun PDX", "PACAOMICS PDX", "TCGA"), # "Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"
    ISRact_bin %in% c("low_ISRact", "high_ISRact", "intermediate_ISRact", "Unknown")
  ) %>%
  ggplot(aes(x = ISRactPCA, y = PC2, color = ISRact_bin, shape = dataset)) +
  geom_point() +
  scale_shape_discrete(limits = c("Sauyeun PDX", "PACAOMICS PDX", "TCGA")) +
  scale_color_discrete(limits = c("low_ISRact", "high_ISRact", "intermediate_ISRact", "Unknown")) + # c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown'))unique(projection_ccle$primary_tissue))
  geom_rect(
    xmin = min(TCGA_PCAspace$ISRactPCA),
    xmax = quantile(TCGA_PCAspace$ISRactPCA, probs = 0.3333),
    ymin = min(TCGA_PCAspace$PC2), ymax = max(TCGA_PCAspace$PC2), linewidth = 0, fill = "red", alpha = 0.002
  ) +
  geom_rect(
    xmin = quantile(TCGA_PCAspace$ISRactPCA, probs = 0.6666),
    xmax = max(TCGA_PCAspace$ISRactPCA),
    ymin = min(TCGA_PCAspace$PC2), ymax = max(TCGA_PCAspace$PC2), linewidth = 0, fill = "green", alpha = 0.002
  )
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Plot comparisons with Basal/Classical and ISRact
## PAMG vs ISRactPCA
correlation_plotter(data = TCGA_ISRact_projection, col1 = "PAMG", col2 = "ISRactPCA", data_name = "TCGA")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Survival curves regarding ISR status
surv_data <- dplyr::rename(TCGA_ISRact_projection, case_id = sample)
## OS
### Kaplan Meyer 33up vs 33down
fit <- survfit(Surv(OS, OS_event) ~ ISRactPCA_bin,
  data = dplyr::filter(surv_data, ISRactPCA_bin != "intermediate_ISRactPCA")
)
print(fit)

#### Change color, linetype by strata, risk.table color by strata
kmos_tcga_ISRact <- ggsurvplot(fit,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw(), # Change ggplot2 theme
  palette = c("tomato3", "seagreen")
)

print(kmos_tcga_ISRact)
ggsave(paste0("results/survplots/FIGURES/kmos_tcga_ISRact.svg"),
  kmos_tcga_ISRact$plot,
  height = 4, width = 4
)
### Cox Proportional hazzards model
cox.mod <- coxph(Surv(OS, OS_event) ~ ISRactPCA,
  data = as.data.frame(surv_data)
)
#### assumption checking
##### Linearity
plot(predict(cox.mod),
  residuals(cox.mod, type = "martingale"),
  xlab = "fitted", ylab = "Martingale residuals",
  main = "Residuals Plot", las = 1
)
abline(h = 0)
lines(
  smooth.spline(
    predict(cox.mod),
    residuals(cox.mod, type = "martingale")
  ),
  col = "red"
)

##### Proportional Hazzards
par(mfrow = c(1, 1))
plot(cox.zph(cox.mod)) # if failed, Proportional
abline(h = 0, col = 2)

#### Model plotting
ggforest(cox.mod)
ggadjustedcurves(cox.mod, data = as.data.frame(surv_data))


## PFS
### Kaplan Meyer 33up vs 33down
fit <- survfit(Surv(PFS, PFS_event) ~ ISRactPCA_bin,
  data = dplyr::filter(surv_data, ISRactPCA_bin != "intermediate_ISRactPCA")
)
print(fit)

#### Change color, linetype by strata, risk.table color by strata
kmpfs_tcga_ISRact <- ggsurvplot(fit,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw(), # Change ggplot2 theme
  palette = c("tomato3", "seagreen")
)

print(kmpfs_tcga_ISRact)
ggsave(paste0("results/survplots/FIGURES/kmpfs_tcga_ISRact.svg"),
  kmpfs_tcga_ISRact$plot,
  height = 4, width = 4
)

### Cox Proportional hazzards model
cox.mod <- coxph(Surv(PFS, PFS_event) ~ ISRactPCA,
  data = as.data.frame(surv_data)
)
#### assumption checking
##### Linearity
plot(predict(cox.mod),
  residuals(cox.mod, type = "martingale"),
  xlab = "fitted", ylab = "Martingale residuals",
  main = "Residuals Plot", las = 1
)
abline(h = 0)
lines(
  smooth.spline(
    predict(cox.mod),
    residuals(cox.mod, type = "martingale")
  ),
  col = "red"
)

##### Proportional Hazzards
par(mfrow = c(1, 1))
plot(cox.zph(cox.mod)) # if failed, Proportional
abline(h = 0, col = 2)

#### Model plotting
ggforest(cox.mod)
ggadjustedcurves(cox.mod, data = as.data.frame(surv_data))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Survival curves regarding PHGDH/CBS status
surv_data <- as_tibble(TCGA_toproj, rownames = "sample") %>%
  dplyr::select(sample, c(PHGDH = "ENSG00000092621", CBS = "ENSG00000160200")) %>%
  inner_join(TCGA_ISRact_projection, by = "sample")

measure <- "mean" # overlap | mean
## Calculate the measure to split samples
if (measure == "mean") {
  ### min mean
  surv_data <- dplyr::mutate(surv_data,
    meanPHGDHCBS = rowMeans(surv_data[c("PHGDH", "CBS")]),
    stratPHGDHCBS = if_else(meanPHGDHCBS < quantile(meanPHGDHCBS, probs = 0.3333), "low_PHGDHCBS",
      if_else(meanPHGDHCBS < quantile(meanPHGDHCBS, probs = 0.6666), "medium_PHGDHCBS", "high_PHGDHCBS")
    )
  )
} else if (measure == "overlap") {
  ### overlap min and overlap max
  surv_data <- dplyr::mutate(surv_data,
    ismaxPHGDH = if_else(PHGDH > median(PHGDH), "high_PHGDHCBS", "low_PHGDHCBS"),
    ismaxCBS = if_else(CBS > median(CBS), "high_PHGDHCBS", "low_PHGDHCBS"),
    stratPHGDHCBS = if_else(ismaxPHGDH == ismaxCBS, ismaxPHGDH, "indetermined")
  )
}
## Visualize stratification
ggplot(surv_data, aes(PHGDH, CBS, color = stratPHGDHCBS)) +
  geom_point()

## Kaplan Meyer of the selected measure
### OS
fit <- survfit(Surv(OS, OS_event) ~ stratPHGDHCBS,
  data = dplyr::filter(surv_data, !stratPHGDHCBS %in% c("medium_PHGDHCBS", "indetermined"))
)
print(fit)

### Change color, linetype by strata, risk.table color by strata
kmos_tcga <- ggsurvplot(fit,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw(), # Change ggplot2 theme
  palette = c("tomato3", "seagreen")
) +
  labs(title = paste0(measure, " split of PHGDH | CBS expression vs OS"))

print(kmos_tcga)
ggsave(paste0("results/survplots/FIGURES/kmos_tcga_", measure, "_phgdhcbs.svg"),
  kmos_tcga$plot,
  height = 4, width = 4
)
### PFS
fit <- survfit(Surv(PFS, PFS_event) ~ stratPHGDHCBS,
  data = dplyr::filter(surv_data, !stratPHGDHCBS %in% c("medium_PHGDHCBS", "indetermined"))
)
print(fit)

### Change color, linetype by strata, risk.table color by strata
kmpfs_tcga <- ggsurvplot(fit,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw(), # Change ggplot2 theme
  palette = c("tomato3", "seagreen")
) +
  labs(title = paste0(measure, " split of PHGDH | CBS expression vs PFS"))

print(kmpfs_tcga)
ggsave(paste0("results/survplots/FIGURES/kmpfs_tcga_", measure, "_phgdhcbs.svg"),
  kmpfs_tcga$plot,
  height = 4, width = 4
)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Somatic mutations
## Parameters
mutation_type <- "germline"
exclude_TP53 <- F
## Load data from m. Guo et al 2022 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9738094/
somatic_BRCANess_raw <- read_xlsx("data/Other/mGuo2022_TableS3.xlsx", skip = 2) %>%
  dplyr::filter(Cancer_Type == "PAAD")

germline_BRCANess_raw <- read_xlsx("data/Other/mGuo2022_TableS4.xlsx", skip = 2) %>%
  dplyr::filter(Norm_Sample_Barcode == "PAAD") # there is a erratum in the excell and this is the column name

order <- dplyr::arrange(TCGA_ISRact_projection, get("ISRactPCA"))
if (mutation_type == "somatic") {
  BRCAness <- dplyr::select(somatic_BRCANess_raw, Gene, `Sample Barcode`) %>%
    dplyr::rename(GeneID = "Gene", sample = "Sample Barcode") %>%
    dplyr::mutate(
      mutated = 1,
      sample = str_replace(sample, "^(.*-01)A.*", "\\1")
    ) %>%
    dplyr::group_by(GeneID, sample) %>%
    summarise_all(sum) %>% # there are 3 genes with multiple mutations (110 vs 107)
    pivot_wider(names_from = "sample", values_from = "mutated") %>%
    replace(is.na(.), 0) %>%
    dplyr::select(GeneID, matches(order$sample)) %>% # From 95 to 84 samples, so 11 non PAAD originally included
    dplyr::filter(if_else(exclude_TP53, GeneID != "TP53", is.character(GeneID)))
} else if (mutation_type == "germline") {
  BRCAness <- dplyr::select(germline_BRCANess_raw, Gene, Sample_Barcode) %>%
    dplyr::rename(GeneID = "Gene", sample = "Sample_Barcode") %>%
    dplyr::mutate(
      mutated = 1,
      sample = str_replace(sample, "^(.*-01)A.*", "\\1")
    ) %>%
    dplyr::group_by(GeneID, sample) %>%
    summarise_all(sum) %>% # there are no multiple mutatins (72 before and after)
    pivot_wider(names_from = "sample", values_from = "mutated") %>%
    replace(is.na(.), 0) %>%
    dplyr::select(GeneID, matches(order$sample)) %>% # From 71 to 65 samples, so 6 non PAAD originally included
    dplyr::filter(if_else(exclude_TP53, GeneID != "TP53", is.character(GeneID)))
}

## pheatmap of BRCANess genes
annot_cols <- dplyr::select(TCGA_ISRact_projection, sample, ISRactPCA, ISRactPCA_bin, PAMG) %>%
  dplyr::mutate(ISRactPCA = scale(ISRactPCA)) %>%
  dplyr::filter(sample %in% names(BRCAness))

annot_colors <- list(
  ISRactPCA = c("seagreen", "white", "tomato3"),
  ISRactPCA_bin = c(`high_ISRactPCA` = "brown", `intermediate_ISRactPCA` = "grey", `low_ISRactPCA` = "#006837"),
  PAMG = c("#FF7F00", "white", "#377DB8")
)

column_to_rownames(BRCAness, "GeneID") %>% pheatmap(
  color = c("white", "grey", "grey", "black"),
  annotation_col = column_to_rownames(annot_cols, "sample"),
  annotation_colors = annot_colors, cluster_cols = F, show_colnames = F,
  main = paste0("BRCAness ", mutation_type, " mutations")
)

#-------------------------------------------------------------------------------
