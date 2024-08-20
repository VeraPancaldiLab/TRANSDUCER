# ISRactPCA projection into ICGC cohort
################################################################################
#' Import different files
#' Filter data and project ISRactPCA
#' Produce projection scatterplot
#' Compare ISRactPCA with PAMG score
#' ISRactPCA survival
#' PHGDH/CBS survival
################################################################################

# Import libraries
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
library("AnnotationDbi")
library("illuminaHumanv4.db")

################################################################################
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")

# Import datasets
## ICGC
### Normalised microarray intensities
ICGC <- read_rds("data/ICGC/ICGC_AU.rds") %>% as_tibble(rownames = "Gene")

### Metadata
ICGC_clinical <- read_rds("data/ICGC/ICGC_AU_survival.rds")

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
################################################################################

# Translate EnsemblID to gene names
## Using microarray data
probe_info <- AnnotationDbi::select(illuminaHumanv4.db,
  ICGC$Gene,
  c("SYMBOL", "ENSEMBL", "PROBEID"),
  keytype = "SYMBOL"
) %>%
  dplyr::rename(
    Gene = SYMBOL,
    EnsemblID = ENSEMBL
  )


# load annotation with Biomart
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
#
# #listAttributes(ensembl75, page="feature_page")
# annot_ensembl <- getBM(attributes = c('ensembl_gene_id',
#                                         'external_gene_name'), mart = ensembl)
#
# gene_to_ensembl = deframe(annot_ensembl[c( "external_gene_name", "ensembl_gene_id")])
# ensembl_to_gene = setNames(names(gene_to_ensembl), gene_to_ensembl)

# Preprocessing
## Deal with Duplicated EnsemblIDs and Gene names
ICGC_ensembl <- left_join(ICGC, probe_info, by = "Gene") %>%
  dplyr::select(-c("Gene", "PROBEID")) %>%
  group_by(EnsemblID) %>%
  summarise_all(mean)

ICGC_gene <- ICGC %>%
  group_by(Gene) %>%
  summarise_all(mean)

## Calculate PAMG
type_pamg <- projectMolGrad(newexp = column_to_rownames(ICGC_gene, "Gene"), geneSymbols = ICGC_gene$Gene) %>%
  as_tibble(rownames = "sample")

sample_info_ICGC <- dplyr::select(type_pamg, sample, ICGCarray) %>%
  dplyr::rename(PAMG = ICGCarray) %>%
  left_join(rownames_to_column(ICGC_clinical, "sample"), "sample")

# Projecting datasets on this PCA
## ICGC
### transpose human data for projection
ICGC_toproj <- ICGC_ensembl %>%
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
pca_pdx$rotation <- filter_pca(pca_pdx$rotation, names(ICGC_toproj)[-1])
pca_pdx$center <- filter_pca(pca_pdx$center, names(ICGC_toproj)[-1])
pca_pdx$scale <- filter_pca(pca_pdx$scale, names(ICGC_toproj)[-1])

ICGC_PCAspace <- predict(pca_pdx, ICGC_toproj) %>%
  as_tibble(rownames = "sample") %>%
  dplyr::rename(ISRactPCA = "PC1")


CPTAC_ISRact_projection <- arrange(ICGC_PCAspace, ISRactPCA) %>%
  mutate(ISRactPCA_bin = if_else(ISRactPCA < quantile(ICGC_PCAspace$ISRactPCA, probs = 0.3333), "low_ISRactPCA",
    if_else(ISRactPCA < quantile(ICGC_PCAspace$ISRactPCA, probs = 0.6666), "intermediate_ISRactPCA", "high_ISRactPCA")
  )) %>%
  dplyr::select(sample, ISRactPCA, ISRactPCA_bin) %>%
  left_join(sample_info_ICGC, by = "sample")

# Plot pca and projections
# Add ISR status to PCA df
pca_training <- pca_pdx[["x"]] %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  dplyr::rename(ISRactPCA = "PC1") %>%
  inner_join(top_samples[, c("sample", "ISRact_bin")])


#-------------------------------------------------------------------------------
# Plot projecion PacaOmics and then projection Proteogenomics

show_projection <- bind_rows(
  as_tibble(pca_training) %>% mutate(dataset = "Sauyeun PDX"),
  as_tibble(ICGC_PCAspace) %>% mutate(
    ISRact_bin = "Unknown",
    dataset = "ICGC"
  )
)


dplyr::filter(
  show_projection, dataset %in% c("Sauyeun PDX", "PACAOMICS PDX", "ICGC"), # "Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"
  ISRact_bin %in% c("low_ISRact", "high_ISRact", "intermediate_ISRact", "Unknown")
) %>%
  ggplot(aes(x = ISRactPCA, y = PC2, color = ISRact_bin, shape = dataset)) +
  geom_point() +
  scale_shape_discrete(limits = c("Sauyeun PDX", "PACAOMICS PDX", "ICGC")) +
  scale_color_discrete(limits = c("low_ISRact", "high_ISRact", "intermediate_ISRact", "Unknown")) + # c('low_ISRact', 'high_ISRact', 'intermediate_ISRact', 'Unknown'))unique(projection_ccle$primary_tissue))
  geom_rect(
    xmin = min(ICGC_PCAspace$ISRactPCA),
    xmax = quantile(ICGC_PCAspace$ISRactPCA, probs = 0.3333),
    ymin = min(ICGC_PCAspace$PC2), ymax = max(ICGC_PCAspace$PC2), linewidth = 0, fill = "red", alpha = 0.002
  ) +
  geom_rect(
    xmin = quantile(ICGC_PCAspace$ISRactPCA, probs = 0.6666),
    xmax = max(ICGC_PCAspace$ISRactPCA),
    ymin = min(ICGC_PCAspace$PC2), ymax = max(ICGC_PCAspace$PC2), linewidth = 0, fill = "green", alpha = 0.002
  )
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Plot comparisons with Basal/Classical and ISRact
## PAMG vs ISRactPCA
correlation_plotter(data = CPTAC_ISRact_projection, col1 = "PAMG", col2 = "ISRactPCA", data_name = "ICGC")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Survival curves regarding ISR status
surv_data <- dplyr::rename(CPTAC_ISRact_projection, case_id = sample)
## OS
### Kaplan Meyer 33up vs 33down
fit <- survfit(Surv(OS, OS_event) ~ ISRactPCA_bin,
  data = dplyr::filter(surv_data, ISRactPCA_bin != "intermediate_ISRactPCA")
)
print(fit)

#### Change color, linetype by strata, risk.table color by strata
kmos_icgc_ISRact <- ggsurvplot(fit,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw(), # Change ggplot2 theme
  palette = c("tomato3", "seagreen")
)

print(kmos_icgc_ISRact)
ggsave("results/survplots/FIGURES/kmos_icgc_ISRact.svg",
  kmos_icgc_ISRact$plot,
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
kmpfs_icgc_ISRact <- ggsurvplot(fit,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw(), # Change ggplot2 theme
  palette = c("tomato3", "seagreen")
)

print(kmpfs_icgc_ISRact)
ggsave("results/survplots/FIGURES/kmpfs_icgc_ISRact.svg",
  kmpfs_icgc_ISRact$plot,
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
surv_data <- as_tibble(ICGC_toproj, rownames = "sample") %>%
  dplyr::select(sample, c(PHGDH = "ENSG00000092621", CBS = "ENSG00000160200")) %>%
  inner_join(CPTAC_ISRact_projection, by = "sample")

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
kmos_icgc <- ggsurvplot(fit,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw(), # Change ggplot2 theme
  palette = c("tomato3", "seagreen")
) +
  labs(title = paste0(measure, " split of PHGDH | CBS expression vs OS"))

print(kmos_icgc)
ggsave(paste0("results/survplots/FIGURES/kmos_icgc_", measure, "_phgdhcbs.svg"),
  kmos_icgc$plot,
  height = 4, width = 4
)
### PFS
fit <- survfit(Surv(PFS, PFS_event) ~ stratPHGDHCBS,
  data = dplyr::filter(surv_data, !stratPHGDHCBS %in% c("medium_PHGDHCBS", "indetermined"))
)
print(fit)

### Change color, linetype by strata, risk.table color by strata
kmpfs_icgc <- ggsurvplot(fit,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw(), # Change ggplot2 theme
  palette = c("tomato3", "seagreen")
) +
  labs(title = paste0(measure, " split of PHGDH | CBS expression vs PFS"))

print(kmpfs_icgc)
ggsave(paste0("results/survplots/FIGURES/kmpfs_icgc_", measure, "_phgdhcbs.svg"),
  kmpfs_icgc$plot,
  height = 4, width = 4
)
#-------------------------------------------------------------------------------
