# ISRactPCA/ISRactICA projection into PACAOMICs Cohort (NOT TESTED IN PACAOMICs_90_raw)

################################################################################
#' Import different files
#' Filter data and project ISRactPCA or ISRactICA
#' Produce projection scatterplot
#' Compare marker gene transcription
#' Explore relation with PAMG sample score
#' Explore mutation status compared with ISRactPCA/ISRactICA
################################################################################

# Import library
library(tidyverse)
library(biomaRt)
library(pdacmolgrad) # devtools::install_github("RemyNicolle/pdacmolgrad")
library(edgeR)
library(Hmisc)
library(car)
library(ggpubr)

## setwd and import local functions
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")
source("../06_Human_Cohort_Projection/The-Molecular-Signature-Generator/R/functions.R")

# Import datasets
## PACAOMICS
### Remy Nicolle's 2017 PDX data
RN2017_raw <- read_delim("data/PACAOMICs/Human-Tumor_rawcount_Transcriptome.tsv",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)
### Extended 73 sample cohort
load("data/PACAOMICs/Alexias_53PDX/PDX_HUMAN_RAW.RData")
PACAOMICs_90_raw <- as_tibble(x, rownames = "EnsemblID")

### inherited sample info and top samples from Sauyeun_PDX
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
PACAOMICS_raw <- RN2017_raw # RN2017_raw | PACAOMICs_90_raw
norm_method <- "upperquartile" # TMM | upperquartile
signature_ISRact <- "ISRactPCA" # ISRactPCA | ISRactICA
signature_bin = paste0(signature_ISRact, "_bin")
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
  "external_gene_id",
  "entrezgene",
  "chromosome_name"
), mart = ensembl75)

translate <- deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

# Processing and Normalization
PACAOMICS_norm_ <- column_to_rownames(PACAOMICS_raw, "EnsemblID") %>%
  DGEList() %>%
  calcNormFactors(method = norm_method) %>%
  cpm(log = TRUE)

PACAOMICS_norm_Ensembl <- as_tibble(PACAOMICS_norm_, rownames = "EnsemblID")

# Create sample_info_PACAOMICS
PACAOMICS_norm_GeneNames <- PACAOMICS_norm_

rownames(PACAOMICS_norm_GeneNames) <- translate[rownames(PACAOMICS_norm_GeneNames)]

type_pamg <- projectMolGrad(newexp = PACAOMICS_norm_GeneNames, geneSymbols = rownames(PACAOMICS_norm_GeneNames)) %>%
  as_tibble(rownames = "sample")

sample_info_PACAOMICS <- dplyr::select(sample_info, -PAMG) %>%
  right_join(dplyr::select(type_pamg, sample, PDX)) %>%
  dplyr::rename(PAMG = PDX)

## ISRactPCA
if (signature_ISRact == "ISRactPCA") {
  # Import the projection
  ## ISRactPCA
  pca_pdx <- read_rds("data/Classifiers/pca_pdx_ENZO.RDS")

  PACAOMICS_norm <- t(PACAOMICS_norm_) %>%
    as_tibble(rownames = "sample")

  # Projection
  filter_pca <- function(.data, objective) {
    as_tibble(.data, rownames = "tmp") %>%
      dplyr::filter(tmp %in% objective) %>%
      column_to_rownames("tmp") %>%
      data.matrix()
  }

  pca_pdx$rotation <- filter_pca(pca_pdx$rotation, names(PACAOMICS_norm)[-1])
  pca_pdx$center <- filter_pca(pca_pdx$center, names(PACAOMICS_norm)[-1])
  pca_pdx$scale <- filter_pca(pca_pdx$scale, names(PACAOMICS_norm)[-1])


  PACAOMICS_PCAspace <- predict(pca_pdx, PACAOMICS_norm) %>%
    as_tibble() %>%
    mutate(sample = PACAOMICS_norm$sample, .before = 1) %>% # manual assignation due to loss of sampleID in predict
    left_join(top_samples[, c("sample", "ISRact_bin")]) %>%
    mutate(
      ISRact_bin = replace(ISRact_bin, is.na(ISRact_bin) & sample %in% sample_info$sample, "intermediate_ISRact"),
      ISRact_bin = replace_na(ISRact_bin, "Unknown")
    )

  PACAOMICS_ISRact_projection <- arrange(PACAOMICS_PCAspace, PC1) %>%
    dplyr::filter(!sample == "!") %>%
    mutate(ISRactPCA_bin = cut(.$PC1, breaks = c(quantile(.$PC1, c(0:3 / 3))), labels = c("low_ISRactPCA", "intermediate_ISRactPCA", "high_ISRactPCA"), include.lowest = TRUE)) %>%
    dplyr::select(sample, PC1, ISRactPCA_bin, ISRact_bin) %>%
    left_join(sample_info_PACAOMICS[, c("sample", "PAMG", "ICA3")], by = "sample") %>%
    dplyr::rename(ISRact = "ICA3", ISRactPCA = "PC1")


  # Plot pca and projections
  ## Add ISR status to PCA df
  pca_training <- pca_pdx[["x"]] %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    inner_join(top_samples[, c("sample", "ISRact", "ISRact_bin")])
} else if (signature_ISRact == "ISRactICA") {
  # Import the projection
  ## ISRactPCA
  ica_pdx <- read_rds("../06_Human_Cohort_Projection/01_PDXTranslation_to_PDXTrascription/ISRactICA_IC3.RDS")

  # filtering
  ## XY
  non_sex <- annot_ensembl75 %>%
    dplyr::filter(!chromosome_name %in% c("X", "Y")) %>%
    pull(ensembl_gene_id)

  PACAOMICS_norm_Ensembl_ <- dplyr::filter(PACAOMICS_norm_Ensembl, EnsemblID %in% non_sex)

  ## mean centring gene wise
  column_to_rownames(PACAOMICS_norm_Ensembl_, "EnsemblID") %>%
    apply(1, function(x) x - mean(x)) %>%
    t() %>%
    data.frame() -> PACAOMICS_norm_Ensembl__

  ## Inter quartile range (measure variability)
  PACAOMICS_norm_Ensembl__ %>% apply(1, IQR) -> iqrs
  mostvar <- iqrs[iqrs > median(iqrs)]
  PACAOMICS_norm_Ensembl__ %>%
    as_tibble(rownames = "EnsemblID") %>%
    dplyr::filter(EnsemblID %in% names(mostvar)) -> PACAOMICS_icaready

  # Provisional projection function
  Best_nc_gw <- function(icas_list,
                         range.comp,
                         metadata,
                         metadata_id = 0,
                         vars,
                         is.categorical = FALSE) {
    to_plot <- tibble(
      nc = character(),
      var = character(),
      IC = character(),
      R = numeric(),
      Pval = numeric()
    )
    
    correlations <- list()
    for (cv in vars) {
      correlations[[cv]] <- list()
      correlations[[cv]][["Pvals"]] <- list() # ! store only Rs?
      correlations[[cv]][["Rs"]] <- list()
    }
    
    range.comp <- names(icas_list$A)
    for (n.comp in range.comp) {
      ica.nc <- icas_list$S[[n.comp]] %>% as_tibble(rownames = "EnsemblID")
      test_cont <- dplyr::rename(metadata, EnsemblID = all_of(metadata_id)) %>%
        dplyr::select(EnsemblID, all_of(c(cv))) %>%
        dplyr::rename(IC_tomatch = cv) %>%
        inner_join(ica.nc, by = "EnsemblID") %>%
        drop_na() %>%
        column_to_rownames("EnsemblID")
      
      res2 <- rcorr(x = as.matrix(test_cont), type = "pearson")
      
      res2$r <- res2$r["IC_tomatch", names(ica.nc)[-1], drop = F]
      res2$P <- res2$P["IC_tomatch", names(ica.nc)[-1], drop = F]
      
      for (cv in vars) {
        correlations[[cv]]$Rs[[n.comp]] <- res2$r["IC_tomatch", ]
        correlations[[cv]]$Pvals[[n.comp]] <- res2$P["IC_tomatch", ]
      }
      
      # Check the best IC and save it as a representative of the nc
      for (cv in vars) {
        best_c <- sort(abs(correlations[[cv]]$Rs[[n.comp]]), decreasing = T)[1] %>% names()
        print(best_c)
        to_plot %>% add_row(
          nc = n.comp,
          var = cv,
          IC = best_c,
          R = abs(correlations[[cv]]$Rs[[n.comp]][[best_c]]),
          Pval = correlations[[cv]]$Pvals[[n.comp]][[best_c]]
        ) -> to_plot
      }
    }
    
    order_bars <- group_by(to_plot, nc) %>%
      summarise(meanR = mean(R)) %>%
      arrange(desc(meanR)) %>%
      dplyr::select(nc)
    
    to_plot$nc <- to_plot$nc %>%
      factor(order_bars$nc)
    
    to_plot %>% ggplot(aes(x = nc, y = R, fill = var)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_text(aes(x = nc, y = 0.1, label = IC), angle = 90, position = position_dodge(width = 0.9)) +
      theme_minimal() +
      ggtitle(paste0("ICA Space correlation with variable/s: ", paste(vars, sep = ", "))) +
      rotate_x_text() +
      rremove("xlab")
  }

  # Find Molecular signature
  icas_list <- Range_ICA(PACAOMICS_icaready, "EnsemblID", 2:20)

  Best_nc_gw(icas_list,
    range.comp = 2:20,
    metadata = as_tibble(ica_pdx$S, rownames = "EnsemblID"),
    metadata_id = "EnsemblID",
    vars = c("IC.3"),
    is.categorical = F
  )

  ISRactICA <- Range_ICA(PACAOMICS_icaready, "EnsemblID", 13)
  # IC.10 is the component

  ISRactICA$A[, "IC.10"] <- ISRactICA$A[, "IC.10"] * (-1)
  ISRactICA$S[, "IC.10"] <- ISRactICA$S[, "IC.10"] * (-1)

  # Create a PACAOMICS_PCA like object to explore the association in depth in the following analyses

  PACAOMICS_PCAspace <- as_tibble(ISRactICA$A, rownames = "sample") %>%
    dplyr::arrange(IC.10) %>%
    dplyr::select(sample, IC.10, IC.2) %>%
    dplyr::rename(ISRactICA = "IC.10", other_comp = "IC.2") %>%
    left_join(top_samples[, c("sample", "ISRact_bin")]) %>%
    mutate(
      ISRact_bin = replace(ISRact_bin, is.na(ISRact_bin) & sample %in% sample_info$sample, "intermediate_ISRact"),
      ISRact_bin = replace_na(ISRact_bin, "Unknown")
    )
  

  PACAOMICS_ISRact_projection <- arrange(PACAOMICS_PCAspace, ISRactICA) %>% 
    dplyr::filter(!sample == "!") %>%
    mutate(ISRactICA_bin = cut(.$ISRactICA, breaks = c(quantile(.$ISRactICA, c(0:3 / 3))), labels = c("low_ISRactICA", "intermediate_ISRactICA", "high_ISRactICA"), include.lowest = TRUE)) %>%
    dplyr::select(sample, ISRactICA, ISRactICA_bin, ISRact_bin) %>%
    inner_join(sample_info_PACAOMICS[, c("sample", "PAMG", "ICA3")], by = "sample") %>% 
    dplyr::rename(ISRact = "ICA3")

  pca_training <- ica_pdx$A %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    dplyr::select(sample, IC.3, IC.2) %>%
    dplyr::rename(ISRactICA = "IC.3", other_comp = "IC.2") %>% 
    left_join(sample_info[, c("sample", "ICA3", "PAMG")]) %>% 
    left_join(PACAOMICS_ISRact_projection[, c("sample", "ISRact_bin")]) %>% 
    dplyr::rename(ISRact = "ICA3")
  }

# ! From here in PCA is good but need checking on ICA
#-------------------------------------------------------------------------------
# Plot projecion PacaOmics and then projection Proteogenomics
show_projection <- bind_rows(
  as_tibble(pca_training) %>% mutate(dataset = "Shin et al. PDX"),
  mutate(PACAOMICS_PCAspace,
    sample = paste0(sample, "_OG"),
    dataset = "PaCaOmics PDX"
  )
)

projection_scatter_pacaomics <- dplyr::mutate(show_projection, ISRact_bin = fct(ISRact_bin, levels = c("low_ISRact", "intermediate_ISRact", "high_ISRact", "Unknown"))) %>%
  dplyr::filter(
    dataset %in% c("PaCaOmics PDX"), # "Shin et al. PDX","PACAOMICS PDX", "CPTAC", "CCLE"
    ISRact_bin %in% c("low_ISRact", "high_ISRact", "intermediate_ISRact", "Unknown")
  ) %>%
  ggplot(aes(x = get(signature_ISRact), y = other_comp, color = ISRact_bin, shape = dataset)) +
  geom_point() +
  scale_shape_manual(values = c(`Shin et al. PDX` = 16, `PaCaOmics PDX` = 17, CPTAC = 15, CCLE = 18)) +
  scale_color_manual(values = c(low_ISRact = "seagreen", high_ISRact = "tomato3", intermediate_ISRact = "grey", Unknown = "#619CFF")) + # c('low_ISRact', 'high_ISRact', 'intermediate_ISRact', 'Unknown'))unique(projection_ccle$primary_tissue))
  # ylim(-250,15) +
  ylab("Other Component") +
  xlab(signature_ISRact) +
  theme_bw()

# ggsave(projection_scatter_pacaomics,
#        filename = "results/Figures/projection_scatter_pacaomics.svg",
#        width = 7,
#        height = 3)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Check ISR marker genes associated with ISRActPCA extreme samples
marker_genes <- c("PHGDH", "CBS", "IL18")

Gene_ISRact_projection_comparison <- as_tibble(PACAOMICS_norm_GeneNames, rownames = "Gene") %>%
  pivot_longer(-Gene, names_to = "sample") %>%
  dplyr::filter(Gene %in% marker_genes) %>%
  inner_join(PACAOMICS_ISRact_projection, by = "sample") %>%
  dplyr::filter(!str_detect(get(signature_bin), "intermediate_ISRact"))

for (mgene in marker_genes) {
  mgene_plot <- dplyr::filter(Gene_ISRact_projection_comparison, Gene == mgene) %>%
    ggplot(aes(y = value, x = get(signature_bin), fill = get(signature_bin))) +
    geom_violin() +
    scale_fill_manual(values = c("seagreen", "tomato3")) +
    rotate_x_text(90) +
    geom_boxplot(width = 0.1, fill = "white") +
    labs(title = mgene) +
    rremove("xlab") +
    rremove("legend.title")

  print(mgene_plot)
  # ggsave(mgene_plot,
  #        filename = paste0("results/Figures/pacaomics_",mgene,".svg"),
  #        width = 3,
  #        height = 2)
}

#-------------------------------------------------------------------------------

# Plot comparisons with Basal/Classical and ISRact
## ISR vs ISRactPCA
scatter_pacaomics_isractvsisractpca <- correlation_plotter(data = PACAOMICS_ISRact_projection, col1 = "ISRact", col2 = signature_ISRact, data_name = "PaCaOmics PDX")

# ggsave(scatter_pacaomics_isractvsisractpca,
#        filename = "results/Figures/scatter_pacaomics_isractvsisractpca.svg",
#        width = 2,
#        height = 2)

## PAMG vs ISRactPCA
scatter_pacaomics_pamgvsisractpca <- correlation_plotter(data = PACAOMICS_ISRact_projection, col1 = "PAMG", col2 = signature_bin, data_name = "PaCaOmics PDX")

# ggsave(scatter_pacaomics_pamgvsisractpca,
#        filename = "results/Figures/scatter_pacaomics_pamgvsisractpca.svg",
#        width = 2,
#        height = 2)

## ISRact vs PAMG
scatter_pacaomics_isractvspamg <- correlation_plotter(data = PACAOMICS_ISRact_projection, col1 = "ISRact", col2 = "PAMG", data_name = "PaCaOmics PDX")

# ggsave(scatter_pacaomics_isractvspamg,
#        filename = "results/Figures/scatter_pacaomics_isractvspamg.svg",
#        width = 2,
#        height = 2)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Comparisons with DNA damage mutations
## Parameters
only_BRCAness <- F

## load gene set and mutation data
BRCAness_gs <- read_tsv("data/Other/mGuo2022_BRCAness.tsv", col_names = F) %>%
  dplyr::filter(X1 != "TP53") %>%
  deframe()
mut_PDX_raw <- read_tsv("data/PACAOMICs/DNA/somaticMutations.tsv") # keep the original file for ref
mut_PDX <- dplyr::select(mut_PDX_raw, Gene.refGene, Sample_Name) %>%
  dplyr::rename(GeneID = "Gene.refGene", sample = "Sample_Name") %>%
  dplyr::mutate(mutated = 1) %>%
  dplyr::distinct() %>% # there are like 40 genes mutated in multiple sites (2,385 vs 2,342)
  pivot_wider(names_from = "sample", values_from = "mutated") %>%
  replace(is.na(.), 0)

order <- dplyr::arrange(PACAOMICS_ISRact_projection, ISRact) # PC1, ICA3
mut_PDX_BRCAness <- dplyr::select(mut_PDX, GeneID, matches(order$sample)) %>%
  dplyr::filter(GeneID %in% BRCAness_gs)

## pheatmap of BRCANess genes
annot_cols <- dplyr::select(PACAOMICS_ISRact_projection, sample, ISRact, ISRact_bin, matches(signature_ISRact), matches(signature_bin), PAMG)

annot_colors <- list(
  ISRactPCA = c("seagreen", "white", "tomato3"),
  ISRactICA = c("seagreen", "white", "tomato3"),
  ISRactPCA_bin = c(`high_ISRactPCA` = "brown", `intermediate_ISRactPCA` = "grey", `low_ISRactPCA` = "#006837"),
  ISRactICA_bin = c(`high_ISRactICA` = "brown", `intermediate_ISRactICA` = "grey", `low_ISRactICA` = "#006837"),
  ISRact_bin = c(`high_ISRact` = "brown", `intermediate_ISRact` = "grey", `low_ISRact` = "#006837", Unknown = "#619CFF"),
  PAMG = c("#FF7F00", "white", "#377DB8"),
  ISRact = c("#FFFFCC", "#006837")
)

column_to_rownames(mut_PDX_BRCAness, "GeneID") %>% pheatmap(
  color = c("white", "black"), annotation_col = column_to_rownames(annot_cols, "sample"),
  annotation_colors = annot_colors, cluster_cols = F
)


## Statistic analysis: ISRactPCA any mutated vs not
### general
if (only_BRCAness) {
  mutated_genes <- dplyr::filter(mut_PDX, GeneID %in% BRCAness_gs) %>%
    dplyr::summarise(across(where(is.double), sum)) %>%
    pivot_longer(everything(), names_to = "sample", values_to = "mutated_genes") %>%
    inner_join(annot_cols, by = "sample")
} else {
  mutated_genes <- dplyr::summarise(mut_PDX, across(where(is.double), sum)) %>%
    pivot_longer(everything(), names_to = "sample", values_to = "mutated_genes") %>%
    inner_join(annot_cols, by = "sample")
}
### correlation
correlation_plotter(mutated_genes, "ISRact", "mutated_genes", "PaCaOmics")
correlation_plotter(mutated_genes, signature_ISRact, "mutated_genes", "PaCaOmics")

## ttest low vs high
mutated_genes_filt <- dplyr::filter(mutated_genes, !str_detect(get(signature_bin), "intermediate_ISRact"))
leveneTest(mutated_genes ~ get(signature_bin), mutated_genes_filt) # unsignificant -> we can use t test

stats <- t.test(mutated_genes ~ get(signature_bin), mutated_genes_filt, alternative = "two.sided", var.equal = T)

stats_text <- paste0("t-test: pval = ", format(stats$p.value, scientific = T))

ggplot(mutated_genes_filt, aes(x = get(signature_bin), y = mutated_genes)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  stat_summary(
    fun.y = median, geom = "point", shape = 18,
    size = 3, color = "red"
  ) +
  # scale_y_log10() +
  labs(subtitle = stats_text) +
  ylab("BRCAness mutated genes") +
  rremove("xlab") +
  theme_bw()
#-------------------------------------------------------------------------------
