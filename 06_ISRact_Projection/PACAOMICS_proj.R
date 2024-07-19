#Import library
library(tidyverse)
library(biomaRt)
library(pdacmolgrad) #devtools::install_github("RemyNicolle/pdacmolgrad")
library(edgeR)
library(Hmisc)
library(car)
#library(ggpubr)
################################################################################
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")
source("../06_Human_Cohort_Projection/The-Molecular-Signature-Generator/R/functions.R")

# Import datasets
## PACAOMICS
### Remy Nicolle's 2017 PDX data
RN2017_raw <- read_delim("data/PACAOMICs/Human-Tumor_rawcount_Transcriptome.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)
### Extended 73 sample cohort
load("data/PACAOMICs/Alexias_53PDX/PDX_HUMAN_RAW.RData")
PACAOMICs_90_raw <- as_tibble(x, rownames ="EnsemblID")

### inherited sample info and top samples from Sauyeun_PDX
sample_info <- read_delim("data/Sauyeun_PDX/sample_info.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

top_samples <-arrange(sample_info, ICA3) %>%
  dplyr::slice( unique(c(1:5, n() - 0:4)) ) %>%
  mutate(ISRact = ifelse(ICA3 < 0, "low_ISRact", "high_ISRact")) %>%
  arrange(sample)

################################################################################
# Provisional projection function
Best_nc_gw <- function(icas_list,
                       range.comp,
                       metadata,
                       metadata_id = 0,
                       vars,
                       is.categorical = FALSE) {
  
  to_plot <- tibble(nc = character(),
                    var = character(),
                    IC = character(),
                    R = numeric(),
                    Pval = numeric())
  
  correlations <- list()
  for (cv in vars){
    correlations[[cv]] <- list()
    correlations[[cv]][["Pvals"]] <- list() #! store only Rs?
    correlations[[cv]][["Rs"]] <- list()
  }
  
  range.comp = names(icas_list$A)
  for (n.comp in range.comp) {
    ica.nc <- icas_list$S[[n.comp]] %>% as_tibble(rownames = "EnsemblID")
    test_cont <- dplyr::rename(metadata, EnsemblID = all_of(metadata_id)) %>%
      dplyr::select(EnsemblID, all_of(c(cv))) %>%
      dplyr::rename(IC_tomatch = cv) %>%
      inner_join(ica.nc, by = "EnsemblID") %>%
      drop_na() %>%
      column_to_rownames("EnsemblID")
    
    res2 <- rcorr(x=as.matrix(test_cont), type = "pearson")
    
    res2$r <- res2$r["IC_tomatch", names(ica.nc)[-1],drop=F]
    res2$P <- res2$P["IC_tomatch", names(ica.nc)[-1],drop=F]
    
    for (cv in vars) { 
      correlations[[cv]]$Rs[[n.comp]] <- res2$r["IC_tomatch",]
      correlations[[cv]]$Pvals[[n.comp]] <- res2$P["IC_tomatch",]
    }
    
    # Check the best IC and save it as a representative of the nc
    for (cv in vars){
      best_c <- sort(abs(correlations[[cv]]$Rs[[n.comp]]),decreasing = T)[1] %>% names()
      print(best_c)
      to_plot %>% add_row(nc = n.comp,
                          var = cv,
                          IC = best_c,
                          R = abs(correlations[[cv]]$Rs[[n.comp]][[best_c]]),
                          Pval = correlations[[cv]]$Pvals[[n.comp]][[best_c]]) -> to_plot
    }
  }
  
  order_bars <- group_by(to_plot, nc) %>%
    summarise(meanR = mean(R)) %>%
    arrange(desc(meanR)) %>%
    dplyr::select(nc)
  
  to_plot$nc <- to_plot$nc %>%
    factor(order_bars$nc)
  
  to_plot %>% ggplot(aes(x=nc, y=R, fill=var)) +
    geom_bar(stat='identity', position=position_dodge()) +
    geom_text(aes(x = nc, y = 0.1, label = IC), angle = 90, position = position_dodge(width = 0.9)) +
    theme_minimal() +
    ggtitle(paste0("ICA Space correlation with variable/s: ", paste(vars, sep = ", "))) +
    rotate_x_text() +
    rremove("xlab")
}
################################################################################
# PARAMETERS
PACAOMICS_raw = RN2017_raw #RN2017_raw | PACAOMICs_90_raw
norm_method = "upperquartile" #TMM | upperquartile
signature_ISRact = "ISRactIPCA"
################################################################################

# Translate EnsemblID to gene names
## Version 75 for PDX data
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)#listAttributes(ensembl75, page="feature_page")

annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id',
                                        'entrezgene',
                                        'chromosome_name'), mart = ensembl75)

translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

# Processing and Normalization
PACAOMICS_norm_ <- column_to_rownames(PACAOMICS_raw, "EnsemblID") %>%
  DGEList() %>%
  calcNormFactors(method= norm_method) %>%
  cpm(log=TRUE)

PACAOMICS_norm_Ensembl <- as_tibble(PACAOMICS_norm_, rownames = "EnsemblID")

# Create sample_info_PACAOMICS
rownames(PACAOMICS_norm_) <- translate[rownames(PACAOMICS_norm_)] #! PROBLEMATIC ASSIGNATION

type_pamg <- projectMolGrad(newexp = PACAOMICS_norm_, geneSymbols = rownames(PACAOMICS_norm_)) %>%
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


  projection_PACAOMICS <- predict(pca_pdx, PACAOMICS_norm) %>%
    as_tibble() %>%
    mutate(sample = PACAOMICS_norm$sample, .before = 1) %>% # manual assignation due to loss of sampleID in predict
    left_join(top_samples[, c("sample", "ISRact")]) %>%
    mutate(
      ISRact = replace(ISRact, is.na(ISRact) & sample %in% sample_info$sample, "intermediate_ICA3"),
      ISRact = replace_na(ISRact, "Unknown"),
      ISRact = str_replace(ISRact, "ICA3", "ISRact")
    )

  PACAOMICS_PC1 <- arrange(projection_PACAOMICS, PC1) %>% # ! wrong! remove the last part of medium or the full object
    dplyr::filter(!sample == "!") %>%
    mutate(PC1status = cut(.$PC1, breaks = c(quantile(.$PC1, c(0:3 / 3))), labels = c("low_ISRact", "intermediate_ISRact", "high_ISRact"), include.lowest = TRUE)) %>%
    dplyr::select(sample, PC1, PC1status) %>%
    left_join(top_samples[, c("sample", "ISRact")]) %>%
    mutate(
      ISRact = replace_na(ISRact, "intermediate_ICA3"),
      ISRact = str_replace(ISRact, "ICA3", "ISRact")
    ) %>%
    inner_join(sample_info_PACAOMICS[, c("sample", "PAMG", "ICA3")], by = "sample")



  # Plot pca and projections
  ## Add ISR status to PCA df
  pca_full_df <- pca_pdx[["x"]] %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    inner_join(top_samples[, c("sample", "ISRact")])
  
} else if (signature_ISRact == "ISRactICA"){
  # Import the projection
  ## ISRactPCA
  ica_pdx <- read_rds("../06_Human_Cohort_Projection/01_PDXTranslation_to_PDXTrascription/ISRactICA_IC3.RDS")
  
  # filtering
  ## XY
  non_sex <- annot_ensembl75 %>% dplyr::filter(!chromosome_name %in% c("X", "Y")) %>%
    pull(ensembl_gene_id)
  
  PACAOMICS_norm_Ensembl_ <- dplyr::filter(PACAOMICS_norm_Ensembl, EnsemblID %in% non_sex)
  
  ## mean centring gene wise
  column_to_rownames(PACAOMICS_norm_Ensembl_,"EnsemblID") %>%
    apply(1, function(x) x - mean(x)) %>% t() %>%
    data.frame() -> PACAOMICS_norm_Ensembl__
  
  ## Inter quartile range (measure variability)
  PACAOMICS_norm_Ensembl__ %>% apply(1, IQR) -> iqrs
  mostvar <- iqrs[iqrs > median(iqrs)]
  PACAOMICS_norm_Ensembl__ %>% as_tibble(rownames = "EnsemblID") %>%
    dplyr::filter(EnsemblID %in% names(mostvar)) -> PACAOMICS_icaready
  
  # Check distribution variability and number of genes of the different steps
  
  
  #Find Molecular signature
  icas_list <- Range_ICA(PACAOMICS_icaready, "EnsemblID", 2:20)
  
  Best_nc_gw(icas_list,
          range.comp = 2:20,
          metadata = as_tibble(ica_pdx$S, rownames ="EnsemblID"),
          metadata_id = "EnsemblID",
          vars = c("IC.3"),
          is.categorical = F)
  
  ISRactICA <- Range_ICA(PACAOMICS_icaready, "EnsemblID", 13) 
    # IC.10 is the component
    
  ISRactICA$A[,"IC.10"] = ISRactICA$A[,"IC.10"]*(-1)
  ISRactICA$S[,"IC.10"] = ISRactICA$S[,"IC.10"]*(-1)
  
# Create a PACAOMICS_PCA like object to explore the association in depth in the following analyses

    projection_PACAOMICS <- as_tibble(ISRactICA$A, rownames = "sample") %>%
    dplyr::arrange(IC.10) %>%
    dplyr::select(sample, IC.10, IC.2) %>%
    dplyr::rename(PC1 = "IC.10", PC2 = "IC.2") %>% 
    left_join(top_samples[, c("sample", "ISRact")]) %>%
    mutate(
      ISRact = replace(ISRact, is.na(ISRact) & sample %in% sample_info$sample, "intermediate_ICA3"),
      ISRact = replace_na(ISRact, "Unknown"),
      ISRact = str_replace(ISRact, "ICA3", "ISRact")
    )
  
  PACAOMICS_PC1 <- arrange(projection_PACAOMICS, PC1) %>% # ! wrong! remove the last part of medium or the full object
    dplyr::filter(!sample == "!") %>%
    mutate(PC1status = cut(.$PC1, breaks = c(quantile(.$PC1, c(0:3 / 3))), labels = c("low_ISRact", "intermediate_ISRact", "high_ISRact"), include.lowest = TRUE)) %>%
    dplyr::select(sample, PC1, PC1status) %>%
    left_join(top_samples[, c("sample", "ISRact")]) %>%
    mutate(
      ISRact = replace_na(ISRact, "intermediate_ICA3"),
      ISRact = str_replace(ISRact, "ICA3", "ISRact")
    ) %>%
    inner_join(sample_info_PACAOMICS[, c("sample", "PAMG", "ICA3")], by = "sample")
  
  pca_full_df <- ica_pdx$A %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    dplyr::select(sample, IC.3, IC.2) %>% 
    dplyr::rename(PC1 = "IC.3", PC2 = "IC.2") %>% 
    inner_join(top_samples[, c("sample", "ISRact")])
  
}



#-------------------------------------------------------------------------------
# Plot projecion PacaOmics and then projection Proteogenomics
show_projection <- bind_rows(as_tibble(pca_full_df) %>% mutate(dataset = "Shin et al. PDX"),
                             mutate(projection_PACAOMICS, sample = paste0(sample, "_OG"),
                                    dataset = "PaCaOmics PDX")) %>%
  dplyr::rename(ISRactPCA = "PC1")

projection_scatter_pacaomics <- dplyr::mutate(show_projection) %>%
  dplyr::filter(dataset %in% c("PaCaOmics PDX"), # "Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"
                ISRact %in% c('low_ISRact', 'high_ISRact', "intermediate_ISRact", 'Unknown')) %>% 
  ggplot(aes(x=ISRactPCA, y=PC2, color = ISRact, shape = dataset)) +
  geom_point() +
  scale_shape_manual(values = c(`Shin et al. PDX` = 16, `PaCaOmics PDX` = 17, CPTAC = 15, CCLE = 18)) +
  scale_color_manual(values = c(low_ISRact = "seagreen", high_ISRact= "tomato3", intermediate_ISRact = "grey", Unknown = "#619CFF")) + #c('low_ISRact', 'high_ISRact', 'intermediate_ISRact', 'Unknown'))unique(projection_ccle$primary_tissue))
  #ylim(-250,15) +
  theme_bw()

# ggsave(projection_scatter_pacaomics,
#        filename = "results/Figures/projection_scatter_pacaomics.svg",
#        width = 7,
#        height = 3)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Check ISR marker genes
marker_genes = c("PHGDH", "CBS", "IL18")

Gene_ISRactPCA_comparison <- as_tibble(PACAOMICS_norm_, rownames = "Gene") %>%
  pivot_longer(-Gene, names_to = "sample") %>% 
  dplyr::filter(Gene %in% marker_genes) %>% 
  inner_join(PACAOMICS_PC1, by = "sample") %>%
  dplyr::filter(PC1status != "intermediate_ISRact")

for (mgene in marker_genes) {
  
  
  mgene_plot <- dplyr::filter(Gene_ISRactPCA_comparison, Gene == mgene) %>%
    ggplot(aes(y = value, x = PC1status, fill = PC1status)) + 
    geom_violin() +
    scale_fill_manual(values = c(low_ISRact = "seagreen", high_ISRact= "tomato3")) +
    rotate_x_text(90) +
    geom_boxplot(width=0.1, fill="white")+
    labs(title = mgene)
  
  print(mgene_plot)
  # ggsave(mgene_plot,
  #        filename = paste0("results/Figures/pacaomics_",mgene,".svg"),
  #        width = 3,
  #        height = 2)
}

#-------------------------------------------------------------------------------

# Plot comparisons with Basal/Classical and ISRact
## ISR vs ISRactPCA
scatter_pacaomics_isractvsisractpca <- dplyr::select(PACAOMICS_PC1, -ISRact) %>% 
  dplyr::rename(ISRactPCA = "PC1", ISRact = "ICA3") %>% 
  correlation_plotter(data = ., col1 = "ISRact", col2 = "ISRactPCA", data_name = "PaCaOmics PDX")

# ggsave(scatter_pacaomics_isractvsisractpca,
#        filename = "results/Figures/scatter_pacaomics_isractvsisractpca.svg",
#        width = 2,
#        height = 2)

## PAMG vs ISRactPCA
scatter_pacaomics_pamgvsisractpca <- dplyr::select(PACAOMICS_PC1, -ISRact) %>% 
  dplyr::rename(ISRactPCA = "PC1", ISRact = "ICA3") %>% 
  correlation_plotter(data = ., col1 = "PAMG", col2 = "ISRactPCA", data_name = "PaCaOmics PDX")

# ggsave(scatter_pacaomics_pamgvsisractpca,
#        filename = "results/Figures/scatter_pacaomics_pamgvsisractpca.svg",
#        width = 2,
#        height = 2)

## ISRact vs PAMG
scatter_pacaomics_isractvspamg <- dplyr::select(PACAOMICS_PC1, -ISRact) %>% 
  dplyr::rename(ISRactPCA = "PC1", ISRact = "ICA3") %>% 
  correlation_plotter(data = ., col1 = "ISRact", col2 = "PAMG", data_name = "PaCaOmics PDX")

# ggsave(scatter_pacaomics_isractvspamg,
#        filename = "results/Figures/scatter_pacaomics_isractvspamg.svg",
#        width = 2,
#        height = 2)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Comparisons with DNA damage mutations
## Parameters
only_BRCAness = F


## load gene set and mutation data
BRCAness_gs <- read_tsv("data/Other/mGuo2022_BRCAness.tsv", col_names = F) %>%
  dplyr::filter(X1 != "TP53") %>% 
  deframe()
mut_PDX_raw <- read_tsv("data/PACAOMICs/DNA/somaticMutations.tsv") # keep the original file for ref
mut_PDX <- dplyr::select(mut_PDX_raw, Gene.refGene, Sample_Name) %>%
  dplyr::rename(GeneID = "Gene.refGene", sample = "Sample_Name") %>%
  dplyr::mutate(mutated = 1) %>% 
  dplyr::distinct() %>% #there are like 40 genes mutated in multiple sites (2,385 vs 2,342)
  pivot_wider(names_from = "sample", values_from = "mutated") %>% 
  replace(is.na(.), 0)

order = dplyr::arrange(PACAOMICS_PC1, get("ICA3")) #PC1, ICA3
mut_PDX_BRCAness <- dplyr::select(mut_PDX, GeneID, matches(order$sample)) %>%
  dplyr::filter(GeneID %in% BRCAness_gs)

##pheatmap of BRCANess genes
annot_cols =  dplyr::select(PACAOMICS_PC1, sample, ICA3, ISRact, PC1, PAMG) %>%
  dplyr::rename(ISRactPCA = "PC1", ISRact_binned = "ISRact", ISRact = "ICA3")

annot_colors = list(ISRactPCA = c("seagreen", "white", "tomato3"),
                    ISRact_binned = c(`high_ISRact` = "brown", `intermediate_ISRact` = "grey", `low_ISRact` = "#006837"),
                    PAMG = c("#FF7F00", "white", "#377DB8"),
                    ISRact  = c("#FFFFCC", "#006837"))
                    
column_to_rownames(mut_PDX_BRCAness, "GeneID") %>% pheatmap(color = c("white", "black"), annotation_col = column_to_rownames(annot_cols, "sample"),
                                            annotation_colors = annot_colors,cluster_cols = F)


## Statistic analysis: ISractPCA any mutated vs not
### general
if (only_BRCAness){
mutated_genes <- dplyr::filter(mut_PDX, GeneID %in% BRCAness_gs) %>%
  dplyr::summarise(across(where(is.double), sum)) %>%
  pivot_longer(everything(),names_to = "sample", values_to = "mutated_genes") %>%
  inner_join(annot_cols, by = "sample")
} else {
  mutated_genes <- dplyr::summarise(mut_PDX, across(where(is.double), sum)) %>%
    pivot_longer(everything(),names_to = "sample", values_to = "mutated_genes") %>%
    inner_join(annot_cols, by = "sample")
}
### correlation
correlation_plotter(mutated_genes, "ISRact", "mutated_genes", "PaCaOmics")
correlation_plotter(mutated_genes, "ISRactPCA", "mutated_genes", "PaCaOmics")

## ttest low vs high
mutated_genes_filt <- dplyr::filter(mutated_genes, ISRact_binned != "intermediate_ISRact") %>%
  dplyr::mutate(ISRact_binned = fct(ISRact_binned, levels= c("low_ISRact", "high_ISRact")))
leveneTest(mutated_genes ~ ISRact_binned, mutated_genes_filt) # unsignificant -> we can use t test

stats <- t.test(mutated_genes ~ ISRact_binned, mutated_genes_filt, alternative = "two.sided", var.equal = T)

stats_text <- paste0("t-test: pval = ", format(stats$p.value, scientific=T))

ggplot(mutated_genes_filt, aes(x=ISRact_binned, y=mutated_genes)) +
  geom_dotplot(binaxis='y', stackdir='center') +
  stat_summary(fun.y=median, geom="point", shape=18,
               size=3, color="red") +
  #scale_y_log10() +
  labs(subtitle = stats_text) +
  ylab("BRCAness mutated genes") +
  theme_bw()
