# ISRactPCA generation and ISRactPCA/ICA gene weight analysis

################################################################################
#' Import different files
#' Filter data to produce ISRactPCA and export it
#' Predict ISRActPCA on the full cohort
#' Produce several plots fort thesis figures and beyond
################################################################################

# Import libraries
library(tidyverse)
library(readxl)
library(biomaRt)
library(pdacmolgrad) #devtools::install_github("RemyNicolle/pdacmolgrad")

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
library(pheatmap)

################################################################################
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")
source("src/formatted_cors.R")
source('src/plot_bar_enrich.R')

# Import datasets
## CPTAC (Proteogenomics)
### upper quartile normalized and log2 data 
CPTAC_tumor_already_norm <- read_delim("data/PDAC_LinkedOmics_Data/mRNA_RSEM_UQ_log2_Tumor.cct", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE) %>%
  dplyr::rename(Gene = ...1)

### Raw counts obtained by email
CPTAC_tumor_raw <- read_delim("data/PDAC_LinkedOmics_Data/CPTAC-PDAC-raw-RSEM-expected-counts-gene-level.txt", 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(idx, matches("tumor")) %>%
  dplyr::rename_all(~str_remove(., '_tumor')) %>%
  dplyr::rename(Gene = idx)


### Metadata
low_purity_samples <- read_delim("data/PDAC_LinkedOmics_Data/low_purity_samples.csv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

clinical_data <- read_excel("data/PDAC_LinkedOmics_Data/mmc1.xlsx", 
                            sheet = "Clinical_data") %>%
  mutate(follow_up_days = as.numeric(follow_up_days)) %>%
  mutate(status = ifelse(vital_status == "Deceased", 2, 1))


Molecular_phenotype_data <- read_excel("data/PDAC_LinkedOmics_Data/mmc1.xlsx", 
                                       sheet = "Molecular_phenotype_data") %>% 
  mutate_at(vars(immune_deconv:`necrosis_(%OF_TUMOR_WITH_NECROSIS)_histology_estimate`, KRAS_VAF), as.numeric) %>%
  mutate_at(vars(Bailey:Moffitt), as.factor) %>% #change type to avoid errors
  dplyr::rename(sample = case_id)

## ISRactPCA projection 
pca_pdx <- read_rds("data/Classifiers/pca_pdx_ENZO.RDS")
### inherited sample info and top samples from Shin PDX training data
sample_info <- read_delim("data/Sauyeun_PDX/sample_info.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

top_samples <- arrange(sample_info, ICA3) %>%
  dplyr::slice( unique(c(1:5, n() - 0:4)) ) %>%
  mutate(ISRact_bin = ifelse(ICA3 < 0, "low_ISRact", "high_ISRact")) %>%
  dplyr::rename(ISRact = "ICA3") %>%
  arrange(sample)


################################################################################
# PARAMETERS
norm_method = "upperquartile" #TMM | upperquartile | upperquartile_ogdata
#Variables for filter function
filter = FALSE # exclude just the low purity if FALSE
estimation_method = "histology" #histology | deconvolution
neoplastic_min = 0.2
acinar_max = 0.05
islet_max = 1
tumor_tissue_min = 0.8
################################################################################

# Choose raw / already normalized data
if (norm_method == "upperquartile_ogdata"){
  CPTAC_tumor__ <- CPTAC_tumor_already_norm
} else {
  CPTAC_tumor__ <- CPTAC_tumor_raw
}
  
# Translate EnsemblID to gene names
## Version 95 for Proteogenimics data
ensembl95 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl")#,
                        #version = 97)#listAttributes(ensembl95, page="feature_page") # 95, but 05/2024 it doesnt work? 07/24 97 doesnt work neither?

annot_ensembl95 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_name'), mart = ensembl95)


#Add a Ensemble ID column to CPTAC_tumor__ (Ensembl release 95 according to the paper)
#https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=genes&hgta_track=mane&hgta_table=mane&hgta_doSchema=describe+table+schema

reverse_translate = deframe(annot_ensembl95[c("external_gene_name", "ensembl_gene_id")])

# Preprocessing
## filtering
if (filter == TRUE){
  CPTAC_tumor_ <- dataFilter(CPTAC_tumor__, clinical_data, Molecular_phenotype_data, estimation_method, neoplastic_min, acinar_max, islet_max, tumor_tissue_min)
} else {
  CPTAC_tumor_ <- dplyr::select(CPTAC_tumor__, -low_purity_samples$case_id)
}

## normalization and filtering
if (norm_method == "upperquartile_ogdata"){
  CPTAC_tumor <- CPTAC_tumor_
} else {
  CPTAC_tumor <- column_to_rownames(CPTAC_tumor_, "Gene") %>% DGEList() %>%
    calcNormFactors(method= norm_method) %>%
    cpm(log=TRUE) %>% as_tibble(rownames = "Gene")
}

## translate gene names
CPTAC_tumor$EnsemblID <- CPTAC_tumor$Gene %>%
  reverse_translate[.] %>%
  make.names(unique = TRUE)


## Create sample_info_CPTAC
type_pamg <- projectMolGrad(newexp = column_to_rownames(dplyr::select(CPTAC_tumor,-EnsemblID), "Gene"),  geneSymbols = CPTAC_tumor$Gene) %>%
  as_tibble(rownames = "sample")

sample_info_CPTAC <- dplyr::select(type_pamg, sample, ICGCrnaseq) %>% 
  dplyr::rename(PAMG = ICGCrnaseq) %>% 
  left_join(dplyr::select(Molecular_phenotype_data, sample, KRAS_VAF), "sample")

# Projecting CPTAC on ISRactPCA
#### transpose normalised data
CPTAC_tumor_toproj <- CPTAC_tumor %>%
  pivot_longer(cols = 2:(length(CPTAC_tumor)-1), names_to = "case_id", values_to = "Expression") %>% 
  dplyr::filter(!str_detect(EnsemblID, 'NA')) %>% #remove eventual NA for gene Ensembl
  dplyr::select(-Gene) %>% 
  pivot_wider(names_from = EnsemblID, values_from = Expression, names_repair = "minimal") %>%
  column_to_rownames("case_id") 

### remove missing genes from the PCA before projection
filter_pca <- function(.data, objective){
  as_tibble(.data, rownames = "tmp") %>% 
    dplyr::filter(tmp %in% objective) %>% 
    column_to_rownames("tmp") %>%
    data.matrix()
}
pca_pdx$rotation <- filter_pca(pca_pdx$rotation, names(CPTAC_tumor_toproj)[-1])
pca_pdx$center <- filter_pca(pca_pdx$center, names(CPTAC_tumor_toproj)[-1])
pca_pdx$scale <- filter_pca(pca_pdx$scale, names(CPTAC_tumor_toproj)[-1])

CPTAC_PCAspace <-  predict(pca_pdx, CPTAC_tumor_toproj) %>% 
  as_tibble(rownames = "sample") %>% 
  left_join(sample_info_CPTAC, by = "sample") %>%
  dplyr::rename(ISRactPCA = "PC1") %>% 
  mutate(
    ISRactPCA_bin = if_else(ISRactPCA < quantile(ISRactPCA, probs = 0.3333), "low_ISRactPCA",
                               if_else(ISRactPCA < quantile(ISRactPCA, probs = 0.6666), "intermediate_ISRactPCA", "high_ISRactPCA")),
    ISRactPCA_bin = fct(ISRactPCA_bin, levels = c("low_ISRactPCA", "intermediate_ISRactPCA", "high_ISRactPCA"))
  )

CPTAC_ISRact_projection <- arrange(CPTAC_PCAspace, ISRactPCA) %>%
  dplyr::select(sample, ISRactPCA,  ISRactPCA_bin, PAMG, KRAS_VAF)

#Plot pca and projections
#Add ISR status to PCA df
pca_training <- pca_pdx[["x"]] %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  dplyr::rename(ISRactPCA = "PC1") %>% 
  inner_join(top_samples[,c("sample","ISRact_bin")])

#-------------------------------------------------------------------------------
# Plot projecion PacaOmics and then projection Proteogenomics

show_projection <- bind_rows(as_tibble(pca_training) %>% mutate(dataset = "Shin et al. PDX"),
                             as_tibble(CPTAC_PCAspace) %>% mutate(ISRact_bin = "Unknown",
                                                                    dataset = "CPTAC"))

projection_scatter_cptac <- dplyr::filter(show_projection, dataset %in% c("Shin et al. PDX","PACAOMICS PDX", "CPTAC", "CCLE"), # "Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"
                ISRact_bin %in% c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) %>% 
  ggplot(aes(x=ISRactPCA, y=PC2, color = ISRact_bin, shape = dataset)) +
  geom_point() +
  scale_shape_discrete(limits = c("Shin et al. PDX", "PACAOMICS PDX", "CPTAC", "CCLE")) +
  scale_color_discrete(c(low_ISRact = "seagreen", high_ISRact= "tomato3", 'Unknown')) + #c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown'))unique(projection_ccle$primary_tissue))
  #ylim(-250,15) +
  theme_bw()

# ggsave(projection_scatter_cptac,
#        filename = "results/Figures/projection_scatter_cptac.svg",
#        width = 7,
#        height = 3)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Check ISR marker genes
marker_genes = c("PHGDH", "CBS", "IL18")

Gene_ISRactPCA_comparison <- dplyr::select(CPTAC_tumor, -EnsemblID) %>%
  pivot_longer(-Gene, names_to = "sample") %>% 
  dplyr::filter(Gene %in% marker_genes) %>% 
  inner_join(CPTAC_ISRact_projection, by = "sample") %>%
  dplyr::filter(ISRactPCA_bin != "intermediate_ISRactPCA")

for (mgene in marker_genes) {
  

  mgene_plot <- dplyr::filter(Gene_ISRactPCA_comparison, Gene == mgene) %>%
    ggplot(aes(y = value, x = ISRactPCA_bin, fill = ISRactPCA_bin)) + 
    geom_violin() +
    scale_fill_manual(values = c(low_ISRactPCA = "seagreen", high_ISRactPCA= "tomato3")) +
    rotate_x_text(90) +
    geom_boxplot(width=0.1, fill="white")+
    labs(title = mgene) +
    rremove("xlab") +
    rremove("legend.title")
  
  # ggsave(mgene_plot,
  #        filename = paste0("results/Figures/cptac_",mgene,".svg"),
  #        width = 3,
  #        height = 2)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Check relation with original PDX samples
## Preparatives
PDX_original <- read_tsv("data/Sauyeun_PDX/PDX_PCA_input_data.tsv")

similarity_df <- inner_join(PDX_original ,dplyr::select(CPTAC_tumor, - Gene)) %>% 
  drop_na()

### Annotation
annot_ref = left_join(sample_info, dplyr::select(top_samples, sample, ISRact_bin)) %>%
  dplyr::select(sample, ISRact_bin, PAMG) %>% 
  dplyr::mutate(ISRact_bin = if_else(is.na(ISRact_bin), "intermediate_ISRact", ISRact_bin),
                df = "Shin et al. PDX")

annot_proj = dplyr::mutate(CPTAC_ISRact_projection, 
                           ISRact_bin = str_remove(ISRactPCA_bin, "PCA"),
                           df = "CPTAC") %>%
  dplyr::select(sample, ISRact_bin, PAMG, df)
annot = bind_rows(annot_ref, annot_proj) %>% column_to_rownames("sample")

annot_colors <- list(ISRact_bin = c(`high_ISRact` = "brown", `intermediate_ISRact` = "grey", `low_ISRact` = "#006837"),
                     PAMG = c("#FF7F00", "white", "#377DB8"),
                     df = c(`Shin et al. PDX` = "grey", CPTAC ="grey2")) 

## Genome wide
gw_corr_hm <- formatted_cors(dplyr::select(similarity_df,-EnsemblID), "pearson", 0) %>% 
  dplyr::select(measure1, measure2, r) %>% 
  pivot_wider(values_from = r, names_from = measure1, id_cols = measure2)

gw_corr_hm %>% column_to_rownames("measure2") %>% 
  pheatmap(scale = "none", color = colorRampPalette(c("#FBFEF9", "#A63446"))(100),
           annotation_row = annot, annotation_col = annot, annotation_colors = annot_colors,
           cluster_rows = T, cluster_cols = T,
           show_colnames = FALSE, show_rownames = FALSE) 

## Subset
subset_corr_hm <- dplyr::filter(similarity_df, EnsemblID %in% rownames(pca_pdx$rotation)) %>%
  dplyr::select(-EnsemblID) %>% 
  formatted_cors("pearson", 0) %>% 
  dplyr::select(measure1, measure2, r) %>%
  pivot_wider(values_from = r, names_from = measure1, id_cols = measure2)

subset_corr_hm %>% column_to_rownames("measure2") %>% 
  pheatmap(scale = "none", color = colorRampPalette(c("#FBFEF9", "#A63446"))(100),
           annotation_row = annot, annotation_col = annot, annotation_colors = annot_colors,
           cluster_rows = T, cluster_cols = T,
           show_colnames = FALSE, show_rownames = FALSE)  

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Plot comparisons with Basal/Classical and ISRact
## PAMG vs PC1
scatter_cptac_isractpcavspamg <- correlation_plotter(data = CPTAC_ISRact_projection, col1 = "PAMG", col2 = "ISRactPCA", data_name = "CPTAC")

# ggsave(scatter_cptac_isractpcavspamg,
#        filename = "results/Figures/scatter_cptac_isractpcavspamg.svg",
#        width = 2,
#        height = 2)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Proteomics DPEA
## load already normalized log transformed and centered data
CPTAC_comparison_metadata <- dplyr::filter(CPTAC_ISRact_projection, ISRactPCA_bin != "intermediate_ISRactPCA")
CPTAC_prot <- read_delim("data/PDAC_LinkedOmics_Data/proteomics_gene_level_MD_abundance_tumor.cct", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE) %>%
  dplyr::select(...1, CPTAC_comparison_metadata$sample) %>%
  dplyr::relocate(...1, CPTAC_comparison_metadata$sample) %>%
  column_to_rownames("...1") %>%
  na.omit()


### Boxplot of intensities distribution
rownames_to_column(CPTAC_prot, "Gene") %>%
  pivot_longer(-Gene, names_to = "samples") %>% 
  ggplot(aes(y = value, x = samples)) +
  geom_boxplot() +
  rotate_x_text(45)

## Define contrast and fit model
cond <- CPTAC_comparison_metadata$ISRactPCA_bin

design = model.matrix(~0+cond)
colnames(design) = gsub("cond","",colnames(design))

x <- c("low_ISRactPCA-high_ISRactPCA")
contrast =  makeContrasts(contrasts=x,levels=design)
fit1 <- lmFit(CPTAC_prot, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2)

### add PTSMs info to fit 3
count_columns = seq(16,34,2)
psm.count.table = read_tsv("data/PDAC_LinkedOmics_Data/CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Proteome.summary.tsv") %>%
  column_to_rownames("Gene")
fit3$count = psm.count.table[rownames(fit3$coefficients),"Distinct Peptides"]
fit4 = spectraCounteBayes(fit3)


## Visualize fit curve
### n=30 limits the boxplot to show only proteins quantified by <= 30 PSMs.
VarianceBoxplot(fit4,n=30,main="CPTAC TMT proteomics",xlab="PSM count")
VarianceScatterplot(fit4,main="CPTAC TMT proteomics")

## Analysis
DEqMS.results = outputResult(fit4,coef_col = 1)

### Export data
# write.table(DEqMS.results,"data/",sep = "\t",
#             row.names = F,quote=F)


### Use ggplot2 to plot volcano
highlight = c("PHGDH","IL18", "IDO1",
              "APP", "LUM")

DEqMS.results$log.sca.pval = -log10(DEqMS.results$sca.P.Value)
dpea_isractpca_volcano <- ggplot(DEqMS.results, aes(x = logFC, y =log.sca.pval )) + 
  geom_point(size=0.5 )+
  theme_bw(base_size = 16) + # change theme
  xlab(expression("log2(ISRactPCA high/low)")) + # x-axis label
  ylab(expression(" -log10(P-value)")) + # y-axis label
  geom_vline(xintercept = c(-1,1), colour = "red") + # Add fold change cutoffs
  geom_hline(yintercept = 3, colour = "red") + # Add significance cutoffs
  geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
  scale_colour_gradient(low = "black", high = "black", guide = FALSE)+
  geom_text_repel(data=subset(DEqMS.results, abs(logFC)>0.5&log.sca.pval > 3),
                  aes( logFC, log.sca.pval ,label=gene))+# add gene label
  geom_point(data= DEqMS.results[rownames(DEqMS.results) %in% highlight,],
             aes( logFC, log.sca.pval), color="red") + 
  geom_label_repel(data= DEqMS.results[rownames(DEqMS.results) %in% highlight,],
                   aes( logFC, log.sca.pval ,label=gene)) # add gene label

print(dpea_isractpca_volcano)
# ggsave(dpea_isractpca_volcano,
#        filename = "results/Figures/dpea_isractpca_volcano.svg",
#        width = 5,
#        height = 5)

### GSEA of fold changes
all_genesets <- msigdbr("Homo sapiens")
#### reactome
set.seed(42)
reactome <- all_genesets %>% filter(gs_subcat %in% c("CP:REACTOME"))
reactome_list = split(x = reactome$gene_symbol, f = reactome$gs_name)
fgseaReactome <- fgsea(pathways = reactome_list, 
                  stats    = deframe(rownames_to_column(DEqMS.results["logFC"])),
                  minSize  = 15,
                  maxSize  = 500)

##### Plot
topPathwaysUp <- fgseaReactome[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaReactome[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(reactome_list[topPathways], deframe(rownames_to_column(DEqMS.results["logFC"])), fgseaReactome, 
              gseaParam=0.5)

barplot_DEP_reactome <- dplyr::filter(fgseaReactome, pathway %in% topPathways) %>%
  arrange(NES) %>% 
  mutate(pathway = str_remove(pathway,"REACTOME_"),
         pathway = as_factor(pathway)) %>% 
  #dplyr::filter(p.adjust < p.adjust_th) %>%
  #slice_tail(n=20) %>%
  ggplot(aes(x = NES, y = pathway, fill = padj)) +
  geom_bar(stat="identity") +
  scale_fill_gradientn(colours = c("red", "white"), guide=guide_colourbar(reverse = TRUE))+
  theme_bw() +
  labs(title = "REACTOME") +
  ylab("")

print(barplot_DEP_reactome)
# ggsave(barplot_DEP_reactome,
#        filename = "results/Figures/barplot_DEP_reactome.svg",
#        width = 8,
#        height = 5)

#### Kegg
kegg <- all_genesets %>% filter(gs_subcat %in% c("CP:KEGG"))
kegg_list = split(x = kegg$gene_symbol, f = kegg$gs_name)
fgseaKegg <- fgsea(pathways = kegg_list, 
                  stats    = deframe(rownames_to_column(DEqMS.results["logFC"])),
                  minSize  = 15,
                  maxSize  = 500)

##### Plot
topPathwaysUp <- fgseaKegg[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaKegg[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(kegg_list[topPathways], deframe(rownames_to_column(DEqMS.results["logFC"])), fgseaKegg, 
              gseaParam=0.5)

barplot_DEP_kegg <- dplyr::filter(fgseaKegg, pathway %in% topPathways) %>%
  arrange(NES) %>% 
  mutate(pathway = str_remove(pathway,"KEGG_"),
         pathway = as_factor(pathway)) %>% 
  #dplyr::filter(p.adjust < p.adjust_th) %>%
  #slice_tail(n=20) %>%
  ggplot(aes(x = NES, y = pathway, fill = padj)) +
  geom_bar(stat="identity") +
  scale_fill_gradientn(colours = c("red", "white", "grey"), guide=guide_colourbar(reverse = TRUE),
                       rescaler = ~ scales::rescale_mid(., mid = 0.05))+
  theme_bw() +
  labs(title = "KEGG") +
  ylab("")

print(barplot_DEP_kegg)
# ggsave(barplot_DEP_kegg,
#        filename = "results/Figures/barplot_DEP_kegg.svg",
#        width = 8,
#        height = 5)
#-------------------------------------------------------------------------------
# Survival curves regarding ISR status
surv_data <- dplyr::rename(CPTAC_ISRact_projection, case_id = sample) %>%
  inner_join(clinical_data, by = "case_id")

## Kaplan Meyer 33up vs 33down
fit <-  survfit(Surv(follow_up_days, status) ~ ISRactPCA_bin, 
                data = dplyr::filter(surv_data, ISRactPCA_bin != "intermediate_ISRactPCA"))
print(fit)

### Change color, linetype by strata, risk.table color by strata
kmos_cptac_ISRact <- ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("tomato3", "seagreen"))

print(kmos_cptac_ISRact)
# ggsave(paste0("results/survplots/FIGURES/kmos_cptac_ISRact.svg"),
#        kmos_cptac_ISRact$plot, height = 4 , width = 4)

## Cox Proportional hazzards model
cox.mod <- coxph(Surv(follow_up_days, status) ~ ISRactPCA, 
                 data = as.data.frame(surv_data))
### assumption checking
#### Linearity
plot(predict(cox.mod),
     residuals(cox.mod, type = "martingale"),
     xlab = "fitted", ylab = "Martingale residuals",
     main = "Residuals Plot", las = 1)
abline(h=0)
lines(smooth.spline(predict(cox.mod),
                    residuals(cox.mod, type = "martingale")),
                    col="red")

#### Proportional Hazzards
par(mfrow = c(1,1))
plot(cox.zph(cox.mod)) # if failed, Proportional
abline(h=0, col =2)

### Model plotting
ggforest(cox.mod)
ggadjustedcurves(cox.mod, data=as.data.frame(surv_data))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Survival curves regarding PHGDH/CBS status
surv_data <- as_tibble(CPTAC_tumor_toproj, rownames = "case_id") %>% 
  dplyr::select(case_id, reverse_translate[c("CBS", "PHGDH")]) %>%
  inner_join(clinical_data, by = "case_id")

measure = "mean" # overlap | mean
## Calculate the measure to split samples
if(measure == "mean"){
### min mean
surv_data <- dplyr::mutate(surv_data, 
                           meanPHGDHCBS = rowMeans(surv_data[c("PHGDH", "CBS")]),
                           stratPHGDHCBS = if_else(meanPHGDHCBS < quantile(meanPHGDHCBS, probs = 0.3333), "low_PHGDHCBS",
                                                   if_else(meanPHGDHCBS < quantile(meanPHGDHCBS, probs = 0.6666), "indetermined", "high_PHGDHCBS")))
} else if(measure == "overlap"){
### overlap min and overlap max
surv_data <- dplyr::mutate(surv_data,
                           ismaxPHGDH = if_else(PHGDH > median(PHGDH), "high_PHGDHCBS", "low_PHGDHCBS"),
                           ismaxCBS = if_else(CBS > median(CBS), "high_PHGDHCBS", "low_PHGDHCBS"),
                           stratPHGDHCBS = if_else(ismaxPHGDH == ismaxCBS, ismaxPHGDH, "indetermined")
)
}
## Visualize stratification
strat_visu <- ggplot(surv_data, aes(PHGDH, CBS, color = stratPHGDHCBS)) +
  scale_color_manual(values = c(high_PHGDHCBS = "tomato3", low_PHGDHCBS = "seagreen", indetermined = "darkgrey")) +
  geom_point() + theme_classic()

# ggsave(paste0("results/survplots/FIGURES/cptac_split", measure ,"_phgdhcbs.svg"),
#        strat_visu, height = 3 , width = 4)
## Kaplan Meyer of the selected measure
fit <-  survfit(Surv(follow_up_days, status) ~ stratPHGDHCBS, 
                data = dplyr::filter(surv_data, !stratPHGDHCBS %in% c("medium_PHGDHCBS", "indetermined")))
print(fit)

### Change color, linetype by strata, risk.table color by strata
kmos_cptac <- ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("tomato3", "seagreen")) +
  labs(title = paste0(measure, " split of PHGDH | CBS expression"))

print(kmos_cptac)
# ggsave(paste0("results/survplots/FIGURES/kmos_cptac_", measure ,"_phgdhcbs.svg"),
#        kmos_cptac$plot, height = 4 , width = 4)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Deconvolution proportions
df_tocorr <- inner_join(Molecular_phenotype_data, CPTAC_ISRact_projection, "sample") %>%
  rename_all( ~ str_replace_all(., " ", ".")) %>% 
  rename_all( ~ str_replace_all(., "-", "."))

corrplot_CPTACmeta <- df_tocorr %>% 
  dplyr::select(ISRactPCA, matches("MCPCounter"), matches("pathway"), matches("deconv")) %>% # , matches("deconv") misses samples!
  formatted_cors("spearman") %>%
  filter(measure1 == "ISRactPCA",
         measure2 != "ISRactPCA") %>% #not square corr
  dplyr::mutate(measure2 = fct(measure2, levels = unique(measure2))) %>% 
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nCorrelation", title="Correlations ISRactPCA vs CPTAC variables",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  ggpubr::rotate_x_text(angle = 0)

# ggsave(file="results/Figures/corrplot_CPTACmeta.svg", plot=corrplot_CPTACmeta, width=4, height=6)



tcell_ISRactPCA <- correlation_plotter(df_tocorr, "ISRactPCA", "T.cells_MCPCounter", "")
# ggsave(tcell_ISRactPCA,
#        filename = "results/Figures/tcell_ISRactPCA.svg",
#        width = 4,
#        height = 4)

cd8tcell_ISRactPCA <-correlation_plotter(df_tocorr, "ISRactPCA", "CD8.T.cells_MCPCounter", "")
# ggsave(cd8tcell_ISRactPCA,
#        filename = "results/Figures/cd8tcell_ISRactPCA.svg",
#        width = 4,
#        height = 4)

cytotcell_ISRactPCA <-correlation_plotter(df_tocorr, "ISRactPCA", "Cytotoxic.lymphocytes_MCPCounter", "")
# ggsave(cytotcell_ISRactPCA,
#        filename = "results/Figures/cytotcell_ISRactPCA.svg",
#        width = 4,
#        height = 4)

Bcell_ISRactPCA <-correlation_plotter(df_tocorr, "ISRactPCA", "B.lineage_MCPCounter", "")
# ggsave(Bcell_ISRactPCA,
#        filename = "results/Figures/Bcell_ISRactPCA.svg",
#        width = 4,
#        height = 4)

jakstat_ISRactPCA <-correlation_plotter(df_tocorr, "ISRactPCA", "JAK.STAT_pathway", "")
# ggsave(jakstat_ISRactPCA,
#        filename = "results/Figures/jakstat_ISRactPCA.svg",
#        width = 4,
#        height = 4)

stromal_deconv <-correlation_plotter(df_tocorr, "ISRactPCA", "stromal_deconv", "")
# ggsave(jakstat_ISRactPCA,
#        filename = "results/Figures/stromal_deconv.svg",
#        width = 4,
#        height = 4)

egfr_pathway <-correlation_plotter(df_tocorr, "ISRactPCA", "EGFR_pathway", "")
# ggsave(egfr_pathway,
#        filename = "results/Figures/egfr_pathway.svg",
#        width = 4,
#        height = 4)
#-------------------------------------------------------------------------------