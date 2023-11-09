#Import library
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

################################################################################
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")
source("src/formatted_cors.R")

## CPTAC (Proteogenomics)
### upper quartile normalized and log2 data 
CPTAC_tumor_already_norm <- read_delim("data/PDAC_LinkedOmics_Data/mRNA_RSEM_UQ_log2_Tumor.cct", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE) %>%
  dplyr::rename(Gene = ...1)

## Raw counts obtained by email
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

### inherited sample info and top samples from Sauyeun_PDX
sample_info <- read_delim("data/Sauyeun_PDX/sample_info.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

top_samples <-arrange(sample_info, ICA3) %>%
  dplyr::slice( unique(c(1:5, n() - 0:4)) ) %>%
  mutate(ISRact = ifelse(ICA3 < 0, "low_ISRact", "high_ISRact")) %>%
  arrange(sample)

# Import the projection
## ISRactPCA
pca_pdx <- read_rds("data/Classifiers/pca_pdx_ENZO.RDS")

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
#Version 95 for Proteogenimics data
ensembl95 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 95)#listAttributes(ensembl95, page="feature_page")

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

## Projecting datasets on this PCA
### CPTAC
#### transpose human data for projection
human_data <- CPTAC_tumor %>%
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
pca_pdx$rotation <- filter_pca(pca_pdx$rotation, names(human_data)[-1])
pca_pdx$center <- filter_pca(pca_pdx$center, names(human_data)[-1])
pca_pdx$scale <- filter_pca(pca_pdx$scale, names(human_data)[-1])

projection_CPTAC <-  predict(pca_pdx, human_data) %>% 
  as_tibble(rownames = "sample")


CPTAC_PC1 <- arrange(projection_CPTAC, PC1) %>%
  mutate(PC1status = if_else(PC1 < quantile(projection_CPTAC$PC1, probs = 0.3333), "low_ISRact",
                             if_else(PC1 < quantile(projection_CPTAC$PC1, probs = 0.6666), "medium_ISRact", "high_ISRact"))) %>%
  dplyr::select(sample, PC1,  PC1status) %>%
  left_join(sample_info_CPTAC, by = "sample")

#Plot pca and projections
#Add ISR status to PCA df
pca_full_df <- pca_pdx[["x"]] %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  inner_join(top_samples[,c("sample","ISRact")])

#-------------------------------------------------------------------------------
# Plot projecion PacaOmics and then projection Proteogenomics

show_projection <- bind_rows(as_tibble(pca_full_df) %>% mutate(dataset = "Sauyeun PDX"),
                             as_tibble(projection_CPTAC) %>% mutate(ISRact = "Unknown",
                                                                    dataset = "CPTAC"))

dplyr::mutate(show_projection, ISRact = str_replace(ISRact, 'ICA3', 'ISRact')) %>%
  dplyr::filter(dataset %in% c("Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"), # "Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"
                ISRact %in% c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) %>% 
  ggplot(aes(x=PC1, y=PC2, color = ISRact, shape = dataset)) +
  geom_point() +
  scale_shape_discrete(limits = c("Sauyeun PDX", "PACAOMICS PDX", "CPTAC", "CCLE")) +
  scale_color_discrete(limits = c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) + #c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown'))unique(projection_ccle$primary_tissue))
  ylim(-250,15)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Check relation with original PDX samples
PDX_original <- read_tsv("data/Sauyeun_PDX/PDX_PCA_input_data.tsv")

## Genome wide
similarity_df <- inner_join(PDX_original ,dplyr::select(CPTAC_tumor, - Gene)) %>% 
  drop_na()
similarity_corr_hm <- formatted_cors(dplyr::select(similarity_df,-EnsemblID), "pearson", 0) %>% 
  pivot_wider(dplyr::select(similarity_corr, measure1, measure2, r), values_from = r, names_from = measure1, id_cols = measure2)

### Annotation

annot_ref = dplyr::select(top_samples, sample, ISRact, PAMG) %>% dplyr::mutate(df = "Sauyeun")
annot_proj = dplyr::rename(CPTAC_PC1,  ISRact = PC1status) %>% dplyr::select(sample, ISRact, PAMG) %>% dplyr::mutate(df = "CPTAC")
annot = bind_rows(annot_ref, annot_proj) %>% column_to_rownames("sample")

annot_colors <- list(ISRact = c(`high_ISRact` = "brown", `medium_ISRact` = "grey", `low_ISRact` = "#006837"),
                     PAMG = c("#FF7F00", "white", "#377DB8"),
                     df = c(Sauyeun = "grey", CPTAC ="grey2")) 

### Plot
heatmap <- similarity_corr_hm %>% column_to_rownames("measure2") %>% 
  pheatmap(scale = "row", color = colorRampPalette(c("#FBFEF9", "#A63446"))(100),
           annotation_row = annot, annotation_col = annot, annotation_colors = annot_colors,
           cluster_rows = T, cluster_cols = T, show_colnames = TRUE) 

## Subset

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Plot comparisons with Basal/Classical and ISRact
## PAMG vs PC1
correlation_plotter(data = CPTAC_PC1, col1 = "PAMG", col2 = "PC1", data_name = "CPTAC")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Proteomics DPEA
## load already normalized log transformed and centered data
CPTAC_comparison_metadata <- dplyr::filter(CPTAC_PC1, PC1status != "medium_ISRact")
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
cond <- fct(CPTAC_comparison_metadata$PC1status)

design = model.matrix(~0+cond)
colnames(design) = gsub("cond","",colnames(design))

x <- c("low_ISRact-high_ISRact")
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
DEqMS.results$log.sca.pval = -log10(DEqMS.results$sca.P.Value)
ggplot(DEqMS.results, aes(x = logFC, y =log.sca.pval )) + 
  geom_point(size=0.5 )+
  theme_bw(base_size = 16) + # change theme
  xlab(expression("log2(ISRactPCA high/low)")) + # x-axis label
  ylab(expression(" -log10(P-value)")) + # y-axis label
  geom_vline(xintercept = c(-1,1), colour = "red") + # Add fold change cutoffs
  geom_hline(yintercept = 3, colour = "red") + # Add significance cutoffs
  geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
  scale_colour_gradient(low = "black", high = "black", guide = FALSE)+
  geom_text_repel(data=subset(DEqMS.results, abs(logFC)>0.5&log.sca.pval > 3),
                  aes( logFC, log.sca.pval ,label=gene)) # add gene label

### same volcano but highlighting some genes
highlight = c("PHGDH","IL18", "IDO1",
              "APP")

ggplot(DEqMS.results, aes(x = logFC, y =log.sca.pval )) + 
  geom_point(size=0.5 )+
  theme_bw(base_size = 16) + # change theme
  xlab(expression("log2(ISRactPCA high/low)")) + # x-axis label
  ylab(expression(" -log10(P-value)")) + # y-axis label
  geom_vline(xintercept = c(-1,1), colour = "red") + # Add fold change cutoffs
  geom_hline(yintercept = 3, colour = "red") + # Add significance cutoffs
  geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
  scale_colour_gradient(low = "black", high = "black", guide = FALSE) +
  geom_point(data= DEqMS.results[rownames(DEqMS.results) %in% highlight,],
             aes( logFC, log.sca.pval), color="red") + 
  geom_label_repel(data= DEqMS.results[rownames(DEqMS.results) %in% highlight,],
                  aes( logFC, log.sca.pval ,label=gene)) # add gene label


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

#-------------------------------------------------------------------------------
# Survival curves reguarding ISR status
surv_data <- dplyr::rename(CPTAC_PC1, case_id = sample) %>%
  inner_join(clinical_data, by = "case_id")

## Kaplan Meyer 33up vs 33down
fit <-  survfit(Surv(follow_up_days, status) ~ PC1status, 
                data = dplyr::filter(surv_data, PC1status != "medium_ISRact"))
print(fit)

### Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Cox Proportional hazzards model
cox.mod <- coxph(Surv(follow_up_days, status) ~ PC1, 
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
