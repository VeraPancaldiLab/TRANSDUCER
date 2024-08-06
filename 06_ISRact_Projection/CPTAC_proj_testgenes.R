# DGEA on CPTAC human cohort with respect to ISRactPCA

################################################################################
#' Import rawcount
#' follow the limma voom workflow for DGEA analysis
#' GSEA on ranking regarding DGEA results
#' Comparative analysis of most correlated TF activities
################################################################################

# Import libraries
library(readxl)
library(edgeR)
library(ggpubr)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(EnhancedVolcano)
library(fgsea)
library(msigdbr)
library(RColorBrewer)
library(dorothea)
library(pheatmap)

## setwd and import local functions
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/plot_PCs.R")
source("src/plotBestPathways.R")



# Import data
## Expression
counts <- read_delim("data/PDAC_LinkedOmics_Data/CPTAC-PDAC-raw-RSEM-expected-counts-gene-level.txt",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
) %>%
  column_to_rownames("idx") %>%
  dplyr::select(contains("tumor")) %>%
  rename_all(funs(str_replace_all(., "_tumor", "")))

## Metadata
clinical_data <- read_excel("data/PDAC_LinkedOmics_Data/mmc1.xlsx",
  sheet = "Clinical_data"
) %>%
  mutate(follow_up_days = as.numeric(follow_up_days)) %>%
  mutate(status = ifelse(vital_status == "Deceased", 2, 1))

Molecular_phenotype_data <- read_excel("data/PDAC_LinkedOmics_Data/mmc1.xlsx",
  sheet = "Molecular_phenotype_data"
) %>%
  mutate_at(vars(immune_deconv:`necrosis_(%OF_TUMOR_WITH_NECROSIS)_histology_estimate`, KRAS_VAF), as.numeric) %>%
  mutate_at(vars(Bailey:Moffitt), as.factor) # change type to avoid errors

## Import Projection
human_ISRactPCA <- read_csv("results/CPTAC/human_PC1.csv") %>% # This comes from 02_TRANSDUCER/06_Enzos_Internship/human_cohort_analysis/Data
  dplyr::rename(ISRactPCA = "PC1", ISRactPCA_bin = "PC1status") %>%
  dplyr::mutate(ISRactPCA_bin = str_replace(ISRactPCA_bin, "PC1", "ISRactPCA"))
## Unify
sample_info <- dplyr::rename(human_ISRactPCA, case_id = sample) %>%
  inner_join(clinical_data, by = "case_id") %>%
  inner_join(Molecular_phenotype_data, by = "case_id") %>%
  mutate(Diabetes = ifelse(str_detect(medical_condition, "Diabetes"), "yes", "no"))


# DGEA following the Limma Voom pipeline from https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
################################################################################

## Processing and object building
### Visualize expected counts by taking 10 random samples
random10 <- sort(sample.int(dim(sample_info)[1], 10))

data <- rownames_to_column(counts[, random10], "Prot") %>%
  pivot_longer(cols = 2:11, names_to = "case_id", values_to = "Expected_count") # pivot_longer data to make them easily plotable with ggplot 2

ggplot(data, mapping = aes(x = case_id, y = Expected_count, fill = case_id)) +
  geom_violin() +
  theme(legend.position = "none") +
  coord_flip() +
  yscale("log10")

### Create DGEList object
d0 <- dplyr::select(counts, sample_info$case_id) %>%
  DGEList()

### Calculate normalization factors
d0 <- calcNormFactors(d0, method = "TMM")
d0

### Filter low-expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop, ]
dim(d) # number of genes left

### Create a new variable “group”
group <- as.factor(human_ISRactPCA$ISRactPCA_bin)


### Specify the model to be fitted
mm <- model.matrix(~ 0 + group)

### Filter by gene expression
keep <- filterByExpr(d0, mm)
d <- d0[keep, ]
dim(d) # number of genes left

## Visualize before DGEA
dnorm <- cpm(d)

data2 <- as.data.frame(dnorm[, random10]) %>%
  rownames_to_column("Prot") %>%
  pivot_longer(cols = 2:11, names_to = "case_id", values_to = "Expected_count_norm") # pivot_longer data to make them easily plotable with ggplot 2

ggplot(data2, mapping = aes(x = case_id, y = Expected_count_norm, fill = case_id)) +
  geom_violin() +
  theme(legend.position = "none") +
  coord_flip() +
  yscale("log10")

### Multidimensional scaling (MDS) plot
plotMDS(d, col = as.numeric(group))

### PCA plot regarding clinical and phenotype data
dnorm.t <- t(dnorm)
dnorm.pca <- PCA(dnorm.t, ncp = 10, graph = FALSE)
summary(dnorm.pca)

fviz_eig(dnorm.pca, addlabels = TRUE, ylim = c(0, 50))
eig.val <- get_eigenvalue(dnorm.pca)
eig.val # 63 dimensions to reach 90% of explained variance

# Join PCA results and more relevant clinical data
pca_fulldf <- dnorm.pca$ind$coord %>%
  as_tibble(rownames = "case_id") %>%
  inner_join(Molecular_phenotype_data[, c(1:14, 15, 16, 17, 21)]) %>%
  inner_join(clinical_data[, c(1, 4, 5, 6, 8, 9, 11, 12, 13, 14, 22)])

for (feat in colnames(pca_fulldf)[-c(1:11)]) {
  plot_PCs(pca_fulldf, feat, 10, 1)
}

## Voom transformation and calculation of variance weights
y <- voom(d, mm, plot = T)

## Fitting linear models in limma
fit <- lmFit(y, mm) # lmFit fits a linear model using weighted least squares for each gene
head(coef(fit))

### define contrast matrix
contr <- makeContrasts(grouphigh_ISRactPCA - grouplow_ISRactPCA, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr) # Estimate contrast for each gene

tmp <- eBayes(tmp) # Empirical Bayes smoothing of standard errors

## Visualize top differentiated genes with vulcano plot
top_diff_gene <- topTreat(tmp, coef = 1, n = Inf) %>%
  rownames_to_column("external_gene_name")
EnhancedVolcano(top_diff_gene,
  lab = top_diff_gene$external_gene_name, x = "logFC", y = "P.Value",
  # selectLab = c("PHGDH","CBS","IL18","IFNG"),
  pCutoffCol = "adj.P.Val", pCutoff = 0.05, title = "Low ISRactPCA vs. high ISRactPCA",
  subtitle = "adj.P.Val < 0.05"
)

### other visualisations of the viggnettes
top.table <- topTable(tmp, sort.by = "P", n = Inf) %>%
  mutate(
    Expression = case_when(
      logFC >= log(2) & adj.P.Val <= 0.05 ~ "Up-regulated",
      logFC <= -log(2) & adj.P.Val <= 0.05 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    )
  ) %>%
  rownames_to_column("Gene")

head(top.table, 20)

dt <- decideTests(fit)
summary(dt)
vennDiagram(dt[, 1:2], circle.col = c("turquoise", "salmon"))

# GSEA on l2FC
## load gene sets
msigdbr_reactome_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  mutate(gs_name = str_remove(gs_name, "REACTOME_"))
msigdbr_hallmark_df <- msigdbr(species = "Homo sapiens", category = "H") %>%
  mutate(gs_name = str_remove(gs_name, "HALLMARK_"))

## extract ranked list
ranked_genes <- top_diff_gene %>%
  as_tibble() %>%
  arrange(desc(logFC))
ranked_genes <- ranked_genes[, c(1, 2)] %>% deframe()

msigdbr_list_reactome <- split(x = msigdbr_reactome_df$gene_symbol, f = msigdbr_reactome_df$gs_name)
msigdbr_list_hallmark <- split(x = msigdbr_hallmark_df$gene_symbol, f = msigdbr_hallmark_df$gs_name)

## analysis
gsea_res_reactome <- fgsea(pathways = msigdbr_list_reactome, stats = ranked_genes, minSize = 15)
gsea_res_hallmark <- fgsea(pathways = msigdbr_list_hallmark, stats = ranked_genes, minSize = 15)

## Plots
### Official viggnette plots
plotBestPathways(gsea_res_reactome, msigdbr_list_reactome, ranked_genes)
plotBestPathways(gsea_res_hallmark, msigdbr_list_hallmark, ranked_genes)

### custom plots thesis
#### Reactome
topPathwaysUp <- gsea_res_reactome[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown <- gsea_res_reactome[ES < 0][head(order(pval), n = 10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

cptacbarplot_DEGEA_reactome <- dplyr::filter(gsea_res_reactome, pathway %in% topPathways) %>%
  arrange(NES) %>%
  mutate(
    pathway = str_remove(pathway, "REACTOME_"),
    pathway = as_factor(pathway)
  ) %>%
  # dplyr::filter(p.adjust < p.adjust_th) %>%
  # slice_tail(n=20) %>%
  ggplot(aes(x = NES, y = pathway, fill = padj)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(
    colours = c("red", "white", "grey"), guide = guide_colourbar(reverse = TRUE),
    rescaler = ~ scales::rescale_mid(., mid = 0.05)
  ) +
  theme_bw() +
  labs(title = "REACTOME") +
  ylab("")

print(cptacbarplot_DEGEA_reactome)
# ggsave(cptacbarplot_DEGEA_reactome,
#        filename = "results/Figures/cptacbarplot_DEGEA_reactome.svg",
#        width = 8,
#        height = 5)

#### Hallmark
topPathwaysUp <- gsea_res_hallmark[ES > 0][head(order(pval), n = 10), pathway]
topPathwaysDown <- gsea_res_hallmark[ES < 0][head(order(pval), n = 10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

cptacbarplot_DEGEA_hallmark <- dplyr::filter(gsea_res_hallmark, pathway %in% topPathways) %>%
  arrange(NES) %>%
  mutate(
    pathway = as_factor(pathway),
    padj = ifelse(padj < 0.1, padj, NA)
  ) %>% # else the colour bar is screwed
  # dplyr::filter(p.adjust < p.adjust_th) %>%
  # slice_tail(n=20) %>%
  ggplot(aes(x = NES, y = pathway, fill = padj)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(
    colours = c("red", "white"), guide = guide_colourbar(reverse = TRUE),
    rescaler = ~ scales::rescale_mid(., mid = 0.05)
  ) +
  theme_bw() +
  labs(title = "HALLMARK") +
  ylab("")

print(cptacbarplot_DEGEA_hallmark)
# ggsave(cptacbarplot_DEGEA_hallmark,
#        filename = "results/Figures/cptacbarplot_DEGEA_hallmark.svg",
#        width = 8,
#        height = 5)

# most correlated TF activities
## Generate TF act
### Select the regulon data
data(dorothea_hs_pancancer, package = "dorothea")
regulons <- dorothea_hs_pancancer %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))

### VIPER Analysis
################################################################################
n <- 30 # n TFs to plot in heatmaps
minsize <- 5
ges.filter <- FALSE
################################################################################
tf_act__ <- dorothea::run_viper(counts, regulons,
  options = list(
    minsize = minsize, eset.filter = ges.filter,
    cores = 1, verbose = FALSE, nes = TRUE
  )
)

# Format for analysis
tf_act_ <- as_tibble(tf_act__, rownames = "TF") %>%
  pivot_longer(cols = -TF, names_to = "bcr_patient_barcode") %>%
  pivot_wider(id_cols = c(bcr_patient_barcode), names_from = TF)

### Visualisation
#### Most variable TFs
##### Calculate the n most variable TFs
tf_50_var <- summarise(tf_act_, across(where(is.double), var)) %>%
  pivot_longer(cols = everything(), names_to = "TF", values_to = "var") %>%
  dplyr::arrange(desc(var)) %>%
  dplyr::slice_head(n = n)

###### Define Annotations
annotation_row_var <- tf_50_var %>%
  column_to_rownames(var = "TF")

annotation_col <- dplyr::select(sample_info, case_id, ISRactPCA, KRAS_VAF, Moffitt, Bailey, Diabetes) %>%
  column_to_rownames(var = "case_id")

ann_colors <- list(
  ISRactPCA = c("seagreen", "white", "tomato3"),
  ISRactPCA_weight = c("white", "darkorange"),
  KRAS_VAF = c("white", "#0B2356"),
  Moffitt = c(`BASAL-LIKE` = "#FF7F00", CLASSICAL = "#0000FC", `NA` = "grey"),
  Bailey = c(ADEX = "#713E23", immunogenic = "#EF2C18", `NA` = "darkgrey", `pancreatic progenitor` = "#0000FC", squamous = "#FF7F00"),
  Diabetes = c(no = "grey", yes = "black"),
  var = c("white", "darkgreen"),
  cor = brewer.pal(n = 11, name = "BrBG")
)

# Create breaks and limits for the continuous expression scale
breaksList <- seq(-3.5, 3.5, by = 0.1)

##### Heatmap
as_tibble(tf_act__, rownames = "TF") %>%
  dplyr::select(TF, human_ISRactPCA$sample) %>%
  dplyr::filter(TF %in% tf_50_var$TF) %>%
  column_to_rownames("TF") %>%
  as.matrix() %>%
  pheatmap(
    annotation_col = annotation_col,
    annotation_row = annotation_row_var,
    annotation_colors = ann_colors, scale = "row", show_colnames = FALSE,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
    breaks = breaksList,
    fontsize = 12,
    annotation_names_row = FALSE,
    main = paste0("Activity of the ", n, " most variable transcription factor on the Sauyeun PDX dataset")
  )

#### Most correlated with ISRact
##### calculate correlation TF ISRact
tf_ISRactPCA_cor <- tf_act_ %>%
  dplyr::rename(sample = bcr_patient_barcode) %>%
  inner_join(dplyr::select(human_ISRactPCA, sample, ISRactPCA), by = "sample") %>% # make the two table in same order
  arrange(ISRactPCA) %>%
  dplyr::select(-ISRactPCA) %>%
  column_to_rownames("sample") %>% # select continuous variables
  as.matrix() %>%
  cor(y = human_ISRactPCA$ISRactPCA, method = "spearman") %>%
  as.data.frame()

##### filter by absolute correlation
top_tf <- mutate(tf_ISRactPCA_cor, V1 = abs(V1)) %>%
  arrange(desc(V1)) %>%
  dplyr::slice(1:n) %>%
  dplyr::rename(cor = V1) %>%
  rownames_to_column(var = "TF")

##### Define annotations
annotation_row_cor <- tf_ISRactPCA_cor %>%
  rownames_to_column("TF") %>%
  inner_join(dplyr::select(top_tf, TF), by = "TF") %>%
  dplyr::rename(cor = V1) %>%
  column_to_rownames(var = "TF")

##### Heatmap
as_tibble(tf_act__, rownames = "TF") %>%
  dplyr::select(TF, human_ISRactPCA$sample) %>%
  dplyr::filter(TF %in% top_tf$TF) %>%
  column_to_rownames("TF") %>%
  as.matrix() %>%
  pheatmap(
    annotation_col = annotation_col,
    annotation_row = annotation_row_cor,
    annotation_colors = ann_colors, scale = "row", show_colnames = FALSE,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
    breaks = breaksList,
    fontsize = 12,
    annotation_names_row = FALSE,
    main = "Activity of the 50 most correlated transcription factor  with ISRactPCA on the CPTAC-PDAC dataset",
    height = 6,
    width = 10,
    filename = "results/Figures/cptac_ISRactPCA_TFact.pdf"
  )
